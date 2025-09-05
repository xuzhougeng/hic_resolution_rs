use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::stdin;
use std::path::PathBuf;

use crate::{coverage, parser, resolution, straw, utils};
use rayon::prelude::*;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(name = "hic_resolution")]
#[command(about = "Fast calculation of Hi-C map resolution")]
pub struct Cli {
    /// Path to merged_nodups file (can be .gz compressed)
    #[arg(value_name = "MERGED_NODUPS")]
    pub nodups: Option<PathBuf>,

    /// Path to chromosome sizes file
    #[arg(short, long, value_name = "CHROM_SIZE")]
    pub chrom_size: Option<PathBuf>,

    /// Total genome size in base pairs
    #[arg(long, default_value = "1000000000")]
    pub genome_size: u64,

    /// Minimum bin size (base pairs)
    #[arg(long, default_value = "50")]
    pub bin_width: u32,

    /// Proportion of bins that must meet coverage threshold
    #[arg(long, default_value = "0.8")]
    pub prop: f64,

    /// Minimum contacts per bin to be considered "good"
    #[arg(long, default_value = "1000")]
    pub count_threshold: u32,

    /// Step size for initial coarse search
    #[arg(long, default_value = "1000")]
    pub step_size: u32,

    /// Number of threads to use (0 = auto)
    #[arg(short, long, default_value = "4")]
    pub threads: usize,

    /// Aggregation chunk size in number of pairs
    /// Default chosen to fit comfortably under ~8 GB RAM
    #[arg(long, value_name = "PAIRS", default_value = "4000000")]
    pub chunk_pairs: usize,

    /// Per-worker subchunk size in number of pairs
    /// Default tuned for throughput under ~8 GB RAM
    #[arg(long, value_name = "PAIRS", default_value = "128000")]
    pub subchunk_pairs: usize,

    /// Optional subcommand. Use `straw` to work with .hic slices.
    #[command(subcommand)]
    pub cmd: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Straw-compatible utilities
    Straw(StrawCli),
}

#[derive(Args, Debug)]
pub struct StrawCli {
    #[command(subcommand)]
    pub cmd: StrawCmd,
}

#[derive(Subcommand, Debug)]
pub enum StrawCmd {
    /// Dump genome-wide observed counts at resolution to a slice file (.slc.gz)
    Dump {
        /// observed/oe/expected (only observed supported)
        matrix_type: String,
        /// NONE/VC/VC_SQRT/KR (only NONE supported)
        norm: String,
        /// Input Hi-C file (.hic)
        input: PathBuf,
        /// Units (BP only supported)
        unit: String,
        /// Bin size / resolution in bp
        binsize: i32,
        /// Output file path (.slc.gz)
        output: PathBuf,
    },
    /// List chromosomes in a .hic file
    List {
        /// Input Hi-C file (.hic)
        input: PathBuf,
    },
    /// Estimate effective resolution / coverage
    Effres {
        /// Input Hi-C file (.hic)
        input: PathBuf,
        /// Chromosome name, e.g. 1 / chr1 / X. Omit to summarize across all chromosomes.
        chromosome: Option<String>,
        /// Minimum contacts per bin to count as covered
        #[arg(long, default_value_t = 1000)]
        thr: i32,
        /// Coverage fraction threshold (0–1)
        #[arg(long, default_value_t = 0.8)]
        pct: f64,
    },
}

pub fn run() -> Result<()> {
    let args = Cli::parse();

    // Subcommands take precedence
    if let Some(Commands::Straw(cli)) = &args.cmd {
        return run_straw(cli);
    }

    // Set thread pool size
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    println!("Hi-C Resolution Calculator (Rust)");
    println!("=================================");

    // Create coverage structure (auto-detect pairtools header if present)
    let chrom_size_path = args.chrom_size.as_ref().map(|p| p.to_str().unwrap());
    let mut pairs_mode = false;
    let mut pairs_chr_map: Option<utils::ChrLookup> = None;
    let genome_names: Vec<String>;
    let genome_lengths: Vec<u32>;

    // Decide source of chromosome names + lengths, and build coverage
    let mut coverage = if let Some(path) = args.nodups.as_ref() {
        if let Ok(Some((map, names, lengths))) = parser::sniff_pairs_header_from_path(path.as_path()) {
            pairs_mode = true;
            pairs_chr_map = Some(map);
            genome_names = names;
            genome_lengths = lengths.clone();
            coverage::Coverage::from_lengths(args.bin_width, lengths)
        } else {
            if let Some(cs) = chrom_size_path {
                let (names, lengths) = utils::read_chrom_sizes_with_names(cs)?;
                genome_names = names;
                genome_lengths = lengths.clone();
                coverage::Coverage::from_lengths(args.bin_width, lengths)
            } else {
                genome_names = utils::get_default_genome_names();
                genome_lengths = utils::get_default_genome_lengths();
                coverage::Coverage::from_lengths(args.bin_width, genome_lengths.clone())
            }
        }
    } else {
        if let Some(cs) = chrom_size_path {
            let (names, lengths) = utils::read_chrom_sizes_with_names(cs)?;
            genome_names = names;
            genome_lengths = lengths.clone();
            coverage::Coverage::from_lengths(args.bin_width, lengths)
        } else {
            genome_names = utils::get_default_genome_names();
            genome_lengths = utils::get_default_genome_lengths();
            coverage::Coverage::from_lengths(args.bin_width, genome_lengths.clone())
        }
    };
    // Now that we have names + lengths, print computed genome info and settings
    let genome_size: u64 = genome_lengths.iter().map(|&x| x as u64).sum();
    println!("Genome size: {} bp", genome_size);
    println!("Bin width: {} bp", args.bin_width);
    println!("Coverage threshold: {} contacts", args.count_threshold);
    println!("Required proportion: {:.1}%", args.prop * 100.0);
    println!("Chromosome lookup: {}", utils::chr_lookup_impl());
    // Top-10 chromosomes by length (descending)
    if !genome_names.is_empty() && !genome_lengths.is_empty() {
        let mut pairs: Vec<(&str, u32)> = genome_names
            .iter()
            .map(|s| s.as_str())
            .zip(genome_lengths.iter().copied())
            .collect();
        pairs.sort_unstable_by(|a, b| b.1.cmp(&a.1));
        let topn = pairs.iter().take(10).collect::<Vec<_>>();
        println!("Top 10 chromosomes by length:");
        for (i, (nm, ln)) in topn.into_iter().enumerate() {
            println!("  {}. {}: {} bp", i + 1, nm, ln);
        }
    }
    println!();
    println!(
        "Initialized coverage tracking for {} chromosomes",
        coverage.bins.len()
    );

    // Set up progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")?
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ "),
    );

    // Parse input file and build coverage
    pb.set_message("Reading merged_nodups file...");
    let pairs_processed = if let Some(path) = args.nodups {
        let file = File::open(&path)?;
        let is_gz = path.extension().map_or(false, |ext| ext == "gz");
        if pairs_mode {
            let chr_map = pairs_chr_map.expect("pairs chr_map should be set");
            if is_gz {
                let iter = parser::open_pairs_file(file, chr_map)?;
                process_pairs(iter, &mut coverage, &pb, args.chunk_pairs, args.subchunk_pairs)?
            } else {
                let iter = parser::open_pairs_file_uncompressed(file, chr_map)?;
                process_pairs(iter, &mut coverage, &pb, args.chunk_pairs, args.subchunk_pairs)?
            }
        } else {
            if is_gz {
                let iter = parser::open_file(file, chrom_size_path)?;
                process_pairs(iter, &mut coverage, &pb, args.chunk_pairs, args.subchunk_pairs)?
            } else {
                let iter = parser::open_file_uncompressed(file, chrom_size_path)?;
                process_pairs(iter, &mut coverage, &pb, args.chunk_pairs, args.subchunk_pairs)?
            }
        }
    } else {
        // Read from stdin
        let iter = parser::open_file(stdin(), chrom_size_path)?;
        process_pairs(iter, &mut coverage, &pb, args.chunk_pairs, args.subchunk_pairs)?
    };

    pb.set_message("Computing resolution...");

    // Find resolution
    let resolution =
        resolution::find_resolution(&coverage, args.prop, args.count_threshold, args.step_size);

    pb.finish_and_clear();

    // Output results
    println!("Processed {} valid pairs", pairs_processed);
    println!();
    println!("Map resolution = {} bp", resolution);

    Ok(())
}

fn process_pairs<I>(
    iter: I,
    coverage: &mut coverage::Coverage,
    pb: &ProgressBar,
    chunk_pairs: usize,
    subchunk_pairs: usize,
) -> Result<u64>
where
    I: Iterator<Item = Result<utils::Pair>>,
{
    let mut count = 0u64;
    let mut buf: Vec<utils::Pair> = Vec::with_capacity(chunk_pairs.min(8_000_000));

    for pair_result in iter {
        let pair = pair_result?;
        buf.push(pair);
        if buf.len() >= chunk_pairs {
            aggregate_pairs_chunk(&buf, coverage, subchunk_pairs);
            buf.clear();
        }
        count += 1;

        if count % 1_000_000 == 0 {
            pb.set_message(format!(
                "Processed {:.1}M pairs...",
                count as f64 / 1_000_000.0
            ));
        }
    }

    if !buf.is_empty() {
        aggregate_pairs_chunk(&buf, coverage, subchunk_pairs);
        buf.clear();
    }

    Ok(count)
}

fn aggregate_pairs_chunk(pairs: &[utils::Pair], coverage: &mut coverage::Coverage, subchunk_pairs: usize) {
    let binw = coverage.bin_width;
    let chr_lens = &coverage.chr_lengths;

    // Process in parallel: for each subchunk, build a vector of (key, count),
    // where key packs (chrom_index, bin_index) into u64; then sort+compress.
    let scl = subchunk_pairs.max(16_000);
    let partials: Vec<Vec<(u64, u32)>> = pairs
        .par_chunks(scl)
        .map(|chunk| {
            #[inline]
            fn pack(ci: usize, b: u32) -> u64 { ((ci as u64) << 32) | (b as u64) }

            let mut vec: Vec<(u64, u32)> = Vec::with_capacity(chunk.len() * 2);
            for p in chunk {
                // First end
                let ci1 = (p.chr1 as usize).saturating_sub(1);
                if ci1 < chr_lens.len() {
                    let pos1 = p.pos1;
                    if pos1 < chr_lens[ci1] {
                        let b1 = pos1 / binw;
                        vec.push((pack(ci1, b1), 1));
                    }
                }
                // Second end
                let ci2 = (p.chr2 as usize).saturating_sub(1);
                if ci2 < chr_lens.len() {
                    let pos2 = p.pos2;
                    if pos2 < chr_lens[ci2] {
                        let b2 = pos2 / binw;
                        vec.push((pack(ci2, b2), 1));
                    }
                }
            }
            // sort by key and run-length compress counts
            vec.sort_unstable_by(|a, b| a.0.cmp(&b.0));
            let mut out: Vec<(u64, u32)> = Vec::with_capacity(vec.len());
            let mut it = vec.into_iter();
            if let Some((mut k, mut v)) = it.next() {
                for (kk, vv) in it {
                    if kk == k { v = v.saturating_add(vv); } else { out.push((k, v)); k = kk; v = vv; }
                }
                out.push((k, v));
            }
            out
        })
        .collect();

    // Merge compressed vectors into dense bins
    for part in partials {
        for (key, v) in part {
            let ci = (key >> 32) as usize;
            let b = (key & 0xFFFF_FFFF) as usize;
            if ci < coverage.bins.len() {
                let row = &mut coverage.bins[ci];
                if b < row.len() {
                    row[b] = row[b].saturating_add(v);
                }
            }
        }
    }
}

fn run_straw(cli: &StrawCli) -> Result<()> {
    match &cli.cmd {
        StrawCmd::Dump {
            matrix_type,
            norm,
            input,
            unit,
            binsize,
            output,
        } => {
            if matrix_type.to_ascii_lowercase() != "observed" {
                anyhow::bail!("Only 'observed' is supported in this Rust port");
            }
            if norm.to_ascii_uppercase() != "NONE" {
                anyhow::bail!("Only 'NONE' normalization is supported in this Rust port");
            }
            if unit.to_ascii_uppercase() != "BP" {
                anyhow::bail!("Only BP units are supported in this Rust port");
            }
            straw::dump_hic_genome_wide(input.as_path(), *binsize, output.as_path())
        }
        StrawCmd::List { input } => straw::list_hic_chromosomes(input.as_path()),
        StrawCmd::Effres {
            input,
            chromosome,
            thr,
            pct,
        } => straw::effres_hic(input.as_path(), chromosome.as_deref(), *thr, *pct),
    }
}
