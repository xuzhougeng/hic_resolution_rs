use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::stdin;
use std::path::PathBuf;

use crate::{coverage, parser, resolution, straw, utils};

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
    #[arg(long, default_value = "2428425688")]
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
    #[arg(short, long, default_value = "0")]
    pub threads: usize,

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
    /// Estimate effective resolution on a chromosome
    Effres {
        /// Input Hi-C file (.hic)
        input: PathBuf,
        /// Chromosome name, e.g. 1 / chr1 / X
        chromosome: String,
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
    println!("Genome size: {} bp", args.genome_size);
    println!("Bin width: {} bp", args.bin_width);
    println!("Coverage threshold: {} contacts", args.count_threshold);
    println!("Required proportion: {:.1}%", args.prop * 100.0);
    println!();

    // Create coverage structure (auto-detect pairtools header if present)
    let chrom_size_path = args.chrom_size.as_ref().map(|p| p.to_str().unwrap());
    let mut pairs_mode = false;
    let mut pairs_chr_map: Option<utils::ChrMap> = None;
    let coverage = if let Some(path) = args.nodups.as_ref() {
        if let Ok(Some((map, lengths))) = parser::sniff_pairs_header_from_path(path.as_path()) {
            pairs_mode = true;
            pairs_chr_map = Some(map);
            coverage::Coverage::from_lengths(args.bin_width, lengths)
        } else {
            coverage::Coverage::new(args.bin_width, chrom_size_path)
        }
    } else {
        coverage::Coverage::new(args.bin_width, chrom_size_path)
    };
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
                process_pairs(iter, &coverage, &pb)?
            } else {
                let iter = parser::open_pairs_file_uncompressed(file, chr_map)?;
                process_pairs(iter, &coverage, &pb)?
            }
        } else {
            if is_gz {
                let iter = parser::open_file(file, chrom_size_path)?;
                process_pairs(iter, &coverage, &pb)?
            } else {
                let iter = parser::open_file_uncompressed(file, chrom_size_path)?;
                process_pairs(iter, &coverage, &pb)?
            }
        }
    } else {
        // Read from stdin
        let iter = parser::open_file(stdin(), chrom_size_path)?;
        process_pairs(iter, &coverage, &pb)?
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

fn process_pairs<I>(iter: I, coverage: &coverage::Coverage, pb: &ProgressBar) -> Result<u64>
where
    I: Iterator<Item = Result<utils::Pair>>,
{
    let mut count = 0u64;

    for pair_result in iter {
        let pair = pair_result?;
        coverage.add_pair(&pair);
        count += 1;

        if count % 1_000_000 == 0 {
            pb.set_message(format!(
                "Processed {:.1}M pairs...",
                count as f64 / 1_000_000.0
            ));
        }
    }

    Ok(count)
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
        } => straw::effres_hic(input.as_path(), chromosome.as_str(), *thr, *pct),
    }
}
