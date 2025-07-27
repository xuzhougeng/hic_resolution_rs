pub mod utils;
pub mod parser;
pub mod coverage;
pub mod resolution;

use clap::Parser;
use anyhow::Result;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::stdin;
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(name = "hic_resolution")]
#[command(about = "Fast calculation of Hi-C map resolution")]
struct Cli {
    /// Path to merged_nodups file (can be .gz compressed)
    #[arg(value_name = "MERGED_NODUPS")]
    nodups: Option<PathBuf>,
    
    /// Total genome size in base pairs
    #[arg(long, default_value = "2428425688")]
    genome_size: u64,
    
    /// Minimum bin size (base pairs)
    #[arg(long, default_value = "50")]
    bin_width: u32,
    
    /// Proportion of bins that must meet coverage threshold
    #[arg(long, default_value = "0.8")]
    prop: f64,
    
    /// Minimum contacts per bin to be considered "good"
    #[arg(long, default_value = "1000")]
    count_threshold: u32,
    
    /// Step size for initial coarse search
    #[arg(long, default_value = "1000")]
    step_size: u32,
    
    /// Number of threads to use (0 = auto)
    #[arg(short, long, default_value = "0")]
    threads: usize,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    
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
    
    // Create coverage structure
    let coverage = coverage::Coverage::new(args.bin_width);
    println!("Initialized coverage tracking for {} chromosomes", coverage.bins.len());
    
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
        if path.extension().map_or(false, |ext| ext == "gz") {
            let iter = parser::open_file(file)?;
            process_pairs(iter, &coverage, &pb)?
        } else {
            let iter = parser::open_file_uncompressed(file)?;
            process_pairs(iter, &coverage, &pb)?
        }
    } else {
        // Read from stdin
        let iter = parser::open_file(stdin())?;
        process_pairs(iter, &coverage, &pb)?
    };
    
    pb.set_message("Computing resolution...");
    
    // Find resolution
    let resolution = resolution::find_resolution(
        &coverage,
        args.prop,
        args.count_threshold,
        args.step_size,
    );
    
    pb.finish_and_clear();
    
    // Output results
    println!("Processed {} valid pairs", pairs_processed);
    println!();
    println!("Map resolution = {} bp", resolution);
    
    Ok(())
}

fn process_pairs<I>(
    iter: I,
    coverage: &coverage::Coverage,
    pb: &ProgressBar,
) -> Result<u64>
where
    I: Iterator<Item = Result<utils::Pair>>,
{
    let mut count = 0u64;
    
    for pair_result in iter {
        let pair = pair_result?;
        coverage.add_pair(&pair);
        count += 1;
        
        if count % 1_000_000 == 0 {
            pb.set_message(format!("Processed {:.1}M pairs...", count as f64 / 1_000_000.0));
        }
    }
    
    Ok(count)
}