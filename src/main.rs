pub mod coverage;
pub mod parser;
pub mod resolution;
pub mod utils;
pub mod straw;
pub mod filter;
mod cli;

use anyhow::Result;

fn main() -> Result<()> {
    cli::run()
}
