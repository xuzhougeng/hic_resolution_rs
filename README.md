# hickit: Hi-C Toolkit (Rust)

Fast Hi-C toolkit in Rust: map resolution estimation (Juicer-compatible), .hic (straw) helpers, and region filters for merged_nodups.

## Overview

This toolkit provides:

1. Resolution: estimate map resolution from merged_nodups or pairtools .pairs
2. Straw: list/dump/effres for .hic files (observed/NONE/BP)
3. Filter: extract merged_nodups lines overlapping a genomic region

## Features

- **Fast**: 10-100x faster than the original Bash/awk implementation
- **Memory efficient**: Uses atomic counters and parallel processing
- **Flexible**: Supports compressed (.gz) and uncompressed input files
- **Standards compliant**: Compatible with Juicer merged_nodups format

## Installation

Download the pre - compiled version from the release, or compile it manually and add it to the environment variables.

```bash
cargo build --release
```

## Usage

Basic usage with a merged_nodups file:
```bash
hickit resolution merged_nodups.txt
```

With compressed input:
```bash
hickit resolution merged_nodups.txt.gz
```

Reading from stdin:
```bash
zcat merged_nodups.txt.gz | hickit resolution
```

### Resolution Options

- `--genome-size <SIZE>`: Total genome size in bp (default: 2428425688 for hg19)
- `--bin-width <WIDTH>`: Base bin width in bp (default: 50)
- `--prop <PROPORTION>`: Required proportion of good bins (default: 0.8)
- `--count-threshold <COUNT>`: Minimum contacts per bin (default: 1000)
- `--step-size <SIZE>`: Step size for coarse search (default: 1000)
- `--threads <NUM>`: Number of threads (default: auto)

### Examples

```bash
# Human hg38 genome
hickit resolution --genome-size 3137161264 merged_nodups.txt

# Custom parameters  
hickit resolution --prop 0.75 --count-threshold 500 merged_nodups.txt

# With 8 threads
hickit resolution --threads 8 merged_nodups.txt.gz
```

### Pairtools .pairs Usage

The tool also accepts pairtools `.pairs` or `.pairs.gz` files with header lines (e.g., `#chromsize:`):

```bash
# Read directly from a .pairs file
hickit resolution data/mapped.pairs

# Read compressed .pairs.gz
hickit resolution data/mapped.pairs.gz
```

- Chrom sizes are auto-derived from the `.pairs` header; `--chrom-size` is not required.
- As a proxy for mapping quality, only rows with `pair_type == UU` are counted.
- Note: Auto-detection relies on reading the file path. If you use stdin piping for `.pairs`, header detection is skipped; prefer passing the file path directly.

### Filter merged_nodups by region

Extract lines from `merged_nodups(.gz)` where either end overlaps a genomic region and print them to stdout.

```bash
# CHR:START-END form
hickit filter data/merged_nodups.txt ptg000001l:23805-33805 > subset.txt

# CHR START-END form
hickit filter data/merged_nodups.txt ptg000001l 23805-33805 > subset.txt

# Read gzip directly (auto-detected by .gz extension)
hickit filter data/merged_nodups.txt.gz ptg000001l:23805-33805 > subset.txt

# From stdin (decompress yourself if needed)
zcat data/merged_nodups.txt.gz | hickit filter - ptg000001l:23805-33805 > subset.txt
```

- `--uniq`: apply the same uniqueness filter as the main parser (requires `mapq1>0 && mapq2>0` and `frag1!=frag2`).
- Region is inclusive `[start, end]`. Separators `-`, `..`, or `_` are accepted; commas in numbers are allowed (e.g., `23,805-33,805`).
- Outputs matching original lines unmodified, suitable for downstream tools.

## Straw (.hic) Utilities

List resolutions and chromosomes in a `.hic` file:

```bash
hickit straw list data/example.hic
# Outputs:
# Resolutions (BP): 25000, 10000, 5000, ...
# Chromosomes (name\tlength):
# chr1   248956422
```

Dump genome-wide observed counts at a resolution to a slice file (gzip):

```bash
hickit straw dump observed NONE data/example.hic BP 10000 out.slc.gz
```

- Supports local `.hic` files; unit must be `BP` and normalization `NONE`.
- Output slice format: magic `HICSLICE`, `i32` resolution, `i32` chrom count, then per-chrom mapping followed by records `(i16 chr1Key, i32 binX, i16 chr2Key, i32 binY, f32 value)`.

Estimate effective resolution per chromosome (Python reference logic):

```bash
hickit straw effres data/example.hic chr1 --thr 1000 --pct 0.8
# Prints coverage by resolution and the first resolution meeting the threshold
```

- Computes, for each available BP resolution in the `.hic`, the fraction of bins on the chromosome with ≥ `thr` contacts (summing both ends of contacts), and reports the minimum resolution where coverage ≥ `pct`.

## Input Format

The tool expects Juicer merged_nodups format with tab-separated fields:
```
str1  chr1  pos1  frag1  str2  chr2  pos2  frag2  mapq1  cigar1  seq1  mapq2
```

Pairs are counted when `mapq1 > 0`, `mapq2 > 0`, and `frag1 != frag2`. Both intra- and inter-chromosomal pairs are included.

Also supports pairtools `.pairs[.gz]` format with header lines (e.g., `#chromsize:`). For `.pairs` input, the tool auto-detects the header, builds chromosome lengths from it, and parses data rows using columns `chrom1 pos1 chrom2 pos2`. As a proxy for mapping quality, only rows with `pair_type == UU` are used.

## Performance

On a typical workstation with 8 cores and 32GB RAM:

| Input Size | Processing Time |
|------------|----------------|
| 100M pairs | ~30 seconds |
| 1B pairs   | ~5 minutes |

Memory usage scales linearly with genome size (~240MB for human genome).

## Comparison with Original

This implementation produces identical results to the original Bash script but with significant performance improvements:

- **Parsing**: 5-10x faster with optimized field extraction
- **Coverage building**: 3-5x faster with parallel atomic operations  
- **Resolution search**: 2-3x faster with efficient bin aggregation

## License

MIT License (same as original Juicer implementation)
