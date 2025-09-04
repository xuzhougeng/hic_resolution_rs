# Repository Guidelines

## Project Structure & Module Organization
- `src/`: Core Rust — `coverage.rs` (binning), `parser.rs` (merged_nodups + pairtools), `resolution.rs` (search), `straw.rs` (.hic tools), `utils.rs`, `main.rs`, `lib.rs`.
- `benches/benchmark.rs`: Criterion benchmarks. `calculate_map_resolution.sh`: original Juicer script.
- CI: `.github/workflows/release.yml` builds Linux-musl/macOS/Windows on `v*` tags.

## Build, Test, and Development Commands
- Build: `cargo build --release` → `target/release/hic_resolution`.
- Resolution (merged_nodups/.pairs): `hic_resolution [OPTIONS] <file>`
  - Examples: `--chrom-size chrom.size demo.txt` | `data/mapped.pairs.gz` (auto-detect header).
- Straw (.hic):
  - List resolutions + chromosomes: `hic_resolution straw list data/example.hic`
  - Dump slice (observed/NONE/BP): `hic_resolution straw dump observed NONE data/example.hic BP 10000 out.slc.gz`
  - Effective resolution: `hic_resolution straw effres data/example.hic chr1 --thr 1000 --pct 0.8`
- Tests: `cargo test` (e.g., `utils` reads `chrom.size`).
- Lint/format: `cargo clippy -- -D warnings` | `cargo fmt --all`.

## Performance Notes & Tuning
- merged_nodups parsing:
  - Uses a zero-copy byte scanner with early filters (frag/mapq first), then parses chr/pos only for passing lines.
  - I/O buffers are increased (`BufReader` 256 KiB) to reduce syscalls.
  - Ingestion batches pairs and performs parallel aggregation with per-worker sparse maps (rayon), then merges into per‑chrom arrays.
- pairs(.gz) parsing:
  - Header is auto-sniffed (`#chromsize:`/`#samheader:`) to obtain names + lengths; fewer fields to parse, typically very fast. Gzip decompression is not a bottleneck compared to parsing and bin updates.
- Bin updates:
  - Bins are plain `u32` (no atomics). Parsing/aggregation is designed so that mutation of the dense bins is single-threaded, while counting is parallel via sparse partials to avoid contention.
- Tuning knobs:
  - Threads: `--threads N` (or `RAYON_NUM_THREADS`) to control parallel aggregation.
  - Chunk sizes: now configurable via CLI
    - `--chunk-pairs <N>`: pairs per aggregation chunk (default 4,000,000; sized to be safe under ~8 GB RAM with typical overhead).
    - `--subchunk-pairs <N>`: per-worker subchunk size (default 128,000).
  - For extremely sparse data, the search step size auto-scales by 10x.

## Feature Flags
- `fast_chrmap` (experimental): switches chromosome name lookup to a custom open-addressing map.
  - Enable with: `cargo build --release --features fast_chrmap`
  - Default remains `FxHashMap`-based lookup, which is generally faster in practice on diverse contig names.

## Chromosome Sizes & Headers
- merged_nodups / stdin / `.pairs` without header:
  - Provide `--chrom-size chrom.size` to supply per‑chrom lengths and build accurate bins. Without it, the program falls back to default hg19 names/lengths (only suitable for human), which will truncate or mis-map non-human contigs (e.g., scaffoldX/contigX/ptg0000X).
- `.pairs` with header:
  - If the file has `#chromsize:` or `#samheader:` lines, names and lengths are auto-detected from the header; no external `chrom.size` needed.
- `.hic` (straw subcommands):
  - Chromosome names and lengths are read from the `.hic` index; no external `chrom.size`.

## Algorithm Notes vs. Juicer Script
- Juicer `calculate_map_resolution.sh` uses a fixed total genome size to form the denominator; this Rust tool uses the sum of provided per‑chrom lengths (`coverage.total_genome_size()`).
- Per‑chrom binning avoids end-of-chromosome artifacts and ensures positions beyond contig length are ignored.
- Resolution search:
  - Coarse step search then binary search; max bin size capped at 10 Mb.
  - Step size inflates automatically for very sparse inputs.

## Troubleshooting
- merged_nodups is slow:
  - Build in release (`cargo build --release`) and run with `--threads` sized to CPU cores.
  - Ensure `--chrom-size` is correct so bins align with your contigs; wrong lengths can cause excessive rejections.
  - If still CPU-bound, consider increasing chunk sizes in `process_pairs` (src/cli.rs) or share a sample for profiling.
- Wrong chromosome mapping:
  - For non-human assemblies with names like `scaffoldX/contigX/ptg0000X`, always supply `--chrom-size` or a `.pairs` file with a proper header.

## Coding Style & Naming Conventions
- Rust 2021; 4-space indent; `snake_case` for items, `CamelCase` for types; doc public APIs with rustdoc.
- Keep modules focused; prefer explicit, descriptive names (e.g., `find_resolution`, `get_matrix_zoom_data`).

## Testing Guidelines
- Unit tests colocated via `#[cfg(test)]`; add integration tests under `tests/` for CLI paths.
- Cover: parser filters, coverage aggregation, search bounds/rounding, straw `.hic` block decode paths.

## Commit & Pull Request Guidelines
- Commits: imperative; ≤72 chars summary; optional scopes `parser:`, `coverage:`, `straw:`, `ci:`.
- PRs: include problem/solution, perf or behavior deltas, sample commands + outputs, and linked issues.

## Security & Configuration Tips
- Resolution CLI: `.pairs` header is auto-detected only from file path (not stdin). For other genomes, pass `--chrom-size`.
- Straw CLI: `.hic` only; unit `BP`, norm `NONE`. Avoid committing large datasets; `.gitignore` excludes common Hi-C files.
