use anyhow::{anyhow, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub struct Region<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
}

impl<'a> Region<'a> {
    pub fn parse(region_or_chrom: &'a str, maybe_range: Option<&'a str>) -> Result<Self> {
        if let Some(span) = maybe_range {
            let (start, end) = parse_span(span)?;
            return Ok(Region { chrom: region_or_chrom, start, end });
        }

        // Accept forms like: chr:start-end, chr:start..end, chr:start_end
        let mut chrom = region_or_chrom;
        let mut start_end: Option<&str> = None;
        if let Some((c, rest)) = region_or_chrom.split_once(':') {
            chrom = c;
            start_end = Some(rest);
        }
        let se = start_end.ok_or_else(|| anyhow!("Region must be CHR:START-END or CHR START-END"))?;
        let (start, end) = parse_span(se)?;
        Ok(Region { chrom, start, end })
    }
}

fn parse_span(se: &str) -> Result<U32Pair> {
    // Support 12345-67890, 12345..67890, 12345_67890
    let (a, b) = if let Some((a, b)) = se.split_once('-') {
        (a, b)
    } else if let Some((a, b)) = se.split_once("..") {
        (a, b)
    } else if let Some((a, b)) = se.split_once('_') {
        (a, b)
    } else {
        return Err(anyhow!("Invalid span: expected START-END"));
    };
    let start: u32 = a.replace(",", "").parse()?;
    let end: u32 = b.replace(",", "").parse()?;
    if start > end {
        Err(anyhow!("Region start > end"))
    } else {
        Ok((start, end))
    }
}

type U32Pair = (u32, u32);

pub struct FilterOptions<'a> {
    pub region: Region<'a>,
    pub require_unique: bool,
}

/// Filter a merged_nodups(.gz) stream, emitting lines where either end overlaps the region.
pub fn filter_merged_nodups_stream<R: Read, W: Write>(
    reader: R,
    opts: &FilterOptions,
    mut out: W,
) -> Result<()> {
    let mut buf_reader = BufReader::with_capacity(256 * 1024, reader);
    let mut line = String::with_capacity(1024);
    let chrom = opts.region.chrom;
    let start = opts.region.start;
    let end = opts.region.end;
    let require_unique = opts.require_unique;

    loop {
        line.clear();
        let n = buf_reader.read_line(&mut line)?;
        if n == 0 { break; }
        if line.trim().is_empty() { continue; }

        if line_matches_region(&line, chrom, start, end, require_unique) {
            out.write_all(line.as_bytes())?;
        }
    }
    out.flush()?;
    Ok(())
}

#[inline]
fn line_matches_region(line: &str, chrom: &str, start: u32, end: u32, require_unique: bool) -> bool {
    // Fast field scanner similar to parser::parse_line_juicer
    let b = line.as_bytes();
    let mut i = 0usize;
    let n = b.len();

    // Fields needed: 1(chr1),2(pos1),3(frag1),5(chr2),6(pos2),7(frag2),8(mapq1),11(mapq2 optional)
    let mut f1: Option<(usize, usize)> = None;
    let mut f2: Option<(usize, usize)> = None;
    let mut f3: Option<(usize, usize)> = None;
    let mut f5: Option<(usize, usize)> = None;
    let mut f6: Option<(usize, usize)> = None;
    let mut f7: Option<(usize, usize)> = None;
    let mut f8: Option<(usize, usize)> = None;
    let mut f11: Option<(usize, usize)> = None;
    let mut tok = 0usize;
    while i < n {
        while i < n {
            let c = b[i];
            if c == b' ' || c == b'\t' || c == b'\n' || c == b'\r' { i += 1; } else { break; }
        }
        if i >= n { break; }
        let s = i;
        while i < n {
            let c = b[i];
            if c == b' ' || c == b'\t' || c == b'\n' || c == b'\r' { break; }
            i += 1;
        }
        let e = i;
        match tok {
            1 => f1 = Some((s, e)),
            2 => f2 = Some((s, e)),
            3 => f3 = Some((s, e)),
            5 => f5 = Some((s, e)),
            6 => f6 = Some((s, e)),
            7 => f7 = Some((s, e)),
            8 => f8 = Some((s, e)),
            11 => { f11 = Some((s, e)); break; }
            _ => {}
        }
        tok += 1;
    }

    let (s1, e1) = match f1 { Some(v) => v, None => return false };
    let (s2, e2) = match f2 { Some(v) => v, None => return false };
    let (s5, e5) = match f5 { Some(v) => v, None => return false };
    let (s6, e6) = match f6 { Some(v) => v, None => return false };

    if require_unique {
        // Apply same early filter as main parser: frag1 != frag2 and mapq1>0 && mapq2>0
        let ok = match (f3, f7, f8) {
            (Some((fs, fe)), Some((gs, ge)), Some((ms, me))) => {
                let frag1 = crate::utils::parse_u32_fast(&b[fs..fe]).unwrap_or(0);
                let frag2 = crate::utils::parse_u32_fast(&b[gs..ge]).unwrap_or(0);
                let mapq1 = crate::utils::parse_u32_fast(&b[ms..me]).unwrap_or(0);
                let mapq2 = if let Some((qs, qe)) = f11 { crate::utils::parse_u32_fast(&b[qs..qe]).unwrap_or(0) } else { 0 };
                mapq1 > 0 && mapq2 > 0 && frag1 != frag2
            }
            _ => false,
        };
        if !ok { return false; }
    }

    // Compare chr names and positions; treat inclusive range
    let chr1 = unsafe { std::str::from_utf8_unchecked(&b[s1..e1]) };
    let chr2 = unsafe { std::str::from_utf8_unchecked(&b[s5..e5]) };
    let pos1 = crate::utils::parse_u32_fast(&b[s2..e2]).unwrap_or(u32::MAX);
    let pos2 = crate::utils::parse_u32_fast(&b[s6..e6]).unwrap_or(u32::MAX);

    (chr1 == chrom && pos1 >= start && pos1 <= end)
        || (chr2 == chrom && pos2 >= start && pos2 <= end)
}

pub fn run_filter_file(input: Option<&Path>, region: Region<'_>, require_unique: bool) -> Result<()> {
    let opts = FilterOptions { region, require_unique };
    let stdout = io::stdout();
    let handle = stdout.lock();
    match input {
        Some(path) => {
            if path.as_os_str() == "-" { 
                let stdin = io::stdin();
                let lock = stdin.lock();
                return filter_merged_nodups_stream(lock, &opts, handle);
            }
            let is_gz = path.extension().and_then(|e| e.to_str()).map(|e| e.eq_ignore_ascii_case("gz")).unwrap_or(false);
            let file = File::open(path)?;
            if is_gz { filter_merged_nodups_stream(MultiGzDecoder::new(file), &opts, handle) }
            else { filter_merged_nodups_stream(file, &opts, handle) }
        }
        None => {
            // stdin (assume plain text)
            let stdin = io::stdin();
            let lock = stdin.lock();
            filter_merged_nodups_stream(lock, &opts, handle)
        }
    }
}
