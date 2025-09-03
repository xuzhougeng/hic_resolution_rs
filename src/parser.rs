use crate::utils::{ChrMap, Pair};
use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::io::Read;
use std::io::{BufRead, BufReader};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

#[derive(Clone, Copy)]
enum ParseMode {
    Juicer,
    Pairs,
}

pub struct PairIterator<R: BufRead> {
    reader: R,
    chr_map: ChrMap,
    buffer: String,
    mode: ParseMode,
}

impl<R: BufRead> PairIterator<R> {
    fn new(reader: R, chr_map: ChrMap, mode: ParseMode) -> Self {
        Self {
            reader,
            chr_map,
            buffer: String::with_capacity(1024),
            mode,
        }
    }
}

impl<R: BufRead> Iterator for PairIterator<R> {
    type Item = Result<Pair>;

    fn next(&mut self) -> Option<Self::Item> {
        static LINE_COUNT: AtomicU64 = AtomicU64::new(0);
        static PARSED_COUNT: AtomicU64 = AtomicU64::new(0);
        static DEBUG_SHOWN: AtomicBool = AtomicBool::new(false);

        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => {
                    if cfg!(debug_assertions) {
                        let line_count = LINE_COUNT.load(Ordering::Relaxed);
                        let parsed_count = PARSED_COUNT.load(Ordering::Relaxed);
                        if line_count > 0 {
                            eprintln!(
                                "Debug: Processed {} lines, parsed {} pairs",
                                line_count, parsed_count
                            );
                        }
                    }
                    return None; // EOF
                }
                Ok(_) => {
                    if let ParseMode::Pairs = self.mode {
                        // Skip header/comment lines
                        if self.buffer.as_bytes().get(0) == Some(&b'#') {
                            continue;
                        }
                    }
                    let line_count = if cfg!(debug_assertions) { LINE_COUNT.fetch_add(1, Ordering::Relaxed) + 1 } else { 0 };
                    if cfg!(debug_assertions) {
                        if !DEBUG_SHOWN.load(Ordering::Relaxed) && line_count <= 3 {
                            eprintln!("Debug line {}: {}", line_count, self.buffer.trim());
                        }
                        if line_count == 3 {
                            DEBUG_SHOWN.store(true, Ordering::Relaxed);
                        }
                    }

                    let parsed = match self.mode {
                        ParseMode::Juicer => parse_line_juicer(&self.buffer, &self.chr_map),
                        ParseMode::Pairs => parse_line_pairs(&self.buffer, &self.chr_map),
                    };

                    if let Some(pair) = parsed {
                        let parsed_count = if cfg!(debug_assertions) { PARSED_COUNT.fetch_add(1, Ordering::Relaxed) + 1 } else { 0 };
                        if cfg!(debug_assertions) {
                            if parsed_count <= 3 {
                                eprintln!(
                                    "Debug: Parsed pair {}: chr{}:{} - chr{}:{}",
                                    parsed_count, pair.chr1, pair.pos1, pair.chr2, pair.pos2
                                );
                            }
                        }
                        return Some(Ok(pair));
                    }
                    // Invalid line, continue to next
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

fn parse_line_juicer(line: &str, chr_map: &ChrMap) -> Option<Pair> {
    // Fast, zero-copy field scanner over ASCII whitespace
    let bytes = line.as_bytes();
    let mut i = 0usize;
    let n = bytes.len();

    // indices we need (0-based tokens):
    // 1(chr1),2(pos1),3(frag1),5(chr2),6(pos2),7(frag2),8(mapq1),11(mapq2 optional)
    let mut f1: Option<(usize, usize)> = None; // chr1
    let mut f2: Option<(usize, usize)> = None; // pos1
    let mut f3: Option<(usize, usize)> = None; // frag1
    let mut f5: Option<(usize, usize)> = None; // chr2
    let mut f6: Option<(usize, usize)> = None; // pos2
    let mut f7: Option<(usize, usize)> = None; // frag2
    let mut f8: Option<(usize, usize)> = None; // mapq1
    let mut f11: Option<(usize, usize)> = None; // mapq2

    let mut tok_idx = 0usize;
    while i < n {
        // skip whitespace (space, tab, CR, LF)
        while i < n {
            let b = bytes[i];
            if b == b' ' || b == b'\t' || b == b'\n' || b == b'\r' { i += 1; } else { break; }
        }
        if i >= n { break; }
        let start = i;
        while i < n {
            let b = bytes[i];
            if b == b' ' || b == b'\t' || b == b'\n' || b == b'\r' { break; }
            i += 1;
        }
        let end = i;
        match tok_idx {
            1 => f1 = Some((start, end)),
            2 => f2 = Some((start, end)),
            3 => f3 = Some((start, end)),
            5 => f5 = Some((start, end)),
            6 => f6 = Some((start, end)),
            7 => f7 = Some((start, end)),
            8 => f8 = Some((start, end)),
            11 => { f11 = Some((start, end)); break; } // we can stop after mapq2
            _ => {}
        }
        tok_idx += 1;
    }

    // Required fields must exist (mapq2 optional, defaults to 0)
    let (s1,e1) = f1?; // chr1
    let (s2,e2) = f2?; // pos1
    let (s3,e3) = f3?; // frag1
    let (s5,e5) = f5?; // chr2
    let (s6,e6) = f6?; // pos2
    let (s7,e7) = f7?; // frag2
    let (s8,e8) = f8?; // mapq1

    // Parse only integers needed for filter first (fast reject path)
    let frag1 = crate::utils::parse_u32_fast(&bytes[s3..e3])?;
    let frag2 = crate::utils::parse_u32_fast(&bytes[s7..e7])?;
    let mapq1 = crate::utils::parse_u32_fast(&bytes[s8..e8])?;
    let mapq2 = if let Some((s,e)) = f11 { crate::utils::parse_u32_fast(&bytes[s..e]).unwrap_or(0) } else { 0 };
    if !(mapq1 > 0 && mapq2 > 0 && frag1 != frag2) {
        return None;
    }

    // Passed filter: now parse chr and positions
    let chr1 = {
        let s = std::str::from_utf8(&bytes[s1..e1]).ok()?;
        *chr_map.get(s)?
    };
    let pos1 = crate::utils::parse_u32_fast(&bytes[s2..e2])?;
    let chr2 = {
        let s = std::str::from_utf8(&bytes[s5..e5]).ok()?;
        *chr_map.get(s)?
    };
    let pos2 = crate::utils::parse_u32_fast(&bytes[s6..e6])?;

    Some(Pair { chr1, pos1, chr2, pos2 })
}

fn parse_line_pairs(line: &str, chr_map: &ChrMap) -> Option<Pair> {
    let line = line.trim_end();
    if line.is_empty() || line.starts_with('#') {
        return None;
    }

    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return None;
    }

    // #columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
    let chr1_str = fields[1];
    let pos1_str = fields[2];
    let chr2_str = fields[3];
    let pos2_str = fields[4];
    let pair_type = fields[7];

    // Heuristic filter to approximate mapq1>0 && mapq2>0: require both uniquely mapped
    if pair_type != "UU" {
        return None;
    }

    let chr1 = chr_map.get(chr1_str).copied()?;
    let pos1 = pos1_str.parse::<u32>().ok()?;
    let chr2 = chr_map.get(chr2_str).copied()?;
    let pos2 = pos2_str.parse::<u32>().ok()?;

    Some(Pair { chr1, pos1, chr2, pos2 })
}

pub fn open_file<R: Read>(
    reader: R,
    chrom_size_file: Option<&str>,
) -> Result<PairIterator<BufReader<MultiGzDecoder<R>>>> {
    let decoder = MultiGzDecoder::new(reader);
    // Larger buffer helps throughput on large text files
    let buf_reader = BufReader::with_capacity(256 * 1024, decoder);
    let chr_map = crate::utils::create_chr_map(chrom_size_file);
    Ok(PairIterator::new(buf_reader, chr_map, ParseMode::Juicer))
}

pub fn open_file_uncompressed<R: Read>(
    reader: R,
    chrom_size_file: Option<&str>,
) -> Result<PairIterator<BufReader<R>>> {
    let buf_reader = BufReader::with_capacity(256 * 1024, reader);
    let chr_map = crate::utils::create_chr_map(chrom_size_file);
    Ok(PairIterator::new(buf_reader, chr_map, ParseMode::Juicer))
}

pub fn open_pairs_file<R: Read>(
    reader: R,
    chr_map: ChrMap,
) -> Result<PairIterator<BufReader<MultiGzDecoder<R>>>> {
    let decoder = MultiGzDecoder::new(reader);
    let buf_reader = BufReader::with_capacity(64 * 1024, decoder);
    Ok(PairIterator::new(buf_reader, chr_map, ParseMode::Pairs))
}

pub fn open_pairs_file_uncompressed<R: Read>(
    reader: R,
    chr_map: ChrMap,
) -> Result<PairIterator<BufReader<R>>> {
    let buf_reader = BufReader::with_capacity(64 * 1024, reader);
    Ok(PairIterator::new(buf_reader, chr_map, ParseMode::Pairs))
}

use std::path::Path;
pub fn sniff_pairs_header_from_path(path: &Path) -> Result<Option<(ChrMap, Vec<u32>)>> {
    use std::fs::File;
    let file = File::open(path)?;
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    if is_gz {
        sniff_pairs_header(MultiGzDecoder::new(file))
    } else {
        sniff_pairs_header(file)
    }
}

fn sniff_pairs_header<R: Read>(reader: R) -> Result<Option<(ChrMap, Vec<u32>)>> {
    let mut reader = BufReader::with_capacity(64 * 1024, reader);
    let mut buf = String::new();
    let mut lengths: Vec<u32> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    use std::collections::HashMap;
    let mut index_of: HashMap<String, usize> = HashMap::new();
    let mut seen_any = false;

    // Read a limited number of header lines to avoid slurping large files
    for _ in 0..2000 {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            break;
        }
        let line = buf.trim_end();
        if !line.starts_with('#') {
            break;
        }
        seen_any = true;
        if let Some(rest) = line.strip_prefix("#chromsize:") {
            let parts: Vec<&str> = rest.trim().split_whitespace().collect();
            if parts.len() >= 2 {
                if let Ok(len) = parts[1].parse::<u32>() {
                    let name = parts[0].to_string();
                    if let Some(&idx) = index_of.get(&name) {
                        lengths[idx] = len;
                    } else {
                        index_of.insert(name.clone(), names.len());
                        names.push(name);
                        lengths.push(len);
                    }
                }
            }
        } else if let Some(rest) = line.strip_prefix("#samheader:") {
            // Example: @SQ SN:Chr1 LN:30427671
            if let Some(sq) = rest.split_whitespace().find(|s| s.starts_with("@SQ")) {
                let _ = sq; // unused
            }
            let mut name_opt: Option<String> = None;
            let mut len_opt: Option<u32> = None;
            for token in rest.split_whitespace() {
                if let Some(v) = token.strip_prefix("SN:") {
                    name_opt = Some(v.to_string());
                }
                if let Some(v) = token.strip_prefix("LN:") {
                    if let Ok(l) = v.parse::<u32>() {
                        len_opt = Some(l);
                    }
                }
            }
            if let (Some(nm), Some(ln)) = (name_opt, len_opt) {
                if let Some(&idx) = index_of.get(&nm) {
                    lengths[idx] = ln;
                } else {
                    index_of.insert(nm.clone(), names.len());
                    names.push(nm);
                    lengths.push(ln);
                }
            }
        }
    }

    if !lengths.is_empty() {
        let mut map = ChrMap::default();
        for (i, nm) in names.into_iter().enumerate() {
            let idx = (i + 1) as u8;
            map.insert(nm, idx);
        }
        Ok(Some((map, lengths)))
    } else if seen_any {
        Ok(None) // header present but no lengths parsed
    } else {
        Ok(None)
    }
}
