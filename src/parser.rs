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
                    let line_count = LINE_COUNT.load(Ordering::Relaxed);
                    let parsed_count = PARSED_COUNT.load(Ordering::Relaxed);
                    if line_count > 0 {
                        eprintln!(
                            "Debug: Processed {} lines, parsed {} pairs",
                            line_count, parsed_count
                        );
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
                    let line_count = LINE_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
                    if !DEBUG_SHOWN.load(Ordering::Relaxed) && line_count <= 3 {
                        eprintln!("Debug line {}: {}", line_count, self.buffer.trim());
                    }
                    if line_count == 3 {
                        DEBUG_SHOWN.store(true, Ordering::Relaxed);
                    }

                    let parsed = match self.mode {
                        ParseMode::Juicer => parse_line_juicer(&self.buffer, &self.chr_map),
                        ParseMode::Pairs => parse_line_pairs(&self.buffer, &self.chr_map),
                    };

                    if let Some(pair) = parsed {
                        let parsed_count = PARSED_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
                        if parsed_count <= 3 {
                            eprintln!(
                                "Debug: Parsed pair {}: chr{}:{} - chr{}:{}",
                                parsed_count, pair.chr1, pair.pos1, pair.chr2, pair.pos2
                            );
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
    let line = line.trim_end();

    // Split by whitespace (spaces, not tabs in this format)
    let fields: Vec<&str> = line.split_whitespace().collect();

    static PARSE_ATTEMPT_COUNT: AtomicU64 = AtomicU64::new(0);
    let parse_count = PARSE_ATTEMPT_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
    if parse_count <= 3 {
        eprintln!(
            "Debug parse_line {}: {} fields, line='{}'",
            parse_count,
            fields.len(),
            line
        );
    }

    if fields.len() < 9 {
        static FIELD_ERROR_COUNT: AtomicU64 = AtomicU64::new(0);
        let error_count = FIELD_ERROR_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
        if error_count <= 3 {
            eprintln!(
                "Debug: Line rejected - only {} fields (need 9+)",
                fields.len()
            );
        }
        return None;
    }

    // Actual field mapping from demo.txt analysis:
    // Field 1: chr1, Field 2: pos1, Field 3: frag1, Field 4: str1
    // Field 5: chr2, Field 6: pos2, Field 7: frag2, Field 8: str2
    // Field 9: mapq1, Field 12: mapq2
    // Original script: chr1=$2, pos1=$3, frag1=$4, chr2=$6, pos2=$7, frag2=$8, mapq1=$9, mapq2=$12

    let chr1_str = fields[1]; // Field 2 in awk = index 1
    let pos1_str = fields[2]; // Field 3 in awk = index 2
    let frag1_str = fields[3]; // Field 4 in awk = index 3
    let chr2_str = fields[5]; // Field 6 in awk = index 5
    let pos2_str = fields[6]; // Field 7 in awk = index 6
    let frag2_str = fields[7]; // Field 8 in awk = index 7
    let mapq1_str = fields[8]; // Field 9 in awk = index 8
    let mapq2_str = if fields.len() > 11 { fields[11] } else { "0" }; // Field 12 in awk = index 11

    // Parse values with detailed error reporting
    let chr1 = match chr_map.get(chr1_str).copied() {
        Some(c) => c,
        None => {
            static CHR_ERROR_COUNT: AtomicU64 = AtomicU64::new(0);
            let error_count = CHR_ERROR_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
            if error_count <= 3 {
                eprintln!("Debug: Failed to parse chr1='{}' (not in map)", chr1_str);
            }
            return None;
        }
    };

    let pos1 = match pos1_str.parse() {
        Ok(p) => p,
        Err(_) => {
            static POS1_ERROR_COUNT: AtomicU64 = AtomicU64::new(0);
            let error_count = POS1_ERROR_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
            if error_count <= 3 {
                eprintln!("Debug: Failed to parse pos1='{}'", pos1_str);
            }
            return None;
        }
    };

    let frag1 = frag1_str.parse::<u32>().ok()?;
    let chr2 = chr_map.get(chr2_str).copied()?;
    let pos2 = pos2_str.parse().ok()?;
    let frag2 = frag2_str.parse::<u32>().ok()?;
    let mapq1 = mapq1_str.parse::<u32>().ok()?;
    let mapq2 = mapq2_str.parse::<u32>().ok()?;

    // Debug: Show what we parsed for first few lines
    static DEBUG_PARSE_COUNT: AtomicU64 = AtomicU64::new(0);
    let debug_count = DEBUG_PARSE_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
    if debug_count <= 3 {
        eprintln!("Debug parse {}: chr1={}, pos1={}, frag1={}, chr2={}, pos2={}, frag2={}, mapq1={}, mapq2={}", 
                 debug_count, chr1_str, pos1_str, frag1_str, chr2_str, pos2_str, frag2_str, mapq1_str, mapq2_str);
        eprintln!(
            "  Filter check: mapq1={}>0? mapq2={}>0? frag1={}!=frag2={}?",
            mapq1, mapq2, frag1, frag2
        );
    }

    // Apply filters from original script: $9>0 && $12>0 && $4!=$8
    if mapq1 > 0 && mapq2 > 0 && frag1 != frag2 {
        static ACCEPTED_COUNT: AtomicU64 = AtomicU64::new(0);
        let accepted_count = ACCEPTED_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
        if accepted_count <= 3 {
            eprintln!("Debug: Accepted pair {}: chr{}:{} - chr{}:{} (mapq1={}, mapq2={}, frag1={}, frag2={})", 
                     accepted_count, chr1, pos1, chr2, pos2, mapq1, mapq2, frag1, frag2);
        }
        Some(Pair {
            chr1,
            pos1,
            chr2,
            pos2,
        })
    } else {
        static FILTERED_COUNT: AtomicU64 = AtomicU64::new(0);
        let filtered_count = FILTERED_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
        if filtered_count <= 10 {
            eprintln!(
                "Debug: Filtered pair {} - mapq1={}>0? mapq2={}>0? frag1={}!=frag2={}?",
                filtered_count, mapq1, mapq2, frag1, frag2
            );
        }
        None
    }
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
    let buf_reader = BufReader::with_capacity(64 * 1024, decoder);
    let chr_map = crate::utils::create_chr_map(chrom_size_file);
    Ok(PairIterator::new(buf_reader, chr_map, ParseMode::Juicer))
}

pub fn open_file_uncompressed<R: Read>(
    reader: R,
    chrom_size_file: Option<&str>,
) -> Result<PairIterator<BufReader<R>>> {
    let buf_reader = BufReader::with_capacity(64 * 1024, reader);
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
