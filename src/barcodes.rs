use anyhow::Result;
use log::info;
use noodles::fastq;
use noodles::fastq::record::Definition;
use polars::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, BufWriter};
use std::path::Path;

// Type aliases to simplify complex types
/// A dynamic FASTQ reader that can handle both compressed and uncompressed files
type FastqReader = fastq::io::Reader<Box<dyn std::io::BufRead>>;
/// A dynamic FASTQ writer that can handle both compressed and uncompressed files
type FastqWriter = fastq::io::Writer<Box<dyn std::io::Write>>;

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
enum BarcodeType {
    BC1,
    BC2,
    BC3,
    BC4,
}

fn fetch_default_barcode_ranges(barcode_type: &BarcodeType) -> std::ops::Range<usize> {
    match barcode_type {
        BarcodeType::BC1 => 76..84,
        BarcodeType::BC2 => 38..46,
        BarcodeType::BC3 => 0..8,
        BarcodeType::BC4 => 114..119,
    }
}

#[derive(Debug)]
struct SplitPoolBarcodes {
    barcode_type: BarcodeType,
    barcode_sequences: Vec<String>,
    barcode_range: std::ops::Range<usize>,
}
impl SplitPoolBarcodes {
    pub fn new(
        barcode_type: BarcodeType,
        barcode_sequences: Vec<String>,
        barcode_range: std::ops::Range<usize>,
    ) -> Self {
        SplitPoolBarcodes {
            barcode_type,
            barcode_sequences,
            barcode_range,
        }
    }
}

struct Barcodes {
    bc1: SplitPoolBarcodes,
    bc2: SplitPoolBarcodes,
    bc3: SplitPoolBarcodes,
    bc4: SplitPoolBarcodes,
}

fn load_barcodes<P: AsRef<Path>>(path: P) -> Result<HashMap<BarcodeType, Vec<String>>> {
    let df = CsvReadOptions::default()
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(path.as_ref().to_path_buf()))?
        .finish()?;

    let barcode_map: HashMap<BarcodeType, Vec<String>> = df
        .partition_by(["barcode_type"], true)?
        .into_iter()
        .map(|df| {
            let barcode_type = df
                .column("barcode_type")
                .expect("Failed to get 'barcode_type' column")
                .str()
                .unwrap()
                .get(0)
                .unwrap()
                .to_string();

            let barcode_type: BarcodeType = match barcode_type.as_str() {
                "BC1" => BarcodeType::BC1,
                "BC2" => BarcodeType::BC2,
                "BC3" => BarcodeType::BC3,
                "BC4" => BarcodeType::BC4,
                _ => panic!("Unknown barcode type"),
            };

            let barcodes: Vec<String> = df
                .column("barcode_sequence")
                .expect("Failed to get 'barcode_sequence' column")
                .str()
                .unwrap()
                .into_no_null_iter()
                .map(|s| s.to_string())
                .collect();
            (barcode_type, barcodes)
        })
        .collect();

    Ok(barcode_map)
}

fn extract_barcode(sequence: &str, barcode_type: &BarcodeType) -> Result<String> {
    let range = fetch_default_barcode_ranges(barcode_type);

    if range.end > sequence.len() {
        return Err(anyhow::anyhow!("Barcode range is out of bounds"));
    };

    Ok(sequence[range].to_string())
}

/// Identify the best matching barcode from a list of valid barcodes using Hamming distance.
/// Returns the best matching barcode and its Hamming distance if found, otherwise None.
/// # Arguments
/// * `barcode_extracted` - The extracted barcode sequence to be matched.
/// * `reference_barcodes` - A slice of valid barcode sequences to match against.
/// * `n_missmatches` - The minimum Hamming distance score to consider a match valid.
/// # Returns
/// An Option containing a tuple of the best matching barcode and its Hamming distance, or None if no match is found.
/// Small helper struct (kept local to function scope) if you want to store intermediate best.
struct BarcodeMatch {
    idx: usize,
    distance: usize,
}

fn identify_best_barcode(
    barcode_extracted: &str,
    barcodes_reference: &[String],
    n_missmatches: usize,
) -> Option<(String, usize)> {
    if barcodes_reference.is_empty() {
        return None;
    }

    let extracted_bytes = barcode_extracted.as_bytes();
    let extracted_len = extracted_bytes.len();

    // best match index into barcodes_reference + distance
    let mut best: Option<BarcodeMatch> = None;

    for (idx, ref_barcode) in barcodes_reference.iter().enumerate() {
        let ref_bytes = ref_barcode.as_bytes();
        let ref_len = ref_bytes.len();
        let max_len = extracted_len.max(ref_len);

        // fast prune: if we already have a perfect match, stop
        if let Some(b) = &best
            && b.distance == 0
        {
            break;
        }

        // Count mismatches, but stop early if >= current best distance
        let mut mismatches = 0usize;
        // compute a prune_upper bound = if we already have a best match, any candidate
        // whose mismatches >= best.distance can be ignored. Use that to break early.
        let prune_limit = best.as_ref().map(|b| b.distance).unwrap_or(usize::MAX);

        for i in 0..max_len {
            let a = if i < extracted_len {
                extracted_bytes[i]
            } else {
                b'N'
            };
            let b = if i < ref_len { ref_bytes[i] } else { b'N' };
            if a != b {
                mismatches += 1;
                if mismatches >= prune_limit {
                    // no hope to beat current best, break early
                    break;
                }
            }
        }

        // update best if better
        if best.is_none() || mismatches < best.as_ref().unwrap().distance {
            best = Some(BarcodeMatch {
                idx,
                distance: mismatches,
            });
            // perfect match -> we can stop searching early
            if mismatches == 0 {
                break;
            }
        }
    }

    // no candidate matched (shouldn't happen because we checked empty earlier)
    let best = best?;

    // Now check threshold using length difference to avoid penalizing for missing tails
    let best_ref = &barcodes_reference[best.idx];
    let length_diff = barcode_extracted.len().abs_diff(best_ref.len());
    let max_allowed = n_missmatches + length_diff;

    if best.distance > max_allowed {
        None
    } else {
        // clone only the chosen barcode once
        Some((best_ref.clone(), best.distance))
    }
}

/// Return a FASTQ reader and a writer based on file extensions.
/// Supports uncompressed `.fastq` / `.fq` and gzipped `.gz`.
pub fn get_reader_and_writer<P1: AsRef<Path>, P2: AsRef<Path>>(
    input: P1,
    output: P2,
) -> Result<(FastqReader, FastqWriter)> {
    let input_path = input.as_ref();
    let output_path = output.as_ref();

    // --- Reader ---
    let reader_inner: Box<dyn std::io::BufRead> = match input_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase())
    {
        Some(ext) if ext == "gz" => {
            // compressed: wrap decoder in a BufReader so we provide BufRead
            let file = File::open(input_path)?;
            let decoder = flate2::read::MultiGzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        }
        Some(ext) if ext == "fastq" || ext == "fq" => {
            let file = File::open(input_path)?;
            Box::new(BufReader::new(file))
        }
        Some(ext) => {
            return Err(anyhow::anyhow!(
                "Unsupported input file extension: {:?}",
                ext
            ));
        }
        None => return Err(anyhow::anyhow!("Input file has no extension")),
    };

    let reader = fastq::io::Reader::new(reader_inner);

    // --- Writer ---
    let writer: Box<dyn std::io::Write> = match output_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase())
    {
        Some(ext) if ext == "gz" => {
            let file = File::create(output_path)?;
            let encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            Box::new(encoder)
        }
        Some(ext) if ext == "fastq" || ext == "fq" => {
            let file = File::create(output_path)?;
            Box::new(BufWriter::new(file))
        }
        Some(ext) => {
            return Err(anyhow::anyhow!(
                "Unsupported output file extension: {:?}",
                ext
            ));
        }
        None => return Err(anyhow::anyhow!("Output file has no extension")),
    };

    let writer = fastq::io::Writer::new(writer);

    Ok((reader, writer))
}

/// Main function to process FASTQ file, extract barcodes, and write results to output file.
/// # Arguments
/// * `path` - Path to the input FASTQ file.
/// * `barcodes` - Path to the CSV file containing valid barcodes.
/// * `outfile` - Path to the output file where results will be written.
/// * `slack` - Optional slack values for barcode extraction.
/// * `minimum_score` - Optional minimum Hamming distance score for barcode matching.
pub fn run<P: AsRef<Path>>(
    read1: P,
    read2: P,
    barcodes: P,
    output_prefix: String,
    n_missmatches: Option<usize>,
) -> Result<()> {
    let n_missmatches = n_missmatches.unwrap_or(0);

    let output_r1 = format!("{output_prefix}_R1.fastq.gz");
    let output_r2 = format!("{output_prefix}_R2.fastq.gz");
    let (mut reader_r1, mut writer_r1) = get_reader_and_writer(read1, output_r1)?;
    let (mut reader_r2, mut writer_r2) = get_reader_and_writer(read2, output_r2)?;

    info!("Loading barcodes from {:?}", barcodes.as_ref());
    let barcode_map = load_barcodes(barcodes)?;

    // Pre-allocate reusable buffers outside the loop
    let mut barcodes_found: HashMap<BarcodeType, String> = HashMap::new();
    let mut barcodes_matched: HashMap<BarcodeType, (String, usize)> = HashMap::new();

    for (ii, (result_r1, result_r2)) in reader_r1.records().zip(reader_r2.records()).enumerate() {
        let record_r1 = result_r1?;
        let record_r2 = result_r2?;
        let seq_r2 = std::str::from_utf8(record_r2.sequence())?;
        let id_r1 = std::str::from_utf8(record_r1.name())?;
        let id_r2 = std::str::from_utf8(record_r2.name())?;

        if ii % 1_000_000 == 0 {
            info!("Processed {ii} reads");
        }

        // Clear and reuse HashMaps
        barcodes_found.clear();
        barcodes_matched.clear();

        for bc_type in &[
            BarcodeType::BC1,
            BarcodeType::BC2,
            BarcodeType::BC3,
            BarcodeType::BC4,
        ] {
            let bc = extract_barcode(seq_r2, bc_type)?;
            barcodes_found.insert(bc_type.clone(), bc);
        }
        for (bc_type, bc_seq) in &barcodes_found {
            if let Some(valid_barcodes) = barcode_map.get(bc_type)
                && let Some((best_match, distance)) =
                    identify_best_barcode(bc_seq, valid_barcodes, n_missmatches)
            {
                barcodes_matched.insert(bc_type.clone(), (best_match, distance));
            }
        }

        // Create a new identifier with matched barcodes in the read name. Format ID:BC1-BC2-BC3-BC4
        let new_id_for_r2 = format!(
            "{}:{}-{}-{}-{}",
            id_r2,
            barcodes_matched
                .get(&BarcodeType::BC1)
                .map(|(bc, score)| format!("{bc}_{score}"))
                .unwrap_or_else(|| String::from("NNNNNNNN")),
            barcodes_matched
                .get(&BarcodeType::BC2)
                .map(|(bc, score)| format!("{bc}_{score}"))
                .unwrap_or_else(|| String::from("NNNNNNNN")),
            barcodes_matched
                .get(&BarcodeType::BC3)
                .map(|(bc, score)| format!("{bc}_{score}"))
                .unwrap_or_else(|| String::from("NNNNNNNN")),
            barcodes_matched
                .get(&BarcodeType::BC4)
                .map(|(bc, score)| format!("{bc}_{score}"))
                .unwrap_or_else(|| String::from("NNNNN")),
        );

        let new_id_for_r1 = new_id_for_r2.clone().replace(id_r2, id_r1);

        // Trim the sequence to remove barcodes
        // We have ~138 bp of barcodes on the 5' end of the read so just remove the first 138 bp
        let r2_trimmed_seq = &seq_r2[138..];
        let quality_scores = &record_r2.quality_scores()[138..];

        let r1_record = fastq::Record::new(
            Definition::new(new_id_for_r1, record_r1.description()),
            record_r1.sequence(),
            record_r1.quality_scores(),
        );
        let r2_record = fastq::Record::new(
            Definition::new(new_id_for_r2, record_r2.description()),
            r2_trimmed_seq,
            quality_scores,
        );

        writer_r1.write_record(&r1_record)?;
        writer_r2.write_record(&r2_record)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barcode_type_ranges() {
        assert_eq!(fetch_default_barcode_ranges(&BarcodeType::BC1), 76..84);
        assert_eq!(fetch_default_barcode_ranges(&BarcodeType::BC2), 38..46);
        assert_eq!(fetch_default_barcode_ranges(&BarcodeType::BC3), 0..8);
        assert_eq!(fetch_default_barcode_ranges(&BarcodeType::BC4), 114..119);
    }

    #[test]
    fn test_extract_barcode_no_slack() {
        let sequence = "ATCGATCG" // BC3: 0..8
            .to_owned() + &"N".repeat(30)  // padding
            + "GCTAGCTA"  // BC2: 38..46
            + &"N".repeat(30)  // padding
            + "TTAATTAA"  // BC1: 76..84
            + &"N".repeat(30)  // padding
            + "CCGGC"; // BC4: 114..119

        let bc1 = extract_barcode(&sequence, &BarcodeType::BC1).unwrap();
        assert_eq!(bc1, "TTAATTAA");

        let bc2 = extract_barcode(&sequence, &BarcodeType::BC2).unwrap();
        assert_eq!(bc2, "GCTAGCTA");

        let bc3 = extract_barcode(&sequence, &BarcodeType::BC3).unwrap();
        assert_eq!(bc3, "ATCGATCG");

        let bc4 = extract_barcode(&sequence, &BarcodeType::BC4).unwrap();
        assert_eq!(bc4, "CCGGC");
    }

    #[test]
    fn test_identify_best_barcode_exact_match() {
        let valid_barcodes = vec![
            "ATCGATCG".to_string(),
            "GCTAGCTA".to_string(),
            "TTAATTAA".to_string(),
        ];

        let result = identify_best_barcode("ATCGATCG", &valid_barcodes, 0);
        assert!(result.is_some());
        let (barcode, distance) = result.unwrap();
        assert_eq!(barcode, "ATCGATCG");
        assert_eq!(distance, 0);
    }

    #[test]
    fn test_identify_best_barcode_with_mismatch() {
        let valid_barcodes = vec![
            "ATCGATCG".to_string(),
            "GCTAGCTA".to_string(),
            "TTAATTAA".to_string(),
        ];

        // One mismatch: "ATCGATCG" vs "ATCGATCA"
        let result = identify_best_barcode("ATCGATCA", &valid_barcodes, 1);
        assert!(result.is_some());
        let (barcode, distance) = result.unwrap();
        assert_eq!(barcode, "ATCGATCG");
        assert_eq!(distance, 1);
    }

    #[test]
    fn test_identify_best_barcode_too_many_mismatches() {
        let valid_barcodes = vec![
            "ATCGATCG".to_string(),
            "GCTAGCTA".to_string(),
            "TTAATTAA".to_string(),
        ];

        let result = identify_best_barcode("ATCGANNN", &valid_barcodes, 1);
        println!("{:?}", result);
        assert!(result.is_none());
    }

    #[test]
    fn test_identify_best_barcode_no_valid_barcodes() {
        let valid_barcodes = vec!["ATCGATCG".to_string()];
        let result = identify_best_barcode("NNNNNNNNNN", &valid_barcodes, 0);
        println!("{:?}", result);
        assert!(result.is_none());
    }

    #[test]
    fn test_identify_best_barcode_chooses_closest() {
        let valid_barcodes = vec![
            "ATCGATCG".to_string(),
            "ATCGATCA".to_string(), // 2 mismatches from test sequence
            "TTAATTAA".to_string(),
        ];

        // "ATCGATAA" vs "ATCGATCA" = 1 mismatch (position 6: T vs C)
        // Should match "ATCGATCA" with distance 1
        let result = identify_best_barcode("ATCGATAA", &valid_barcodes, 2);
        assert!(result.is_some());
        let (barcode, distance) = result.unwrap();
        assert_eq!(barcode, "ATCGATCA");
        assert_eq!(distance, 1);
    }
}
