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
    Umi,
}

fn fetch_default_barcode_ranges(barcode_type: &BarcodeType) -> std::ops::Range<usize> {
    match barcode_type {
        BarcodeType::BC1 => 76..84,
        BarcodeType::BC2 => 38..46,
        BarcodeType::BC3 => 0..8,
        BarcodeType::BC4 => 114..119,
        BarcodeType::Umi => 119..127,
    }
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
/// * `barcodes_reference` - A slice of valid barcode sequences to match against.
/// * `n_missmatches` - The maximum number of allowed mismatches.
/// * `ignore_ns` - If true, 'N' characters in the extracted barcode are ignored during matching.
fn identify_best_barcode(
    barcode_extracted: &str,
    barcodes_reference: &[String],
    n_missmatches: usize,
    ignore_ns: bool,
) -> Result<Option<(String, usize)>> {
    if barcodes_reference.is_empty() {
        return Err(anyhow::anyhow!("Reference barcodes list is empty"));
    }

    let extracted_bytes = barcode_extracted.as_bytes();
    let extracted_len = extracted_bytes.len();

    // Only consider reference barcodes with the same length.
    let mut best_idx = None;
    let mut best_distance = usize::MAX;

    for (idx, ref_barcode) in barcodes_reference.iter().enumerate() {
        if ref_barcode.len() != extracted_len {
            continue;
        }

        let ref_bytes = ref_barcode.as_bytes();
        let mut mismatches = 0usize;

        // Early break if mismatches exceed best so far or allowed mismatches
        let prune_limit = best_distance.min(n_missmatches + 1);

        for i in 0..extracted_len {
            if extracted_bytes[i] != ref_bytes[i] {
                if ignore_ns && extracted_bytes[i] == b'N' {
                    continue;
                }
                mismatches += 1;
                if mismatches >= prune_limit {
                    break;
                }
            }
        }

        if mismatches < best_distance {
            best_distance = mismatches;
            best_idx = Some(idx);

            // Perfect match, no need to continue
            if best_distance == 0 {
                break;
            }
        }
    }

    match best_idx {
        Some(idx) if best_distance <= n_missmatches => {
            Ok(Some((barcodes_reference[idx].clone(), best_distance)))
        }
        _ => Ok(None),
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

fn collate_barcodes(
    barcodes_identified: &HashMap<BarcodeType, (String, usize)>,
    allow_missing_barcodes: bool,
) -> Option<String> {
    let mut parts = Vec::with_capacity(4);
    // Order and per-type fallback lengths
    let order = [
        (&BarcodeType::BC1, "NNNNNNNN"),
        (&BarcodeType::BC2, "NNNNNNNN"),
        (&BarcodeType::BC3, "NNNNNNNN"),
        (&BarcodeType::BC4, "NNNNN"),
    ];

    for (bc_type, fallback) in order {
        if let Some((bc, _)) = barcodes_identified.get(bc_type) {
            parts.push(bc.clone());
        } else if allow_missing_barcodes {
            parts.push(fallback.to_string());
        } else {
            return None;
        }
    }

    Some(parts.join("-"))
}

/// Main function to process FASTQ file, extract barcodes, and write results to output file.
/// # Arguments
/// * `read1` - Path to the input FASTQ file for read 1.
/// * `read2` - Path to the input FASTQ file for read 2.
/// * `barcodes` - Path to the CSV file containing valid barcodes.
/// * `output_prefix` - Prefix for the output FASTQ files.
/// * `n_missmatches` - Optional number of allowed mismatches.
/// * `ignore_ns` - Optional flag to ignore 'N' characters in barcode matching.
/// * `allow_missing_barcodes` - Optional flag to allow missing barcodes.
pub fn run<P: AsRef<Path>>(
    read1: P,
    read2: P,
    barcodes: P,
    output_prefix: String,
    n_missmatches: usize,
    ignore_ns: bool,
    allow_missing_barcodes: bool,
) -> Result<()> {
    info!("Processing reads from {:?}", read1.as_ref());
    info!("Processing reads from {:?}", read2.as_ref());
    info!("Using barcodes from {:?}", barcodes.as_ref());
    info!("Output prefix: {}", output_prefix);
    info!("Number of mismatches allowed: {}", n_missmatches);
    info!("Ignore Ns in barcode matching: {}", ignore_ns);
    info!(
        "Will filter reads with missing barcodes: {}",
        !allow_missing_barcodes
    );

    let output_r1 = format!("{output_prefix}_R1.fastq.gz");
    let output_r2 = format!("{output_prefix}_R2.fastq.gz");

    // Make sure the output directory exists
    if let Some(parent) = Path::new(&output_prefix).parent() {
        std::fs::create_dir_all(parent)?;
    }

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
                && let Ok(Some((best_match, distance))) =
                    identify_best_barcode(bc_seq, valid_barcodes, n_missmatches, ignore_ns)
            {
                barcodes_matched.insert(bc_type.clone(), (best_match, distance));
            }
        }

        let cell_barcode = collate_barcodes(&barcodes_matched, allow_missing_barcodes);

        if cell_barcode.is_none() {
            // Skip reads without complete barcodes
            continue;
        }

        let umi = extract_barcode(seq_r2, &BarcodeType::Umi)?;
        let cell_barcode = cell_barcode.unwrap();
        let new_id_for_r2 = format!("{}|{}|{}", id_r2, cell_barcode, umi);
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

    // #[test]
    // fn test_barcode_type_ranges() {
    //     assert_eq!(fetch_default_barcode_ranges(&BarcodeType::BC1), 76..84);
    //     assert_eq!(fetch_default_barcode_ranges(&BarcodeType::BC2), 38..46);
    //     assert_eq!(fetch_default_barcode_ranges(&BarcodeType::BC3), 0..8);
    //     assert_eq!(fetch_default_barcode_ranges(&BarcodeType::BC4), 114..119);
    // }

    // #[test]
    // fn test_extract_barcode_no_slack() {
    //     let sequence = "ATCGATCG" // BC3: 0..8
    //         .to_owned() + &"N".repeat(30)  // padding
    //         + "GCTAGCTA"  // BC2: 38..46
    //         + &"N".repeat(30)  // padding
    //         + "TTAATTAA"  // BC1: 76..84
    //         + &"N".repeat(30)  // padding
    //         + "CCGGC"; // BC4: 114..119

    //     let bc1 = extract_barcode(&sequence, &BarcodeType::BC1).unwrap();
    //     assert_eq!(bc1, "TTAATTAA");

    //     let bc2 = extract_barcode(&sequence, &BarcodeType::BC2).unwrap();
    //     assert_eq!(bc2, "GCTAGCTA");

    //     let bc3 = extract_barcode(&sequence, &BarcodeType::BC3).unwrap();
    //     assert_eq!(bc3, "ATCGATCG");

    //     let bc4 = extract_barcode(&sequence, &BarcodeType::BC4).unwrap();
    //     assert_eq!(bc4, "CCGGC");
    // }

    // #[test]
    // fn test_identify_best_barcode_exact_match() {
    //     let valid_barcodes = vec![
    //         "ATCGATCG".to_string(),
    //         "GCTAGCTA".to_string(),
    //         "TTAATTAA".to_string(),
    //     ];

    #[test]
    fn test_identify_best_barcode_exact_match() {
        let valid_barcodes = vec![
            "ATCGATCG".to_string(),
            "GCTAGCTA".to_string(),
            "TTAATTAA".to_string(),
        ];

        let result = identify_best_barcode("ATCGATCG", &valid_barcodes, 0, false);
        assert!(result.is_ok());
        let (barcode, distance) = result.unwrap().unwrap();
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
        let result = identify_best_barcode("ATCGATCA", &valid_barcodes, 1, false);
        assert!(result.is_ok());
        let (barcode, distance) = result.unwrap().unwrap();
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

        let result = identify_best_barcode("ATCGANNN", &valid_barcodes, 1, false);
        assert!(result.is_ok());
        assert!(result.unwrap().is_none());
    }

    #[test]
    fn test_identify_best_barcode_no_valid_barcodes() {
        let valid_barcodes = vec!["ATCGATCG".to_string()];
        let result = identify_best_barcode("NNNNNNNNNN", &valid_barcodes, 0, false);
        assert!(result.is_ok());
        assert!(result.unwrap().is_none());
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
        let result = identify_best_barcode("ATCGATAA", &valid_barcodes, 2, false);
        assert!(result.is_ok());
        let (barcode, distance) = result.unwrap().unwrap();
        assert_eq!(barcode, "ATCGATCA");
        assert_eq!(distance, 1);
    }

    #[test]
    fn test_identify_best_barcode_with_ns_ignored() {
        let example = "NAATCTGAGTGGCCGATGTTTCGCATCGGCGTACGACTAGATCGCAGGATTCGAGGAGCGTGTGCGAACTCCTGTCTCTTATACACATCTGACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAGGGG";
        let extracted = extract_barcode(example, &BarcodeType::BC3).unwrap();
        let valid_barcodes = vec!["GAATCTGA".to_string()];

        println!("Extracted barcode: {}", extracted);
        println!("Valid barcodes: {:?}", valid_barcodes);

        let result = identify_best_barcode(&extracted, &valid_barcodes, 0, true);

        assert!(result.is_ok());
        assert!(result.as_ref().unwrap().is_some());
        let (barcode, distance) = result.unwrap().unwrap();
        println!("Barcode: {}, Distance: {}", barcode, distance);
        assert_eq!(barcode, "GAATCTGA");
        assert_eq!(distance, 0);
    }

    #[test]
    fn tests_extracts_barcode() {
        let example = "GTACGCAAGTGGCCGATGTTTCGCATCGGCGTACGACTCGAACTTAGGATTCGAGGAGCGTGTGCGAACTCCTGTCTCTTATACACATCTGACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAGGGG";
        let extracted = extract_barcode(example, &BarcodeType::BC3).unwrap();
        assert_eq!(extracted, "GTACGCAA");
        println!("Extracted barcode: {}", extracted);
    }
}
