use anyhow::Result;
use bio::alignment::distance::hamming;
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

pub struct Slack {
    left: usize,
    right: usize,
}

impl Slack {
    pub fn new(left: usize, right: usize) -> Self {
        Slack { left, right }
    }
}

fn get_barcode_range(barcode_type: &BarcodeType) -> std::ops::Range<usize> {
    match barcode_type {
        BarcodeType::BC1 => 76..84,
        BarcodeType::BC2 => 38..46,
        BarcodeType::BC3 => 0..8,
        BarcodeType::BC4 => 114..119,
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

fn extract_barcode(sequence: &str, barcode_type: &BarcodeType, slack: &Slack) -> Result<String> {
    let mut range = get_barcode_range(barcode_type);

    if slack.left > range.start || range.end + slack.right > sequence.len() {
        return Err(anyhow::anyhow!("Barcode range with slack is out of bounds"));
    };
    range.start -= slack.left;
    range.end += slack.right;
    Ok(sequence[range].to_string())
}

struct BarcodeMatch {
    barcode: String,
    distance: usize,
}

/// Identify the best matching barcode from a list of valid barcodes using Hamming distance.
/// Returns the best matching barcode and its Hamming distance if found, otherwise None.
/// # Arguments
/// * `barcode_extracted` - The extracted barcode sequence to be matched.
/// * `reference_barcodes` - A slice of valid barcode sequences to match against.
/// * `n_missmatches` - The minimum Hamming distance score to consider a match valid.
/// # Returns
/// An Option containing a tuple of the best matching barcode and its Hamming distance, or None if no match is found.
fn identify_best_barcode(
    barcode_extracted: &str,
    barcodes_reference: &[String],
    n_missmatches: usize,
) -> Option<(String, usize)> {
    let extracted_length = barcode_extracted.len();

    let mut best_match: Option<BarcodeMatch> = None;
    for barcode_reference in barcodes_reference {
        // Check if the lengths are the same, if not pad the shorter one with Ns
        // This is a simple way to handle length differences without introducing complex alignment logic
        let reference_length = barcode_reference.len();
        let barcode_reference = if reference_length < extracted_length {
            let padding = "N".repeat(extracted_length - reference_length);
            format!("{}{}", barcode_reference, padding)
        } else if reference_length > extracted_length {
            let padding = "N".repeat(reference_length - extracted_length);
            format!("{}{}", barcode_extracted, padding)
        } else {
            barcode_reference.clone()
        };

        let distance = hamming(barcode_extracted.as_bytes(), barcode_reference.as_bytes()) as usize;
        if best_match.is_none() || distance < best_match.as_ref()?.distance {
            best_match = Some(BarcodeMatch {
                barcode: barcode_reference.clone(),
                distance,
            });
        }
    }

    // Right found our best match, now check if it meets the number of missmatches criteria
    // Get the length of the extracted barcode and the matched barcode, if they are not the same length we should not penalize
    // For example if we extract a 9bp barcode but the reference is 8bp we should not penalize for the extra base
    let barcode_extracted_length = barcode_extracted.len();
    let barcode_match_length = barcodes_reference[0].len();
    let length_diff = barcode_extracted_length.abs_diff(barcode_match_length);
    let max_distance = n_missmatches + length_diff;

    let best_match = best_match?;
    if best_match.distance > max_distance {
        None
    } else {
        Some((best_match.barcode, best_match.distance))
    }
}

/// Return a FASTQ reader and a writer based on file extensions.
/// Supports uncompressed `.fastq` / `.fq` and gzipped `.gz`.
pub fn get_reader_and_writer<P: AsRef<Path>>(
    input: P,
    output: P,
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
    path: P,
    barcodes: P,
    outfile: P,
    slack: Option<Slack>,
    n_missmatches: Option<usize>,
) -> Result<()> {
    let slack = slack.unwrap_or(Slack { left: 0, right: 0 });
    let n_missmatches = n_missmatches.unwrap_or(0);

    let (mut reader, mut writer) = get_reader_and_writer(path, outfile)?;

    info!("Loading barcodes from {:?}", barcodes.as_ref());
    let barcode_map = load_barcodes(barcodes)?;

    for (ii, result) in reader.records().enumerate() {
        let record = result?;
        let seq = std::str::from_utf8(record.sequence())?;
        let id = record.name().to_string();

        if ii % 1_000_000 == 0 {
            info!("Processed {ii} reads");
        }

        let mut barcodes_found: HashMap<BarcodeType, String> = HashMap::new();
        let mut barcodes_matched: HashMap<BarcodeType, (String, usize)> = HashMap::new();
        for bc_type in &[
            BarcodeType::BC1,
            BarcodeType::BC2,
            BarcodeType::BC3,
            BarcodeType::BC4,
        ] {
            let bc = extract_barcode(seq, bc_type, &slack)?;
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
        let new_id = format!(
            "{}__{}-{}-{}-{}",
            id,
            barcodes_matched
                .get(&BarcodeType::BC1)
                .map(|(bc, _)| bc.as_str())
                .unwrap_or("NNNNNNNN"),
            barcodes_matched
                .get(&BarcodeType::BC2)
                .map(|(bc, _)| bc.as_str())
                .unwrap_or("NNNNNNNN"),
            barcodes_matched
                .get(&BarcodeType::BC3)
                .map(|(bc, _)| bc.as_str())
                .unwrap_or("NNNNNNNN"),
            barcodes_matched
                .get(&BarcodeType::BC4)
                .map(|(bc, _)| bc.as_str())
                .unwrap_or("NNNNN"),
        );

        // Trim the sequence to remove barcodes
        // We have ~138 bp of barcodes on the 5' end of the read so just remove the first 138 bp
        let trimmed_seq = &seq[138..];
        let quality_scores = &record.quality_scores()[138..];
        let new_record = fastq::Record::new(
            Definition::new(new_id, record.description()),
            trimmed_seq,
            quality_scores,
        );
        writer.write_record(&new_record)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barcode_type_ranges() {
        assert_eq!(get_barcode_range(&BarcodeType::BC1), 76..84);
        assert_eq!(get_barcode_range(&BarcodeType::BC2), 38..46);
        assert_eq!(get_barcode_range(&BarcodeType::BC3), 0..8);
        assert_eq!(get_barcode_range(&BarcodeType::BC4), 114..119);
    }

    #[test]
    fn test_slack_creation() {
        let slack = Slack::new(2, 3);
        assert_eq!(slack.left, 2);
        assert_eq!(slack.right, 3);
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

        let slack = Slack::new(0, 0);

        let bc1 = extract_barcode(&sequence, &BarcodeType::BC1, &slack).unwrap();
        assert_eq!(bc1, "TTAATTAA");

        let bc2 = extract_barcode(&sequence, &BarcodeType::BC2, &slack).unwrap();
        assert_eq!(bc2, "GCTAGCTA");

        let bc3 = extract_barcode(&sequence, &BarcodeType::BC3, &slack).unwrap();
        assert_eq!(bc3, "ATCGATCG");

        let bc4 = extract_barcode(&sequence, &BarcodeType::BC4, &slack).unwrap();
        assert_eq!(bc4, "CCGGC");
    }

    #[test]
    fn test_extract_barcode_with_slack() {
        let sequence = "ATCGATCG" // BC3: 0..8
            .to_owned()
            + &"N".repeat(30)
            + "GCTAGCTA"
            + &"N".repeat(30)
            + "TTAATTAA"
            + &"N".repeat(30)
            + "CCGGC";

        let slack = Slack::new(1, 1);

        // With slack, BC3 should include one base before and after
        let bc3 = extract_barcode(&sequence, &BarcodeType::BC3, &slack);
        assert!(bc3.is_err()); // Should fail because slack.left > range.start (1 > 0)
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
