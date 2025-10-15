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

/// Identify the best matching barcode from a list of valid barcodes using Hamming distance.
/// Returns the best matching barcode and its Hamming distance if found, otherwise None.
/// # Arguments
/// * `extracted` - The extracted barcode sequence to be matched.
/// * `valid_barcodes` - A slice of valid barcode sequences to match against.
/// * `minimum_score` - The minimum Hamming distance score to consider a match valid.
/// # Returns
/// An Option containing a tuple of the best matching barcode and its Hamming distance, or None if no match is found.
fn identify_best_barcode(
    extracted: &str,
    valid_barcodes: &[String],
    n_missmatches: usize,
) -> Option<(String, usize)> {
    let mut best_match: Option<(String, usize)> = None;
    let mut min_distance = usize::MAX;

    for valid in valid_barcodes {
        let distance = hamming(extracted.as_bytes(), valid.as_bytes()) as usize;
        if distance < min_distance {
            min_distance = distance;
            best_match = Some((valid.clone(), distance));
        }
    }

    let max_score = valid_barcodes[0].len();
    let min_score = max_score.saturating_sub(n_missmatches);

    if let Some((_, distance)) = &best_match
        && *distance > min_score {
            return None;
        }
    best_match
}

/// Return a FASTQ reader and a writer based on file extensions.
/// Supports uncompressed `.fastq` / `.fq` and gzipped `.gz`.
pub fn get_reader_and_writer<P: AsRef<Path>>(
    input: P,
    output: P,
) -> Result<(
    fastq::io::Reader<Box<dyn std::io::BufRead>>,
    fastq::io::Writer<Box<dyn std::io::Write>>,
)> {
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

    
    for (ii, result) in reader.records().enumerate()  {
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
