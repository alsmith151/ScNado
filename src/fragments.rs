use anyhow::{Context, Result};
use log::{info, warn};
use noodles::bam;
use noodles::sam::alignment::record::data::field::tag;
use noodles::sam::alignment::record::data::field::value::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

fn get_barcode_from_tag(record: &bam::record::Record) -> Option<String> {
    let data = record.data();
    let tag = data.get(&tag::Tag::CELL_BARCODE_ID);

    match tag {
        Some(Ok(Value::String(s))) => Some(s.to_string()),
        Some(Ok(Value::Hex(h))) => Some(String::from_utf8_lossy(h).to_string()),
        _ => None,
    }
}

fn get_barcode_from_read_name(
    record: &bam::record::Record,
    regex: &regex::Regex,
) -> Option<String> {
    let read_name = record.name()?.to_string();
    let caps = regex.captures(&read_name)?;
    caps.get(1).map(|m| m.as_str().to_string())
}

pub fn run<P: AsRef<Path>>(
    bam_file: P,
    barcode_extraction_regex: Option<regex::Regex>,
    outfile: P,
    shift_plus: i64,
    shift_minus: i64,
    fragment_length_extension: Option<usize>,
) -> Result<()> {
    info!("Processing BAM file: {:?}", bam_file.as_ref());
    let mut reader = bam::io::reader::Builder.build_from_path(bam_file)?;

    let header = reader.read_header()?;
    let reference_sequence = header.reference_sequences();

    // Create output file direcrtory if it does not exist
    if let Some(parent) = outfile.as_ref().parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut writer = std::io::BufWriter::new(File::create(outfile)?);
    let mut fragments = HashMap::new();

    for result in reader.records() {
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }
        let barcode = if let Some(regex) = &barcode_extraction_regex {
            get_barcode_from_read_name(&record, regex)
        } else {
            get_barcode_from_tag(&record)
        };

        if barcode.is_none() {
            warn!(
                "Record {} has no barcode",
                record.name().unwrap_or_default()
            );
            continue;
        }
        let barcode = barcode.unwrap();

        let ref_id = record
            .reference_sequence_id()
            .context("No reference sequence ID")??;
        let chrom = reference_sequence
            .get_index(ref_id)
            .map(|(_name, _len)| String::from_utf8_lossy(_name).to_string())
            .unwrap_or_else(|| "chr?".to_string());

        // Get 0 based start position
        let start = record
            .alignment_start()
            .context("No alignment start")??
            .get()
            - 1;
        let end = record.sequence().len() + start - 1;

        // compute insertion coordinate with Tn5 shift
        let insertion_pos: i64 = if record.flags().is_reverse_complemented() {
            end as i64 - shift_minus
        } else {
            start as i64 + shift_plus
        };

        // Clip to non-negative coordinates
        let insertion_pos = std::cmp::max(0, insertion_pos);

        // Define fragment start and end positions
        let (frag_start, frag_end) = if let Some(length) = fragment_length_extension {
            let (s, e) = if record.flags().is_reverse_complemented() {
                (insertion_pos - length as i64 + 1, insertion_pos + 1)
            } else {
                (insertion_pos, insertion_pos + length as i64)
            };

            (s, e)
        } else {
            (insertion_pos, insertion_pos + 1)
        };

        // Check that insertion_pos is within chromosome bounds
        let chrom_length = reference_sequence
            .get_index(ref_id)
            .map(|(_name, len)| len.length().into())
            .unwrap_or(0);
        if frag_end >= chrom_length as i64 {
            warn!(
                "Insertion position {} is out of bounds for chromosome {} with length {}",
                frag_end, chrom, chrom_length
            );
            continue;
        }

        fragments
            .entry((
                chrom.clone(),
                frag_start as usize,
                frag_end as usize,
                barcode.clone(),
            ))
            .and_modify(|c| *c += 1)
            .or_insert(1);
    }

    for ((chrom, start, end, barcode), count) in fragments.into_iter() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            chrom, start, end, barcode, count
        )?;
    }

    Ok(())
}
