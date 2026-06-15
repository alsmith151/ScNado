use anyhow::{Context, Result};
use log::{info, warn};
use noodles::bam;
use noodles::sam::alignment::record::data::field::tag;
use noodles::sam::alignment::record::data::field::value::Value;
use std::collections::{HashMap, HashSet};
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

fn extract_barcode_and_umi(
    record: &bam::record::Record,
    barcode_regex: Option<&regex::Regex>,
    umi_regex: Option<&regex::Regex>,
) -> (Option<String>, Option<String>) {
    let read_name = record.name().map(|n| n.to_string());

    let mut barcode = None;
    let mut umi = None;

    if let Some(rx) = barcode_regex
        && let Some(name) = &read_name
        && let Some(caps) = rx.captures(name)
    {
        barcode = caps
            .name("barcode")
            .map(|m| m.as_str().to_string())
            .or_else(|| caps.get(1).map(|m| m.as_str().to_string()));

        if umi_regex.is_none() {
            umi = caps
                .name("UMI")
                .or(caps.name("umi"))
                .map(|m| m.as_str().to_string());
        }
    }

    if let Some(rx) = umi_regex
        && let Some(name) = &read_name
        && let Some(caps) = rx.captures(name)
    {
        umi = caps
            .name("UMI")
            .or(caps.name("umi"))
            .map(|m| m.as_str().to_string())
            .or_else(|| caps.get(1).map(|m| m.as_str().to_string()));
    }

    if barcode.is_none() {
        barcode = get_barcode_from_tag(record);
    }

    (barcode, umi)
}

struct FragmentCounter {
    simple_counts: HashMap<(String, usize, usize, String), u32>,
    umi_counts: HashMap<(String, usize, usize, String), HashSet<String>>,
    use_umi: bool,
}

impl FragmentCounter {
    fn new(use_umi: bool) -> Self {
        Self {
            simple_counts: HashMap::new(),
            umi_counts: HashMap::new(),
            use_umi,
        }
    }

    fn add(&mut self, key: (String, usize, usize, String), umi: Option<String>) {
        if self.use_umi {
            // If UMI is missing, treat it as empty string (collapsing all missing-UMI reads at this locus)
            let u = umi.unwrap_or_default();
            self.umi_counts.entry(key).or_default().insert(u);
        } else {
            *self.simple_counts.entry(key).or_default() += 1;
        }
    }

    fn write_results<W: Write>(&self, writer: &mut W) -> Result<()> {
        if self.use_umi {
            for ((chrom, start, end, barcode), umis) in &self.umi_counts {
                if !umis.is_empty() {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}",
                        chrom,
                        start,
                        end,
                        barcode,
                        umis.len()
                    )?;
                }
            }
        } else {
            for ((chrom, start, end, barcode), count) in &self.simple_counts {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    chrom, start, end, barcode, count
                )?;
            }
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run<P: AsRef<Path>>(
    bam_file: P,
    barcode_regex: Option<regex::Regex>,
    umi_regex: Option<regex::Regex>,
    outfile: P,
    shift_plus: i64,
    shift_minus: i64,
    fragment_length_extension: Option<usize>,
    check_cancel: Option<&dyn Fn() -> Result<()>>,
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

    // Determine use_umi
    let use_umi = umi_regex.is_some()
        || barcode_regex.as_ref().is_some_and(|r| {
            r.capture_names()
                .any(|n| n == Some("UMI") || n == Some("umi"))
        });

    let mut counter = FragmentCounter::new(use_umi);

    for (i, result) in reader.records().enumerate() {
        if let Some(check) = check_cancel
            && i % 10000 == 0
        {
            check()?;
        }
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }

        let (barcode, umi) =
            extract_barcode_and_umi(&record, barcode_regex.as_ref(), umi_regex.as_ref());

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

        counter.add(
            (chrom, frag_start as usize, frag_end as usize, barcode),
            umi,
        );
    }

    counter.write_results(&mut writer)?;

    Ok(())
}
