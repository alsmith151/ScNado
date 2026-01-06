pub mod barcodes;
pub mod fragments;

#[cfg(feature = "extension-module")]
use pyo3::{prelude::*, wrap_pyfunction};

#[cfg(feature = "extension-module")]
mod python_bindings {
    use super::*;
    use pyo3::{PyResult, exceptions::PyRuntimeError};
    use std::path::PathBuf;

    #[pyfunction]
    pub fn find_barcodes(
        r1: String,
        r2: String,
        barcodes: String,
        output_prefix: String,
        n_missmatches: usize,
        enable_n_to_match: bool,
    ) -> PyResult<()> {
        let r1_path = PathBuf::from(r1);
        let r2_path = PathBuf::from(r2);
        let barcodes_path = PathBuf::from(barcodes);

        crate::barcodes::run(
            &r1_path,
            &r2_path,
            &barcodes_path,
            output_prefix,
            n_missmatches,
            enable_n_to_match,
            false,
        )
        .map_err(|e| PyRuntimeError::new_err(format!("Barcode finding failed: {}", e)))
    }

    #[pyfunction(signature = (bam, output, barcode_regex=None, shift_plus=0, shift_minus=0, fragment_length_extension=None))]
    pub fn fragments(
        bam: String,
        output: String,
        barcode_regex: Option<String>,
        shift_plus: i64,
        shift_minus: i64,
        fragment_length_extension: Option<usize>,
    ) -> PyResult<()> {
        let bam_path = PathBuf::from(bam);
        let output_path = PathBuf::from(output);

        let regex = if let Some(pattern) = barcode_regex {
            Some(
                regex::Regex::new(&pattern)
                    .map_err(|e| PyRuntimeError::new_err(format!("Invalid regex: {}", e)))?,
            )
        } else {
            None
        };

        crate::fragments::run(
            &bam_path,
            regex,
            &output_path,
            shift_plus,
            shift_minus,
            fragment_length_extension,
        )
        .map_err(|e| PyRuntimeError::new_err(format!("Fragment generation failed: {}", e)))
    }
}

#[cfg(feature = "extension-module")]
#[pymodule]
fn _scnado(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Register PyO3-exposed functions; signature annotation above handles optional params ordering.
    m.add_function(wrap_pyfunction!(python_bindings::find_barcodes, m)?)?;
    m.add_function(wrap_pyfunction!(python_bindings::fragments, m)?)?;
    // Suppress unused variable warning when `py` is not referenced directly.
    let _ = py;
    Ok(())
}
