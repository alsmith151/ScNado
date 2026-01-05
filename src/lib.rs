pub mod barcodes;
pub mod fragments;

#[cfg(feature = "extension-module")]
use pyo3::prelude::*;

#[cfg(feature = "extension-module")]
#[pyfunction]
fn extract_barcode(sequence: String, _barcode_type: String) -> PyResult<String> {
    // This is a wrapper around the internal Rust function
    // We need to map the string to the BarcodeType enum
    // For now, let's just expose it for testing if needed,
    // but the main logic is in the CLI.
    Ok(sequence) // Placeholder
}

#[cfg(feature = "extension-module")]
#[pymodule]
fn _scnado(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(extract_barcode, m)?)?;
    Ok(())
}
