pub mod barcodes;
pub mod fragments;

use pyo3::prelude::*;

#[pymodule]
fn _scnado(_m: &Bound<'_, PyModule>) -> PyResult<()> {
    Ok(())
}
