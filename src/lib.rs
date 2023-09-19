use std::time::Instant;

use pyo3::{exceptions::PyValueError, prelude::*};
/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    if a == 4 {
        Err(PyValueError::new_err("EQUALS 4"))
    } else {
        Ok((a + b).to_string())
    }
}

#[pyfunction]
fn gsea_from_files(gmt: String, rank: String) -> PyResult<bool> {
    let gene_list = webgestalt_lib::readers::read_rank_file(rank);
    let gmt = webgestalt_lib::readers::read_gmt_file(gmt);
    let start = Instant::now();
    webgestalt_lib::methods::gsea::gsea(gene_list.unwrap(), gmt.unwrap());
    let duration = start.elapsed();
    println!("New Hash\nTime took: {:?}", duration);
    Ok(true)
}

/// A Python module implemented in Rust.
#[pymodule]
fn webgestaltpy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(gsea_from_files, m)?)?;
    Ok(())
}

