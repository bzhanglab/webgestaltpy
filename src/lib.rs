use std::time::Instant;

use pyo3::prelude::*;
use pyo3::types::PyDict;
use webgestalt_lib::methods::gsea::FullGSEAResult;

fn result_to_dict(obj: FullGSEAResult, py: Python) -> &PyDict {
    let res = vec![
        ("set".to_object(py), obj.set.to_object(py)),
        ("p".to_object(py), obj.p.to_object(py)),
        ("fdr".to_object(py), obj.fdr.to_object(py)),
        ("es".to_object(py), obj.es.to_object(py)),
        ("nes".to_object(py), obj.nes.to_object(py)),
        ("leading_edge".to_object(py), obj.leading_edge.to_object(py)),
    ]
    .to_object(py);
    PyDict::from_sequence(py, res).unwrap()
}

/// Run single-omic GSEA with files at provided paths.
///
/// # Parameters
/// - `gmt` - `String` of the path to the gmt file of interest
/// - `rank` - `String` of the path to the rank file of interest. Tab separated.
///
/// # Returns
///
/// Returns a [`PyResult<Vec<GSEAResultWrapper>>`] containing the GSEA results for every set.
///
/// # Panics
///
/// Panics if the GMT or the rank file is malformed.
#[pyfunction]
fn gsea_from_files(py: Python, gmt: String, rank: String) -> PyResult<Vec<&PyDict>> {
    let gene_list = webgestalt_lib::readers::read_rank_file(rank);
    let gmt = webgestalt_lib::readers::read_gmt_file(gmt);
    let start = Instant::now();
    let res: Vec<FullGSEAResult> =
        webgestalt_lib::methods::gsea::gsea(gene_list.unwrap(), gmt.unwrap());
    let new_res: Vec<&PyDict> = res.into_iter().map(|x| result_to_dict(x, py)).collect();
    let duration = start.elapsed();
    println!("New Hash\nTime took: {:?}", duration);
    Ok(new_res)
}

/// High performance enrichment methods implemented in Rust, with Python bindings.
#[pymodule]
fn webgestaltpy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gsea_from_files, m)?)?;
    Ok(())
}
