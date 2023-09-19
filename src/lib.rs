use std::time::Instant;

use pyo3::{exceptions::PyValueError, prelude::*};
use webgestalt_lib::methods::gsea::FullGSEAResult;
/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    if a == 4 {
        Err(PyValueError::new_err("EQUALS 4"))
    } else {
        Ok((a + b).to_string())
    }
}

#[pyclass]
#[pyo3(get_all, set_all, name = "GSEAResult")]
struct GSEAResultWrapper {
    pub set: String,
    pub p: f64,
    pub fdr: f64,
    pub es: f64,
    pub nes: f64,
    pub leading_edge: i32,
}

impl GSEAResultWrapper {
    fn from_lib(obj: FullGSEAResult) -> GSEAResultWrapper {
        let newobj = obj.clone();
        GSEAResultWrapper {
            set: newobj.set,
            p: newobj.p,
            fdr: newobj.fdr,
            es: newobj.es,
            nes: newobj.nes,
            leading_edge: newobj.leading_edge,
        }
    }
}

#[pyfunction]
fn gsea_from_files(gmt: String, rank: String) -> PyResult<Vec<GSEAResultWrapper>> {
    let gene_list = webgestalt_lib::readers::read_rank_file(rank);
    let gmt = webgestalt_lib::readers::read_gmt_file(gmt);
    let start = Instant::now();
    let res: Vec<FullGSEAResult> =
        webgestalt_lib::methods::gsea::gsea(gene_list.unwrap(), gmt.unwrap());
    let new_res: Vec<GSEAResultWrapper> = res
        .into_iter()
        .map(|x| GSEAResultWrapper::from_lib(x))
        .collect();
    let duration = start.elapsed();
    println!("New Hash\nTime took: {:?}", duration);
    Ok(new_res)
}

/// A Python module implemented in Rust.
#[pymodule]
fn webgestaltpy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(gsea_from_files, m)?)?;
    Ok(())
}
