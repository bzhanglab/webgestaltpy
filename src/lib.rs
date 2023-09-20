use std::time::Instant;

use pyo3::prelude::*;
use webgestalt_lib::methods::gsea::FullGSEAResult;

#[pyclass(mapping, dict)]
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
    /// Clone then convert [`webgestalt_lib::methods::gsea::FullGSEAResult`] to [`GSEAResultWrapper`].
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
fn gsea_from_files(gmt: String, rank: String) -> PyResult<Vec<GSEAResultWrapper>> {
    let gene_list = webgestalt_lib::readers::read_rank_file(rank);
    let gmt = webgestalt_lib::readers::read_gmt_file(gmt);
    let start = Instant::now();
    let res: Vec<FullGSEAResult> =
        webgestalt_lib::methods::gsea::gsea(gene_list.unwrap(), gmt.unwrap());
    let new_res: Vec<GSEAResultWrapper> = res
        .into_iter()
        .map(GSEAResultWrapper::from_lib)
        .collect();
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
