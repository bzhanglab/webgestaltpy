use std::time::Instant;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use webgestalt_lib::methods::gsea::{GSEAConfig, GSEAResult};
use webgestalt_lib::methods::multiomics::{multiomic_ora, ORAJob};
use webgestalt_lib::methods::ora::{ORAConfig, ORAResult};

fn result_to_dict(obj: GSEAResult, py: Python) -> Result<&PyDict, PyErr> {
    let dict = PyDict::new(py);
    dict.set_item("set".to_object(py), obj.set.to_object(py))?;
    dict.set_item("p".to_object(py), obj.p.to_object(py))?;
    dict.set_item("fdr".to_object(py), obj.fdr.to_object(py))?;
    dict.set_item("es".to_object(py), obj.es.to_object(py))?;
    dict.set_item("nes".to_object(py), obj.nes.to_object(py))?;
    dict.set_item("leading_edge".to_object(py), obj.leading_edge.to_object(py))?;
    Ok(dict)
}

fn ora_result_to_dict(obj: ORAResult, py: Python) -> Result<&PyDict, PyErr> {
    let dict = PyDict::new(py);
    dict.set_item("set".to_object(py), obj.set.to_object(py))?;
    dict.set_item("p".to_object(py), obj.p.to_object(py))?;
    dict.set_item("fdr".to_object(py), obj.fdr.to_object(py))?;
    dict.set_item("overlap".to_object(py), obj.overlap.to_object(py))?;
    dict.set_item("expected".to_object(py), obj.expected.to_object(py))?;
    dict.set_item(
        "enrichment_ratio".to_object(py),
        obj.enrichment_ratio.to_object(py),
    )?;
    Ok(dict)
}

/// Run single-omic GSEA with files at provided paths.
///
/// # Parameters
/// - `gmt` - `String` of the path to the gmt file of interest
/// - `rank` - `String` of the path to the rank file of interest. Tab separated.
///
/// # Returns
///
/// Returns a list containing the GSEA results for every set.
///
/// # Panics
///
/// Panics if the GMT or the rank file is malformed or not at specified path.
///
/// # Example
///
/// ```python
/// import webgestaltpy
///
/// res = webgestaltpy.gsea_from_files("kegg.gmt", "example_ranked_list.rnk")
///
/// print(res[0:2]) # print first two results
/// ```
///
/// **Output**
///
/// ```
/// [
///   {
///     'set': 'hsa00010',
///     'p': 0.353,
///     'fdr': 1.276048474073356,
///     'es': 0.40653028852961814,
///     'nes': 1.07659486501464,
///     'leading_edge': 24
///   },
///   {
///     'set': 'hsa00020',
///     'p': 0,
///     'fdr': 0.028834551777982824,
///     'es': 0.6216527702210619,
///     'nes': 1.5721004858071521,
///     'leading_edge': 20
///   }
/// ]
/// ```
#[pyfunction]
fn gsea(py: Python, gmt: String, rank: String) -> PyResult<Vec<&PyDict>> {
    let analyte_list = webgestalt_lib::readers::read_rank_file(rank);
    let gmt = webgestalt_lib::readers::read_gmt_file(gmt);
    let start = Instant::now();
    let res: Vec<GSEAResult> = webgestalt_lib::methods::gsea::gsea(
        analyte_list.unwrap(),
        gmt.unwrap(),
        GSEAConfig::default(),
        None,
    );
    let new_res: Vec<&PyDict> = res
        .into_iter()
        .map(|x| result_to_dict(x, py).unwrap())
        .collect();
    let duration = start.elapsed();
    println!("GSEA\nTime took: {:?}", duration);
    Ok(new_res)
}

/// Run a single-omic ORA with files at the provided paths.
///
/// # Parameters
/// - `gmt_path` - `String` of the path to the gmt file of interest
/// - `analyte_list_path` - `String` of the path to the analyte file of interest.
/// - `reference_list_path`
///
/// # Returns
///
/// Returns a list of dictionaries with the results containing the ORA results for every set.
///
/// # Panics
///
/// Panics if the any file is malformed or not at specified path.
///
/// # Example
///
/// ```python
/// import webgestaltpy
///
/// res = webgestaltpy.ora_from_files("kegg.gmt", "gene_list.txt", "reference.txt")
///
/// print(res[0:2]) # print first two results
/// ```
///
/// **Output**
///
/// ```
/// [
///   {
///     'set': 'hsa00010',
///     'p': 0.7560574551180973,
///     'fdr': 1,
///     'overlap': 2,
///     'expected': 2.6840874707743088,
///     'enrichment_ratio': 0.7451321992211519
///   },
///   {
///     'set': 'hsa00020',
///     'p': 0.7019892669020903,
///     'fdr': 0.9981116297866582,
///     'overlap': 1,
///     'expected': 1.1841562371063128,
///     'enrichment_ratio': 0.8444831591173054
///   }
/// ]
/// ```
#[pyfunction]
fn ora(
    py: Python,
    gmt_path: String,
    analyte_list_path: String,
    reference_list_path: String,
) -> PyResult<Vec<&PyDict>> {
    let (gmt, analyte_list, reference) =
        webgestalt_lib::readers::read_ora_files(gmt_path, analyte_list_path, reference_list_path);
    let start = Instant::now();
    let res: Vec<ORAResult> =
        webgestalt_lib::methods::ora::get_ora(&analyte_list, &reference, gmt, ORAConfig::default());
    let new_res: Vec<&PyDict> = res
        .into_iter()
        .map(|x| ora_result_to_dict(x, py).unwrap())
        .collect();
    let duration = start.elapsed();
    println!("ORA\nTime took: {:?}", duration);
    Ok(new_res)
}

/// Run a meta-analysis ORA with files at the provided paths.
///
/// # Parameters
/// - `gmt_path` - `String` of the path to the gmt file of interest
/// - `analyte_list_paths` -  Lists of `String`s of the path to the analyte files of interest.
/// - `reference_list_paths` - Lists of `String`s of the paths to reference lists.
///
/// # Returns
///
/// Returns a list of dictionaries with the results containing the ORA results for every set.
///
/// # Panics
///
/// Panics if the any file is malformed or not at specified path.
///
/// # Example
///
/// ```python
/// import webgestaltpy
///
/// res = webgestaltpy.ora_from_files("kegg.gmt", "gene_list.txt", "reference.txt")
///
/// print(res[0:2]) # print first two results
/// ```
///
/// **Output**
///
/// ```
/// [
///   {
///     'set': 'hsa00010',
///     'p': 0.7560574551180973,
///     'fdr': 1,
///     'overlap': 2,
///     'expected': 2.6840874707743088,
///     'enrichment_ratio': 0.7451321992211519
///   },
///   {
///     'set': 'hsa00020',
///     'p': 0.7019892669020903,
///     'fdr': 0.9981116297866582,
///     'overlap': 1,
///     'expected': 1.1841562371063128,
///     'enrichment_ratio': 0.8444831591173054
///   }
/// ]
/// ```
#[pyfunction]
fn meta_ora(
    py: Python,
    gmt_path: String,
    analyte_list_paths: Vec<String>,
    reference_list_paths: Vec<String>,
) -> PyResult<Vec<Vec<&PyDict>>> {
    if analyte_list_paths.len() != reference_list_paths.len() {
        // Verify list sizes
        Err(PyValueError::new_err(format!(
            "Number of gene lists ({0}) and reference lists ({1}) don't match!",
            analyte_list_paths.len(),
            reference_list_paths.len()
        )))
    } else {
        let mut jobs: Vec<ORAJob> = Vec::new();
        for (i, analyte_list_path) in analyte_list_paths.iter().enumerate() {
            let (gmt, analyte_list, reference) = webgestalt_lib::readers::read_ora_files(
                gmt_path.clone(),
                analyte_list_path.clone(),
                reference_list_paths[i].clone(),
            );
            let new_job: ORAJob = ORAJob {
                gmt: gmt.clone(),
                interest_list: analyte_list.clone(),
                reference_list: reference.clone(),
                config: ORAConfig::default(),
            };
            jobs.push(new_job);
        }
        let rust_result = multiomic_ora(
            jobs,
            webgestalt_lib::methods::multiomics::MultiOmicsMethod::Meta(
                webgestalt_lib::methods::multiomics::MetaAnalysisMethod::Stouffer,
            ),
        );
        let mut final_results: Vec<Vec<&PyDict>> = Vec::new();
        for res in rust_result {
            let converted = res
                .into_iter()
                .map(|x| ora_result_to_dict(x, py).unwrap())
                .collect();
            final_results.push(converted);
        }
        Ok(final_results)
    }
}

/// High performance enrichment methods implemented in Rust, with Python bindings.
#[pymodule]
fn webgestaltpy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gsea, m)?)?;
    m.add_function(wrap_pyfunction!(ora, m)?)?;
    m.add_function(wrap_pyfunction!(meta_ora, m)?)?;
    Ok(())
}
