#[pyo3::pymodule]
mod webgestaltpy {
    use ahash::AHashSet;
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;
    use pyo3::types::PyDict;
    use webgestalt_lib::methods::gsea::{GSEAConfig, GSEAResult, RankListItem};
    use webgestalt_lib::methods::multilist::{multilist_gsea, multilist_ora, GSEAJob, ORAJob};
    use webgestalt_lib::methods::nta::{NTAConfig, NTAResult};
    use webgestalt_lib::methods::ora::{ORAConfig, ORAResult};
    use webgestalt_lib::readers::utils::Item;

    /// Enum of the NTA Methods supported by WebGestalt
    ///
    /// # Enum Values
    ///
    /// - `Prioritization` - Finds the N seeds (input analytes) that are most likely to be encountered with a random walk
    /// - `Expansion` - Finds the N non-seed (non-input analytes) nodes that are most likely to be encountered with a random walk
    ///
    /// # Example
    ///
    /// ```python
    /// import webgestaltpy
    ///
    /// method = webgestaltpy.NTAMethod.Expansion
    /// ```
    #[pyclass]
    pub enum NTAMethod {
        /// Finds the N seeds (input analytes) that are most likely to be encountered with a random walk
        Prioritization,
        /// Finds the N non-seed (non-input analytes) nodes that are most likely to be encountered with a random walk
        Expansion,
    }

    fn gsea_result_to_dict<'a>(
        obj: GSEAResult,
        py: Python<'a>,
    ) -> Result<pyo3::Bound<'a, PyDict>, PyErr> {
        let dict = PyDict::new(py);
        dict.set_item("set", obj.set)?;
        dict.set_item("p", obj.p)?;
        dict.set_item("fdr", obj.fdr)?;
        dict.set_item("es", obj.es)?;
        dict.set_item("nes", obj.nes)?;
        dict.set_item("leading_edge", obj.leading_edge)?;
        Ok(dict)
    }

    fn ora_result_to_dict<'a>(
        obj: ORAResult,
        py: Python<'a>,
    ) -> Result<pyo3::Bound<'a, PyDict>, PyErr> {
        let dict = PyDict::new(py);
        dict.set_item("set", obj.set)?;
        dict.set_item("p", obj.p)?;
        dict.set_item("fdr", obj.fdr)?;
        dict.set_item("overlap", obj.overlap)?;
        dict.set_item("expected", obj.expected)?;
        dict.set_item("enrichment_ratio", obj.enrichment_ratio)?;
        Ok(dict)
    }

    fn nta_result_to_dict<'a>(
        obj: NTAResult,
        py: Python<'a>,
    ) -> Result<pyo3::Bound<'a, PyDict>, PyErr> {
        let dict = PyDict::new(py);
        dict.set_item("candidates", obj.candidates)?;
        dict.set_item("scores", obj.scores)?;
        dict.set_item("neighborhood", obj.neighborhood)?;
        Ok(dict)
    }

    /// Run single-omic NTA (Network-topology based analysis) with files at the provided paths
    ///
    /// # Parameters
    /// - `edge_list_path` - `String` of the path to the edge list file of the network. See below for details.
    /// - `analyte_list_path` - `String` of the path to the seed nodes, with entries separated by new lines
    /// - `nta_method` - a `NTAMethod` object specifying the NTA method for the analysis.
    /// - `n` - the number of seeds or nodes to identify according to `nta_method`
    ///
    /// # Returns
    ///
    /// Returns a dictionary object containing the `candidates` (seed nodes when using prioritization), `scores` (random-walk probabilities), and `neighborhood` (identified nodes)
    ///
    /// # Panics
    ///
    /// Panics if the network or the analyte file is malformed or not at specified path. Will also panic if `nta_method` is not specified correctly
    ///
    /// # Example
    ///
    /// ```python
    /// import webgestaltpy
    ///
    /// nta_method = webgestaltpy.NTAMethod.Prioritization
    /// y = webgestaltpy.nta_from_files("data/hsapiens_network_CPTAC_Proteomics_OV_entrezgene.net", "data/net_genes.txt", nta_method, 5)
    /// print(y)
    /// ```
    ///
    /// **Output**
    ///
    /// ```
    /// {
    ///   'candidates': [
    ///     'ACTA1',
    ///     'ACTA2',
    ///     'ACTB',
    ///     'ACTG1'
    ///   ],
    ///   'scores': [
    ///     0.015611545101449542,
    ///     0.015611545101449542,
    ///     0.015227515228472441,
    ///     0.015227515228472441,
    ///     0.015105514420304793
    ///   ],
    ///   'neighborhood': [
    ///     'ACTA1',
    ///     'ACTA2',
    ///     'ACTB',
    ///     'ACTG1',
    ///     'ACTG2'
    ///   ]
    /// }
    /// ```
    #[pyfunction]
    fn nta_from_files<'a>(
        py: Python<'a>,
        edge_list_path: String,
        analyte_list_path: String,
        nta_method: &'a NTAMethod,
        n: usize,
    ) -> PyResult<pyo3::Bound<'a, PyDict>> {
        let net_file = webgestalt_lib::readers::read_edge_list(edge_list_path);
        let analytes = webgestalt_lib::readers::read_single_list(analyte_list_path);
        let method = match nta_method {
            NTAMethod::Expansion => webgestalt_lib::methods::nta::NTAMethod::Expand(n),
            NTAMethod::Prioritization => webgestalt_lib::methods::nta::NTAMethod::Prioritize(n),
        };
        let res = webgestalt_lib::methods::nta::get_nta(NTAConfig {
            edge_list: net_file,
            seeds: analytes.into_iter().collect(),
            method: Option::Some(method),
            ..Default::default()
        });
        let new_res = nta_result_to_dict(res, py)?;
        Ok(new_res)
    }

    /// Run single-omic NTA (Network-topology based analysis) with edge list and list of starting analytes/nodes
    ///
    /// # Parameters
    /// - `edge_list` - `list[list[str]]` of the network which is a list of lists where each entry is a node.
    /// - `analyte_list` - `list[str]` of analytes for starting the NTA with.
    /// - `nta_method` - a `NTAMethod` object specifying the NTA method for the analysis.
    /// - `n` - the number of seeds or nodes to identify according to `nta_method`
    ///
    /// # Returns
    ///
    /// Returns a dictionary object containing the `candidates` (seed nodes when using prioritization), `scores` (random-walk probabilities), and `neighborhood` (identified nodes)
    ///
    /// # Panics
    ///
    /// Panics if the network or the analyte file is malformed or not at specified path. Will also panic if `nta_method` is not specified correctly
    ///
    /// # Example
    ///
    /// ```python
    /// import webgestaltpy
    ///
    /// def file_to_list(file_path: str) -> list[str]:
    ///     """Convert file to list[str]"""
    ///     with open(file_path, "r") as r:
    ///         return list(r.readlines())
    ///
    /// def read_network(file_path: str) -> list[list[str]]:
    ///     """Read .net file to format for nta()"""
    ///     with open(file_path, "r") as r:
    ///         lines = r.readlines()
    ///     net = []
    ///     for line in lines:
    ///         if "\t" in line:
    ///             net.append(line.split("\t"))
    ///     return net
    ///
    /// nta_method = webgestaltpy.NTAMethod.Prioritization
    /// y = webgestaltpy.nta(read_network("data/hsapiens_network_CPTAC_Proteomics_OV_entrezgene.net"),
    ///                      file_to_list("data/net_genes.txt"(), nta_method, 5)
    /// print(y)
    /// ```
    ///
    /// **Output**
    ///
    /// ```
    /// {
    ///   'candidates': [
    ///     'ACTA1',
    ///     'ACTA2',
    ///     'ACTB',
    ///     'ACTG1'
    ///   ],
    ///   'scores': [
    ///     0.015611545101449542,
    ///     0.015611545101449542,
    ///     0.015227515228472441,
    ///     0.015227515228472441,
    ///     0.015105514420304793
    ///   ],
    ///   'neighborhood': [
    ///     'ACTA1',
    ///     'ACTA2',
    ///     'ACTB',
    ///     'ACTG1',
    ///     'ACTG2'
    ///   ]
    /// }
    /// ```
    #[pyfunction]
    fn nta<'a>(
        py: Python<'a>,
        edge_list: Vec<Vec<String>>,
        analyte_list: Vec<String>,
        nta_method: &'a NTAMethod,
        n: usize,
    ) -> PyResult<pyo3::Bound<'a, PyDict>> {
        let analytes: AHashSet<String> = analyte_list.into_iter().collect();
        let method = match nta_method {
            NTAMethod::Expansion => webgestalt_lib::methods::nta::NTAMethod::Expand(n),
            NTAMethod::Prioritization => webgestalt_lib::methods::nta::NTAMethod::Prioritize(n),
        };
        let res = webgestalt_lib::methods::nta::get_nta(NTAConfig {
            edge_list,
            seeds: analytes.into_iter().collect(),
            method: Option::Some(method),
            ..Default::default()
        });
        let new_res = nta_result_to_dict(res, py)?;
        Ok(new_res)
    }

    /// Run single-omic GSEA with files at provided paths.
    ///
    /// # Parameters
    /// - `gmt_path` - `String` of the path to the gmt file of interest
    /// - `rank_file_path` - `String` of the path to the rank file of interest. Tab separated.
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
    ///     'set': 'has00010',
    ///     'p': 0.353,
    ///     'fdr': 1.0,
    ///     'es': 0.40653028852961814,
    ///     'nes': 1.07659486501464,
    ///     'leading_edge': 24
    ///   },
    ///   {
    ///     'set': 'has00020',
    ///     'p': 0,
    ///     'fdr': 0.028834551777982824,
    ///     'es': 0.6216527702210619,
    ///     'nes': 1.5721004858071521,
    ///     'leading_edge': 20
    ///   }
    /// ]
    /// ```
    #[pyfunction]
    fn gsea_from_files<'a>(
        py: Python<'a>,
        gmt_path: String,
        rank_file_path: String,
    ) -> PyResult<Vec<pyo3::Bound<'a, PyDict>>> {
        let analyte_list = webgestalt_lib::readers::read_rank_file(rank_file_path);
        let gmt = webgestalt_lib::readers::read_gmt_file(gmt_path);
        let res: Vec<GSEAResult> = webgestalt_lib::methods::gsea::gsea(
            analyte_list.unwrap(),
            gmt.unwrap(),
            GSEAConfig::default(),
            None,
        );
        let new_res: Vec<pyo3::Bound<PyDict>> = res
            .into_iter()
            .map(|x| gsea_result_to_dict(x, py).unwrap())
            .collect();
        Ok(new_res)
    }

    /// Run single-omic GSEA with a gmt and a rank list in form of `list[tuple[str, float]]`.
    ///
    /// # Parameters
    /// - `gmt_path` - `String` of the path to the gmt file of interest
    /// - `rank_list` - `list[tuple[str, float]]` of the path to the rank file of interest. Tab separated.
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
    /// def create_rank_list(file_path: str) -> list[tuple[str, float]]:
    ///     """Function that converts rnk file to format for gsea().
    ///
    ///     This is for demo purposes. gsea_from_files would be more convenient in this case.
    ///     However, for when the scores are calculated in a Python script, gsea() allows you to directly
    ///     use the results without having to save to a file first.
    ///
    ///     """
    ///
    ///     res = []
    ///     with open(file_path, "r") as r:
    ///         lines = r.readlines()
    ///     for line in lines:
    ///         if "\t" in line:
    ///             vals = line.split("\t")
    ///             res.append((vals[0], float(vals[1])))
    ///     return res
    ///
    /// res = webgestaltpy.gsea("kegg.gmt", create_rank_list("example_ranked_list.rnk"))
    ///
    /// print(res[0:2]) # print first two results
    /// ```
    ///
    /// **Output**
    ///
    /// _Your results may vary depending on random permutations_
    ///
    /// ```python
    /// [
    ///   {
    ///     'set': 'has00010',
    ///     'p': 0.353,
    ///     'fdr': 1.0,
    ///     'es': 0.40653028852961814,
    ///     'nes': 1.07659486501464,
    ///     'leading_edge': 24
    ///   },
    ///   {
    ///     'set': 'has00020',
    ///     'p': 0,
    ///     'fdr': 0.028834551777982824,
    ///     'es': 0.6216527702210619,
    ///     'nes': 1.5721004858071521,
    ///     'leading_edge': 20
    ///   }
    /// ]
    /// ```
    #[pyfunction]
    fn gsea<'a>(
        py: Python<'a>,
        gmt_path: String,
        rank_file: Vec<(String, f64)>,
    ) -> PyResult<Vec<pyo3::Bound<'a, PyDict>>> {
        let analyte_list = rank_file
            .iter()
            .map(|(analyte, value)| RankListItem {
                analyte: analyte.clone(),
                rank: *value,
            })
            .collect();
        let gmt = webgestalt_lib::readers::read_gmt_file(gmt_path);
        let res: Vec<GSEAResult> = webgestalt_lib::methods::gsea::gsea(
            analyte_list,
            gmt.unwrap(),
            GSEAConfig::default(),
            None,
        );
        let new_res: Vec<pyo3::Bound<PyDict>> = res
            .into_iter()
            .map(|x| gsea_result_to_dict(x, py).unwrap())
            .collect();
        Ok(new_res)
    }

    /// Run a meta-analysis GSEA with files at the provided paths.
    ///
    /// # Parameters
    /// - `gmt_path` - `String` of the path to the gmt file of interest
    /// - `rank_files` -  Lists of `String`s of the paths to the rank files of interest. Tab separated.
    ///
    /// # Returns
    ///
    /// Returns a list of a list of dictionaries with the results containing the GSEA results for every set.
    ///
    /// The first list contains the results of the meta-analysis. The following lists are the results for each list individually.
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
    /// res = webgestaltpy.meta_gsea_from_files("kegg.gmt", ["rank_list1.txt", "rank_list2.txt"])
    /// ```
    ///
    /// `res` would be a list containing the results of the meta-analysis and each list run
    /// individually. In this example, `res[0]` would look be the results of the meta-analysis.
    /// `res[1]` would be the results from `rank_list1.txt`, `res[2]` would be the results from `rank_list2.txt`, and so on.
    ///
    /// See the documentation for [`webgestaltpy.gsea`](./gsea.md) for specifics about the format of the results.
    #[pyfunction]
    fn meta_gsea_from_files<'a>(
        py: Python<'a>,
        gmt: String,
        rank_files: Vec<String>,
    ) -> PyResult<Vec<Vec<pyo3::Bound<'a, PyDict>>>> {
        let mut jobs: Vec<GSEAJob> = Vec::new();
        let gmt_vec: Vec<Item> = webgestalt_lib::readers::read_gmt_file(gmt).unwrap();
        for rank_file in rank_files {
            let analyte_list_result = webgestalt_lib::readers::read_rank_file(rank_file.clone());
            if analyte_list_result.is_ok() {
                let analyte_list = analyte_list_result.unwrap();
                let new_job = GSEAJob {
                    gmt: gmt_vec.clone(),
                    rank_list: analyte_list.clone(),
                    config: GSEAConfig::default(),
                };
                jobs.push(new_job);
            } else {
                return Err(PyValueError::new_err(format!(
                    "Error when reading rank file at: {}",
                    rank_file
                )));
            }
        }

        let rust_result = multilist_gsea(
            jobs,
            webgestalt_lib::methods::multilist::MultiListMethod::Meta(
                webgestalt_lib::methods::multilist::MetaAnalysisMethod::Stouffer,
            ),
            webgestalt_lib::stat::AdjustmentMethod::BH,
        );
        let mut final_results: Vec<Vec<pyo3::Bound<PyDict>>> = Vec::new();
        for res in rust_result {
            let converted: Vec<pyo3::Bound<PyDict>> = res
                .into_iter()
                .map(|x| gsea_result_to_dict(x, py).unwrap())
                .collect();
            final_results.push(converted);
        }
        Ok(final_results)
    }

    /// Run a meta-analysis GSEA with provided rank files.
    ///
    /// # Parameters
    /// - `gmt_path` - `String` of the path to the gmt file of interest
    /// - `rank_lists` -  Lists of `list[tuple[str, float]]`s of the rank lists of interest.
    ///
    /// # Returns
    ///
    /// Returns a list of a list of dictionaries with the results containing the GSEA results for every set.
    ///
    /// The first list contains the results of the meta-analysis. The following lists are the results for each list individually.
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
    /// def create_rank_list(file_path: str) -> list[tuple[str, float]]:
    ///     """Function that converts rnk file to format for gsea().
    ///
    ///     This is for demo purposes. gsea_from_files would be more convenient in this case.
    ///     However, for when the scores are calculated in a Python script, gsea() allows you to directly
    ///     use the results without having to save to a file first.
    ///
    ///     """
    ///
    ///     res = []
    ///     with open(file_path, "r") as r:
    ///         lines = r.readlines()
    ///     for line in lines:
    ///         if "\t" in line:
    ///             vals = line.split("\t")
    ///             res.append((vals[0], float(vals[1])))
    ///     return re/
    /// res = webgestaltpy.meta_gsea("kegg.gmt", [create_rank_list("rank_list1.txt"), create_rank_list("rank_list2.txt")])
    /// ```
    ///
    /// `res` would be a list containing the results of the meta-analysis and each list run
    /// individually. In this example, `res[0]` would look be the results of the meta-analysis.
    /// `res[1]` would be the results from `rank_list1.txt`, `res[2]` would be the results from `rank_list2.txt`, and so on.
    ///
    /// See the documentation for [`webgestaltpy.gsea`](./gsea.md) for specifics about the format of the results.
    #[pyfunction]
    fn meta_gsea<'a>(
        py: Python<'a>,
        gmt: String,
        rank_lists: Vec<Vec<(String, f64)>>,
    ) -> PyResult<Vec<Vec<pyo3::Bound<'a, PyDict>>>> {
        let mut jobs: Vec<GSEAJob> = Vec::new();
        let gmt_vec: Vec<Item> = webgestalt_lib::readers::read_gmt_file(gmt).unwrap();
        for rank_file in rank_lists {
            let analyte_list = rank_file
                .iter()
                .map(|(analyte, value)| RankListItem {
                    analyte: analyte.clone(),
                    rank: *value,
                })
                .collect();
            let new_job = GSEAJob {
                gmt: gmt_vec.clone(),
                rank_list: analyte_list,
                config: GSEAConfig::default(),
            };
            jobs.push(new_job);
        }

        let rust_result = multilist_gsea(
            jobs,
            webgestalt_lib::methods::multilist::MultiListMethod::Meta(
                webgestalt_lib::methods::multilist::MetaAnalysisMethod::Stouffer,
            ),
            webgestalt_lib::stat::AdjustmentMethod::BH,
        );
        let mut final_results: Vec<Vec<pyo3::Bound<PyDict>>> = Vec::new();
        for res in rust_result {
            let converted = res
                .into_iter()
                .map(|x| gsea_result_to_dict(x, py).unwrap())
                .collect();
            final_results.push(converted);
        }
        Ok(final_results)
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
    ///     'set': 'has00010',
    ///     'p': 0.7560574551180973,
    ///     'fdr': 1,
    ///     'overlap': 2,
    ///     'expected': 2.6840874707743088,
    ///     'enrichment_ratio': 0.7451321992211519
    ///   },
    ///   {
    ///     'set': 'has00020',
    ///     'p': 0.7019892669020903,
    ///     'fdr': 0.9981116297866582,
    ///     'overlap': 1,
    ///     'expected': 1.1841562371063128,
    ///     'enrichment_ratio': 0.8444831591173054
    ///   }
    /// ]
    /// ```
    #[pyfunction]
    fn ora_from_files<'a>(
        py: Python<'a>,
        gmt_path: String,
        analyte_list_path: String,
        reference_list_path: String,
    ) -> PyResult<Vec<pyo3::Bound<'a, PyDict>>> {
        let (gmt, analyte_list, reference) = webgestalt_lib::readers::read_ora_files(
            gmt_path,
            analyte_list_path,
            reference_list_path,
        );
        let res: Vec<ORAResult> = webgestalt_lib::methods::ora::get_ora(
            &analyte_list,
            &reference,
            gmt,
            ORAConfig::default(),
        );
        let new_res: Vec<pyo3::Bound<PyDict>> = res
            .into_iter()
            .map(|x| ora_result_to_dict(x, py).unwrap())
            .collect();
        Ok(new_res)
    }

    /// Run a single-omic ORA with gmt file and list of analytes and reference provided as strings.
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
    /// res = webgestaltpy.ora("kegg.gmt", gene_list, reference) # gene_list and reference are both list[str]
    ///
    /// print(res[0:2]) # print first two results
    /// ```
    ///
    /// **Output**
    ///
    /// ```
    /// [
    ///   {
    ///     'set': 'has00010',
    ///     'p': 0.7560574551180973,
    ///     'fdr': 1,
    ///     'overlap': 2,
    ///     'expected': 2.6840874707743088,
    ///     'enrichment_ratio': 0.7451321992211519
    ///   },
    ///   {
    ///     'set': 'has00020',
    ///     'p': 0.7019892669020903,
    ///     'fdr': 0.9981116297866582,
    ///     'overlap': 1,
    ///     'expected': 1.1841562371063128,
    ///     'enrichment_ratio': 0.8444831591173054
    ///   }
    /// ]
    /// ```
    #[pyfunction]
    fn ora<'a>(
        py: Python<'a>,
        gmt_path: String,
        analyte_list: Vec<String>,
        reference_list: Vec<String>,
    ) -> PyResult<Vec<pyo3::Bound<'a, PyDict>>> {
        let gmt = webgestalt_lib::readers::read_gmt_file(gmt_path).unwrap();
        let reference: AHashSet<String> = reference_list.into_iter().collect();
        let analyte_list: AHashSet<String> = analyte_list.into_iter().collect();
        let res: Vec<ORAResult> = webgestalt_lib::methods::ora::get_ora(
            &analyte_list,
            &reference,
            gmt,
            ORAConfig::default(),
        );
        let new_res: Vec<pyo3::Bound<PyDict>> = res
            .into_iter()
            .map(|x| ora_result_to_dict(x, py).unwrap())
            .collect();
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
    /// Returns a list of a list of dictionaries with the results containing the ORA results for every set.
    ///
    /// The first list contains the results of the meta-analysis. The following lists are the results for each list individually.
    ///
    /// # Panics
    ///
    /// Panics if the any file is malformed or not at specified path.
    ///
    /// # Example
    ///
    ///```python
    /// import webgestaltpy
    ///
    /// def file_to_list(file_path: str) -> list[str]:
    ///     """Convert file to list[str]"""
    ///     with open(file_path, "r") as r:
    ///         return list(r.readlines())
    ///
    /// res = webgestaltpy.meta_ora(
    ///     "kegg.gmt",
    ///     [file_to_list("gene_list1.txt"), file_to_list("gene_list2.txt")],
    ///     [file_to_list("reference.txt"), file_to_list("reference.txt")],
    /// )
    /// ```
    ///
    /// `res` would be a list containing the results of the meta-analysis and each list run
    /// individually. In this example, `res[0]` would look be the results of the meta-analysis.
    /// `res[1]` would be the results from `gene_list1.txt`, `res[2]` would be the results from `gene_list2.txt`, and so on.
    ///
    /// See the documentation for [`webgestaltpy.ora`](./ora.md) for specifics about the format of the results.
    #[pyfunction]
    fn meta_ora<'a>(
        py: Python<'a>,
        gmt_path: String,
        analyte_lists: Vec<Vec<String>>,
        reference_lists: Vec<Vec<String>>,
    ) -> PyResult<Vec<Vec<pyo3::Bound<'a, PyDict>>>> {
        if analyte_lists.len() != reference_lists.len() {
            // Verify list sizes
            Err(PyValueError::new_err(format!(
                "Number of analyte lists ({0}) and reference lists ({1}) don't match!",
                analyte_lists.len(),
                reference_lists.len()
            )))
        } else {
            let mut jobs: Vec<ORAJob> = Vec::new();
            for (i, analyte_list_vec) in analyte_lists.iter().enumerate() {
                let gmt: Vec<Item> =
                    webgestalt_lib::readers::read_gmt_file(gmt_path.clone()).unwrap();
                let analyte_list: AHashSet<String> = analyte_list_vec.iter().cloned().collect();
                let reference: AHashSet<String> = reference_lists[i].iter().cloned().collect();
                let new_job: ORAJob = ORAJob {
                    gmt: gmt.clone(),
                    interest_list: analyte_list.clone(),
                    reference_list: reference.clone(),
                    config: ORAConfig::default(),
                };
                jobs.push(new_job);
            }
            let rust_result = multilist_ora(
                jobs,
                webgestalt_lib::methods::multilist::MultiListMethod::Meta(
                    webgestalt_lib::methods::multilist::MetaAnalysisMethod::Stouffer,
                ),
                webgestalt_lib::stat::AdjustmentMethod::BH,
            );
            let mut final_results: Vec<Vec<pyo3::Bound<PyDict>>> = Vec::new();
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

    /// Run a meta-analysis ORA with files at the provided paths.
    ///
    /// # Parameters
    /// - `gmt_path` - `String` of the path to the gmt file of interest
    /// - `analyte_list_paths` -  Lists of `String`s of the path to the analyte files of interest.
    /// - `reference_list_paths` - Lists of `String`s of the paths to reference lists.
    ///
    /// # Returns
    ///
    /// Returns a list of a list of dictionaries with the results containing the ORA results for every set.
    ///
    /// The first list contains the results of the meta-analysis. The following lists are the results for each list individually.
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
    /// res = webgestaltpy.meta_ora_from_files("kegg.gmt", ["gene_list1.txt", "gene_list2.txt"], ["reference.txt", "reference.txt"])
    /// ```
    ///
    /// `res` would be a list containing the results of the meta-analysis and each list run
    /// individually. In this example, `res[0]` would look be the results of the meta-analysis.
    /// `res[1]` would be the results from `gene_list1.txt`, `res[2]` would be the results from `gene_list2.txt`, and so on.
    ///
    /// See the documentation for [`webgestaltpy.ora`](./ora.md) for specifics about the format of the results.
    #[pyfunction]
    fn meta_ora_from_files<'a>(
        py: Python<'a>,
        gmt_path: String,
        analyte_list_paths: Vec<String>,
        reference_list_paths: Vec<String>,
    ) -> PyResult<Vec<Vec<pyo3::Bound<'a, PyDict>>>> {
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
            let rust_result = multilist_ora(
                jobs,
                webgestalt_lib::methods::multilist::MultiListMethod::Meta(
                    webgestalt_lib::methods::multilist::MetaAnalysisMethod::Stouffer,
                ),
                webgestalt_lib::stat::AdjustmentMethod::BH,
            );
            let mut final_results: Vec<Vec<pyo3::Bound<PyDict>>> = Vec::new();
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
}
