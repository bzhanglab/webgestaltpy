# Usage

For the following examples, the example data can be found here: [data.tar.gz :octicons-download-16:](../data.tar.gz). These examples are derived from <a href="https://www.webgestalt.org/" target="_blank">webgestalt.org :octicons-link-external-16:</a>.


## ORA Example

```python title="ora_test.py"
import WebGestaltPy

res = WebGestaltPy.ora_from_files("kegg.gmt", "genelist.txt", "reference.txt")
```

## Meta-analysis ORA Example


```python title="meta_ora_test.py"
import WebGestaltPy

res = WebGestaltPy.meta_ora_from_files("kegg.gmt", ["genelist.txt", "second_genelist.txt"],
                      ["reference.txt", "reference.txt"])
```

## GSEA Example

```python title="gsea_test.py"
import WebGestaltPy

res = WebGestaltPy.gsea_from_files("kegg.gmt", "test.rnk")
```

## Meta-analysis GSEA Example

```python title="meta_gsea_test.py"
import WebGestaltPy

res = WebGestaltPy.meta_gsea_from_files("kegg.gmt", ["test.rnk", "second_test.rnk"])
```

## NTA Example

```python title="nta_test.py"
import WebGestaltPy
nta_method = webgestaltpy.NTAMethod.Prioritization
res = webgestaltpy.nta_from_files("data/hsapiens_network_CPTAC_Proteomics_OV_entrezgene.net", "data/net_genes.txt", nta_method, 5)
```
