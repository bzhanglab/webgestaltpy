# Usage

For the following examples, the example data can be found here: [data.tar.gz :octicons-download-16:](../data.tar.gz). These examples are derived from <a href="https://www.webgestalt.org/" target="_blank">webgestalt.org :octicons-link-external-16:</a>.


## ORA Example

```python title="ora_test.py"
import WebGestaltPy

WebGestaltPy.ora("kegg.gmt", "genelist.txt", "reference.txt")
```

## Meta-analysis ORA Example


```python title="meta_ora_test.py"
import WebGestaltPy

WebGestaltPy.meta_ora("kegg.gmt", ["genelist.txt", "second_genelist.txt"],
                      ["reference.txt", "reference.txt"])
```

## GSEA Example

```python title="ora_test.py"
import WebGestaltPy

WebGestaltPy.gsea("kegg.gmt", "test.rnk")
```
