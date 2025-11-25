import webgestaltpy

x = webgestaltpy.gsea_from_files("data/kegg.gmt", "data/test.rnk")
print(x[0:2])
print(
    "GSEA: There are {0} pathways with significant fdr values".format(
        len([i for i in x if i["fdr"] < 0.05])
    )
)


def create_rank_list(file_path: str) -> list[tuple[str, float]]:
    """Function that converts rnk file to format for gsea().

    This is for demo purposes. gsea_from_files would be more convenient in this case.
    However, for when the scores are calculated in a Python script, gsea() allows you to directly
    use the results without having to save to a file first.

    """

    res = []
    with open(file_path, "r") as r:
        lines = r.readlines()
    for line in lines:
        if "\t" in line:
            vals = line.split("\t")
            res.append((vals[0], float(vals[1])))
    return res


res = webgestaltpy.gsea("data/kegg.gmt", create_rank_list("data/test.rnk"))

print(res[0:2])  # print first two results

y = webgestaltpy.ora_from_files(
    "data/kegg.gmt", "data/genelist.txt", "data/reference.txt"
)
print(y[0:2])
print(
    "ORA: There are {0} pathways with significant fdr values".format(
        len([i for i in y if i["fdr"] < 0.05])
    )
)

# y = webgestaltpy.meta_ora(
#     "data/kegg.gmt",
#     ["data/genelist.txt", "data/genelist.txt"],
#     ["data/reference.txt", "data/reference.txt"],
# )
# print(y[0][0:2])

# nta_method = webgestaltpy.NTAMethod.Prioritization

# y = webgestaltpy.nta("data/hsapiens_network_CPTAC_Proteomics_OV_entrezgene.net", "data/net_genes.txt", nta_method, 5)
# print(y)
