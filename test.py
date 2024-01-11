import webgestaltpy

x = webgestaltpy.gsea_from_files("data/ktest.gmt", "data/test.rnk")
print(x[0:2])
print(
    "GSEA: There are {0} pathways with significant fdr values".format(
        len([i for i in x if i["fdr"] < 0.05])
    )
)

y = webgestaltpy.ora_from_files(
    "data/ktest.gmt", "data/genelist.txt", "data/reference.txt"
)
print(y[0:2])
print(
    "ORA: There are {0} pathways with significant fdr values".format(
        len([i for i in y if i["fdr"] < 0.05])
    )
)
