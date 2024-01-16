import webgestaltpy

# x = webgestaltpy.gsea("data/ktest.gmt", "data/test.rnk")
# print(x[0:2])
# print(
#     "GSEA: There are {0} pathways with significant fdr values".format(
#         len([i for i in x if i["fdr"] < 0.05])
#     )
# )

y = webgestaltpy.ora("data/ktest.gmt", "data/genelist.txt", "data/reference.txt")
print(y[0:2])
print(
    "ORA: There are {0} pathways with significant fdr values".format(
        len([i for i in y if i["fdr"] < 0.05])
    )
)

y = webgestaltpy.meta_ora(
    "data/ktest.gmt",
    ["data/genelist.txt", "data/genelist.txt"],
    ["data/reference.txt", "data/reference.txt"],
)
print(y[0][0:2])
