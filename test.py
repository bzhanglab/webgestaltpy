import webgestaltpy

# x = webgestaltpy.gsea("data/ktest.gmt", "data/test.rnk")
# print(x[0:2])
# print(
#     "GSEA: There are {0} pathways with significant fdr values".format(
#         len([i for i in x if i["fdr"] < 0.05])
#     )
# )

# y = webgestaltpy.ora("data/kegg.gmt", "data/genelist.txt", "data/reference.txt")
# print(y[0:2])
# print(
#     "ORA: There are {0} pathways with significant fdr values".format(
#         len([i for i in y if i["fdr"] < 0.05])
#     )
# )

# y = webgestaltpy.meta_ora(
#     "data/kegg.gmt",
#     ["data/genelist.txt", "data/genelist.txt"],
#     ["data/reference.txt", "data/reference.txt"],
# )
# print(y[0][0:2])

nta_method = webgestaltpy.NTAMethod.Prioritization

y = webgestaltpy.nta("data/hsapiens_network_CPTAC_Proteomics_OV_entrezgene.net", "data/net_genes.txt", nta_method, 5)
print(y)