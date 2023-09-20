import webgestaltpy

x = webgestaltpy.gsea_from_files("data/ktest.gmt", "data/test.rnk")
print("There are {0} pathways with significant p values".format(len([i for i in x if i.p < 0.05])))
