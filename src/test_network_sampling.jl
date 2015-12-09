include("network_sampling.jl")

adj_mat = readcsv("network.csv")
data = readtable("network_attrs.csv")
GSN(adj_mat, data, "t", 10000, 5000, 1000)
