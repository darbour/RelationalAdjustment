include("network_sampling.jl")

adj_mat = readcsv("network.csv")
data = readtable("network_attrs.csv")
samples = GSN(adj_mat, data, "t", 10000, 50000, 2000)
