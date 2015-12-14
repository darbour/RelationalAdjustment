using RCall

rcopy("source('generate_data.R')")

## Parameters
n_subjects = 100
burn_in = 10000
samples = 5000
thinning = 1000

# generate some data
data, adj_mat = rcopy("generate.data($(n_subjects))")

# now simulate from the joint exposure distribution
include("network_sampling.jl")
samples = GSN(adj_mat, data, "t", burn_in, samples, thinning)

# calculate the exposure probs
(obs_exposures, exposure_probs) = get_exposures(data[:, :t], samples, adj_mat)
for exposure in 1:20
    probs = exposure_probs[exposure]
    # convert to a tabular format for writing
    (rowIdx, colIdx) = findn(probs)
    output = DataFrame([Int, Int, Float64], [:row, :col, :prob], size(rowIdx, 1))
    for i in 1:size(rowIdx, 1)
        output[i, :row] = rowIdx[i]
        output[i, :col] = colIdx[i]
        output[i, :prob] = probs[rowIdx[i], colIdx[i]]
    end
    writetable("network1000_exposure$(exposure).csv", output)
end
data[:exposure] = obs_exposures

writetable("network_attrs1000e.csv", data)
