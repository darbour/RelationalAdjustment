
n_subjects = 1000
n_samples = 5000
thinning = n_subjects
burn_in = n_subjects

include("network_sampling.jl")
include("ht_estimator.jl")
include("compute_exposure.jl")

#println("Running Add Health simulation ...")
#include("generate_addhealth_data.jl")
#adj_mat, dat = generate_addhealth_data()
#samples = GSN(adj_mat, dat, "t", burn_in, n_samples, thinning)
#(obs_exposures, exposure_probs) = get_exposures(dat[:, :t], samples, adj_mat, get_addhealth_exposure, 4)
#means, variances = ht(obs_exposures, dat[:o], exposure_probs)

println("Running logistic simulation ...")
include("generate_data.jl")

# generate some data
println("Generating data for $n_subjects subjects ...")
adj_mat, dat = generate_data(n_subjects)

# now simulate from the joint exposure distribution
println("Gathering $n_samples samples ...")
samples = GSN(adj_mat, dat, "t", burn_in, n_samples, thinning)


# calculate the exposure probs
(obs_exposures, exposure_probs) = get_exposures(dat[:, :t], samples, adj_mat, exist_exposure)

# output data, samples, and exposures probabilities
mkpath("results")
for exposure in keys(exposure_probs)
    probs = exposure_probs[exposure]

    # convert to a tabular format for output
    (rowIdx, colIdx) = findn(probs)
    output = DataFrame([Int, Int, Float64], [:row, :col, :prob], size(rowIdx, 1))
    for i in 1:size(rowIdx, 1)
        output[i, :row] = rowIdx[i]
        output[i, :col] = colIdx[i]
        output[i, :prob] = probs[rowIdx[i], colIdx[i]]
    end
    writetable("results/exposure_probs_$(exposure).csv", output)
end

writetable("results/samples.csv", DataFrame(samples))

dat[:exposure] = obs_exposures
writetable("results/data.csv", dat)
writetable("results/network.csv", DataFrame(adj_mat))


println("Estimating means and variances ...")
means, variances = ht(obs_exposures, dat[:o], exposure_probs)

println("Means:")
println(means)

println("Variances:")
println(variances)
