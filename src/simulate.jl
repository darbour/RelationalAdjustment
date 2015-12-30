
n_subjects = 1000
burn_in = 10000
n_samples = 50000
thinning = 1000

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
(obs_exposures, exposure_probs) = get_exposures(dat[:, :t], samples, adj_mat)

println("Estimating means and variances ...")
means, variances = ht(obs_exposures, dat[:o], exposure_probs)

println("Means:")
println(means)

println("Variances:")
println(variances)
