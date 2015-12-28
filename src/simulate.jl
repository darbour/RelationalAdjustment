#n_subjects = 1000
#burn_in = 10000
#samples = 5000
#thinning = 1000

n_subjects = 1000
burn_in = 10
samples = 500
thinning = 10

include("generate_data.jl")
include("network_sampling.jl")
include("ht_estimator.jl")
include("compute_exposure.jl")

# generate some data
println("Generating data for $n_subjects subjects ...")
adj_mat, dat = generate_data(n_subjects)

# now simulate from the joint exposure distribution
println("Gathering $samples samples ...")
samples = GSN(adj_mat, dat, "t", burn_in, samples, thinning)


# calculate the exposure probs
(obs_exposures, exposure_probs) = get_exposures(dat[:, :t], samples, adj_mat)

println("Estimating means and variances ...")
means, variances = ht(obs_exposures, dat[:o], exposure_probs)

println("Means:")
println(means)

println("Variances:")
println(variances)
