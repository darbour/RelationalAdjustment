using GLM 
using DataFrames

function create_aggregate(adj_mat, values, aggregate_function)
    vals = zeros(size(adj_mat, 1))
    for i in 1:size(adj_mat, 1)
        vals[i] = aggregate_function(values[map(Bool, adj_mat[i,:]), :])
    end
    return vals
end

function create_covariates(adj_mat, data, aggregate_functions)
    # create the data frames
    columns = names(data)
    for column in columns
        if column == :o
            continue
        end
        for agg_function in aggregate_functions
            agg = create_aggregate(adj_mat, data[column], agg_function)
            orig_names = names(data)
            data = hcat(data, agg)
            new_names = vcat(orig_names, parse("$(column)_$(agg_function)"))
            names!(data, new_names)
        end   
    end
    return data
end

function update_covariates!(adj_mat, data, update_cov )
    agg_functions = [mean]
    for agg_function in agg_functions
        data[parse("$(update_cov)_$(agg_function)")] = create_aggregate(adj_mat, data[parse(update_cov)], agg_function)
    end
end

function RPSM(adj_mat, data, treatment)
    agg_functions = [mean]
    covariates = create_covariates(adj_mat, data, agg_functions)
    a = []
    for column in names(covariates)
        if column != :t && column != :o
            append!(a, [column])
        end
    end
    fm = Formula(parse("$treatment"), Expr(:call, :+, a))
    propensity_score_model = glm(fm, covariates, Binomial(), LogitLink())
    println(propensity_score_model)
    
    return (propensity_score_model, covariates)
end

function GSN(adj_mat, data, treatment, burnin, nsamples, thinning)
    (psm, covariates) = RPSM(adj_mat, data, treatment)
    N = size(adj_mat, 1)
    cur_sample = rand(Binomial(), N)
    loop_data = copy(covariates)
    update_covariates!(adj_mat, loop_data, treatment)
    samples = zeros(nsamples, N)
    for i in 1:(burnin+(nsamples*thinning))
        if i % 100 == 0
            println(i)
        end
        # choose the next sample at random
        next_var = sample(1:N, 1)
        prediction = rand(Binomial(1, predict(psm, loop_data[next_var,:])[1]), 1)
        loop_data[:t][next_var] = prediction 
        update_covariates!(adj_mat, loop_data, treatment)
        if i > burnin && i % thinning == 0
            samples[i-burnin, :] = loop_data[:t]
        end
    end
    println(samples[1:10, 1:10])
    return samples
end

adj_mat = readcsv("network.csv")
data = readtable("network_attrs.csv")
GSN(adj_mat, data, "t", 50, 300, 1)
