using GLM 
using DataFrames

function create_aggregate(adj_mat, values, aggregate_function, update_rows)
    vals = zeros(length(update_rows))
    for (v,i) in enumerate(update_rows)
        (_, indices, _) = findnz(adj_mat[i,:])
        vals[v] = aggregate_function(values[indices])
    end
    return vals
end

function create_aggregate(adj_mat, values, aggregate_function)
    update_rows = 1:size(adj_mat, 1)
    return create_aggregate(adj_mat, values, aggregate_function, update_rows)
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

function update_covariates!(adj_mat, data, update_cov, update_indices)
    agg_functions = [mean]
    for agg_function in agg_functions
        data[parse("$(update_cov)_$(agg_function)")][update_indices] = create_aggregate(adj_mat, data[parse(update_cov)], agg_function, update_indices)
    end
end

function update_covariates!(adj_mat, data, update_cov)
    update_indices = 1:size(adj_mat, 1)
    return update_covariates!(adj_mat, data, update_cov, update_indices)
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
    cur_sample_idx = 1
    all_idxs = 1:N
    for i in 1:(burnin+(nsamples*thinning))
        if i % 1000 == 0
            println(i)
        end
        # choose the next sample at random
        next_var = sample(1:N, 1)
        # draw a value for the treatment
        prediction = rand(Binomial(1, predict(psm, loop_data[next_var,:])[1]), 1)
        # update the treatment vector
        loop_data[next_var, :t] = prediction 
        # grab the indices of neighbors
        (adj_vals, _, _) = findnz(adj_mat[:, next_var])
        # include the current node in the set of aggregates to be updated
        update_vals = [next_var;adj_vals]
        # udpate covariates
        update_covariates!(adj_mat, loop_data, treatment, update_vals)
        # save our sample? 
        if i > burnin && i % thinning == 0
            # Yes. Do it.
            samples[cur_sample_idx, :] = copy(loop_data[:t])
            cur_sample_idx += 1
        end
    end
    println(samples[1:10, 1:10])
    return samples
end

adj_mat = readcsv("network.csv")
data = readtable("network_attrs.csv")
GSN(adj_mat, data, "t", 50, 300, 1000)
