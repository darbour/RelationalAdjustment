using GLM 

function create_aggregate(adj_mat, values, aggregate_function)
    vals = zeros(size(adj_mat, 1))
    for i in 1:size(adj_mat, 1)
        vals[i] = aggregate_function(values[map(Bool, adj_mat[i,:])])
    end
    return vals
end

function create_covariates(adj_mat, data, aggregate_functions)
    # create the data frames
    columns = names(data)
    for column in columns
        for agg_function in aggregate_functions
            data["$(column)_$(agg_function)"] = create_aggregate(adj_mat, data[column], agg_function)
        end   
    end
    return covariates
end

function update_covariates!(adj_mat, data, update_cov)
    agg_functions = [mean]
    for agg_function in aggregate_functions
        data["$(update_cov)_$(agg_function)"] = create_aggregate(adj_mat, data[update_cov], agg_function)
    end
end

function RPSM(adj_mat, data, treatment)
    agg_functions = [mean]
    covariates = create_covariates(adj_mat, data, agg_functions)
    formula = "$treatment ~ "
    idx = 0
    for column in columns
        if idx > 0
            formula = string(formula, '+')
        end
        idx = idx + 1
        if column != treatment
            formula = string(formula, column)
        end
        for agg_function in agg_functions
            formula = string(formula, '+', "$(column)_$(agg_function)")
        end
    end
    propensity_score_model = GLM(parse(formula), data=covariates, Binomial(), LogitLink())
    
    return propensity_score_model
end

function GSN(adj_mat, data, treatment, burnin, samples, thinning)
    psm = RPSM(adj_mat, data, treatment)
    N = size(adj_mat, 1)
    cur_sample = rand(Binomial(), N)
    loop_data = copy(data)
    update_covariates!(adj_mat, loop_data, treatment)
    samples = zeros(samples, N)
    for i in 1:(burnin+(samples*thinning))
        # choose the next sample at random
        next_var = sample(1:N, 1)
        prediction = rand(Binomial(1, predict(loop_data[next_var,:])[1]), 1)
        loop_data[treatment][next_var] = prediction 
        update_covariates!(adj_mat, loop_data, treatment)
        if i > burnin && i % thinning == 0
            samples[i, :] = loop_data[treatment]
        end
    end
    return samples
end
