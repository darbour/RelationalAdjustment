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
