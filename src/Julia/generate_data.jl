using Graphs;
using Distributions;
using StatsBase;
using GLM;
using DataFrames;

function generate_data(nsubjects)
    g = watts_strogatz_graph(simple_graph(nsubjects, is_directed=false), nsubjects, 4, 0.2)
    adj_mat = adjacency_matrix(g)
    degrees = mapslices(sum, adj_mat, 2)[:, 1]

    stdNormal = Normal()
    c1 = rand(stdNormal, nsubjects)
    c2 = adj_mat * c1 ./ degrees + rand(stdNormal, nsubjects)
    c3 = adj_mat * c2 ./ degrees + rand(stdNormal, nsubjects)

    t_prob = 1 ./ (1 + exp(-zscore(c1 + c2 + c3 + rand(stdNormal, nsubjects))))
    treatment = zeros(nsubjects)
    for i in 1:size(treatment, 1)
        treatment[i] = rand(Binomial(1, t_prob[i]), 1)[1]
    end

    burn_in = 100 
    iter = 0 
    while(iter <= burn_in) 
        iter = iter + 1 
        cur_t_friends = (adj_mat * treatment) ./ degrees

        t_prob = 1 ./ (1 + exp(-zscore(c1 + c2 + c3 + rand(stdNormal, nsubjects)) - zscore(cur_t_friends)))
        for i in 1:size(treatment, 1)
            treatment[i] = rand(Binomial(1, t_prob[i]), 1)[1]
        end 
    end 
    t_friends = (adj_mat * treatment) ./ degrees
    
    covariates = DataFrame(t = treatment, c1 = c1, c2 = c2, c3 = c3, t_friends=t_friends)
    println(covariates[1:10, :])
    model = glm(t ~ c1 + c2 + c3 + t_friends, covariates, Binomial(), LogitLink())
    println(model)

    o = covariates[:c1] + covariates[:c2] + covariates[:c3] + covariates[:t] + covariates[:t_friends] + rand(stdNormal, nsubjects)
    covariates[:o] = o
    println(lm(o ~ c1 + c2 + c3 + t_friends + t, covariates))
    dat = DataFrame(t=treatment, c1=c1, c2=c2, c3=c3, o=o)

    return (adj_mat, dat)
end


