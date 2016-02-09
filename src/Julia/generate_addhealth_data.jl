using Distributions;
using Graphs;
using DataFrames;

function get_addhealth_exposure(treatment, node, neighbor_indices)
    bool_t = map(Bool, treatment)
    myt = bool_t[node]
    any_treated_friend = sum(bool_t[neighbor_indices]) .> 0

    if any_treated_friend && myt
        return 4
    elseif myt
        return 3
    elseif any_treated_friend
        return 2
    else
        return 1
    end
end

function generate_addhealth_data()
    nsubjects=626
    g = watts_strogatz_graph(simple_graph(nsubjects, is_directed=false), nsubjects, 4, 0.4)
    adj_mat = adjacency_matrix(g)
    bd = Binomial(1, 0.1)
    z = rand(bd, nsubjects)

    any_treated_friend = mapslices(s -> dot(z, s), adj_mat, 1) .> 0
    outcome = zeros(nsubjects)
    odist = Poisson(2.14)
    for i in 1:nsubjects
        base_outcome = rand(odist, 1)[1]
        myt = map(Bool, z[i])
        if any_treated_friend[i] && myt
            outcome[i] = base_outcome * 2
        elseif myt
            outcome[i] = base_outcome * 1.5
        elseif any_treated_friend[i]
            outcome[i] = base_outcome * 1.25
        else
            outcome[i] = base_outcome
        end
    end

    return(adj_mat, DataFrame(t=z, o=outcome))
end 
