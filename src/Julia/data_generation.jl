using RCall; 

function simple_conditional(idx, adj_mat, observation)
    
end

"create a sample network given a conditional model"
function gibbs_sampler(adj_mat, conditional_function, first_obs, burnin)
    num_nodes = size(adj_mat, 1)
    cur_obs = copy(first_obs)
    for iter in 1:(burnin+1)
        for variable in 1:num_nodes
            cur_obs[variable, :] = conditional_function(variable, adj_mat, cur_obs)
        end
    end
    return cur_obs
end
