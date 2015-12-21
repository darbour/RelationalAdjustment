#exposures = ["0,[0,0.1]" , "0,(0.1,0.2]" , "0,(0.2,0.3]" , 
#             "0,(0.3,0.4]" , "0,(0.4,0.5]" , "0,(0.5,0.6]" , 
#             "0,(0.6,0.7]" , "0,(0.7,0.8]" , "0,(0.8,0.9]" , 
#             "0,(0.9,1]" , 
#             "1,[0,0.1]" , "1,(0.1,0.2]" , "1,(0.2,0.3]" , 
#             "1,(0.3,0.4]" , "1,(0.4,0.5]" , "1,(0.5,0.6]" , 
#             "1,(0.6,0.7]" , "1,(0.7,0.8]" , "1,(0.8,0.9]" , 
#             "1,(0.9,1]"]
using DataFrames;

function prop_exposure(treatment, node, neighbor_indices)
    tfrac = sum(treatment[neighbor_indices]) / size(neighbor_indices, 1)
    exposure_idx = (treatment[node] * 10) + floor(min(tfrac, 0.99) / 0.1)
    return(exposure_idx+1)
end

function get_exposures(obs_treatment, treatment_samples, adj_mat, exposurefun=prop_exposure, num_exposures=20)

    nsubjects = size(treatment_samples, 2)
    nsamples = size(treatment_samples, 1)

    exposure_mat = zeros(UInt8, nsamples, nsubjects)
    obs_exposures = zeros(Int, nsubjects)
    for s in 1:nsubjects
        (_, neighbors, _) = findnz(adj_mat[s, :])
        obs_exposures[s] = exposurefun(obs_treatment, s, neighbors)
    end

    # map nodes to exposures
    println("Performing exposure mapping ...")
    for col_idx in 1:nsubjects
        (_, neighbors, _) = findnz(adj_mat[col_idx, :])
        for row_idx in 1:nsamples
            exposure_idx = exposurefun(treatment_samples[row_idx, :], col_idx, neighbors)
            exposure_mat[row_idx, col_idx] = exposure_idx
        end
    end

    # check whether nodes share an exposure
    exposure_probs = Dict()
    for exposure in 1:num_exposures
        exposure_probs[exposure] = spzeros(nsubjects, nsubjects)
    end

    increment = 1.0 / nsamples
    split_point = ceil(Int, nsubjects / 2.0)
    for sample in 1:nsamples
        if sample % 100 == 0
            println("Computing joint exposure probs, sample $sample/$nsamples ...")
        end

        # populate symmetric matrix
        for i in 1:split_point
            for j in i:nsubjects 
                jval = exposure_mat[sample, j] 
                if jval == exposure_mat[sample, i]
                    exposure_probs[jval][i, j] += increment
                    if i != j
                        exposure_probs[jval][j, i] += increment
                    end
                end
            end
        end
    end

    return(obs_exposures, exposure_probs)
end

samples = readcsv("samples1000.csv")
adj_mat = readcsv("network1000.csv")
data = readtable("network_attrs1000.csv")
(obs_exposures, exposure_probs) = get_exposures(data[:, :t], samples, adj_mat)

for exposure in 1:20
    probs = exposure_probs[exposure]
    # convert to a tabular format for writing
    (rowIdx, colIdx) = findn(probs)
    output = DataFrame([Int, Int, Float64], [:row, :col, :prob], size(rowIdx, 1))
    for i in 1:size(rowIdx, 1)
        output[i, :row] = rowIdx[i]
        output[i, :col] = colIdx[i]
        output[i, :prob] = probs[rowIdx[i], colIdx[i]]
    end
    writetable("network1000_exposure$(exposure).csv", output)
end
data[:exposure] = obs_exposures

writetable("network_attrs1000e.csv", data)

