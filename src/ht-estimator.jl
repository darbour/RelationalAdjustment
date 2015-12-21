function ht(exposure, outcome, eprob_fun) 
    subject_id = 1:length(exposure)
    # get all of the unique exposures
    exposures = unique(exposure, 1) 
    for cur_exposure_idxs in size(exposures, 1):
        exposure_comp(idx) = all(exposure[idx, :] .== exposures[cur_exposure_idxs, :])
        subjects = subject_id(map(exposure_comp, subject_id))
    end
end
