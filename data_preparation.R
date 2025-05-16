## add NA ethnicity variable if missing
if (!"ethnicity" %in% names(data_validation)) {
  data_validation[, ethnicity := NA]
}

## Logarithm of iuratio and creatinine
data_validation[,
     c("log_iuratio", "log_creatinine") := list(log(iuratio),
                                                log(creatinine))]

## indicators for whether the variables needed for the pangea models are all observed
data_validation[, pangea_bm_ok := !(is.na(plasmacells) | 
                                      is.na(mspike) | 
                                      is.na(iuratio) | 
                                      is.na(creatinine) | 
                                      is.na(age))]

data_validation[, pangea_no_bm_ok := !(is.na(mspike) | 
                                       is.na(iuratio) | 
                                       is.na(creatinine) | 
                                       is.na(age))]

## patient-specific ids for serial pangea bm and pangea no bm observations
data_validation[pangea_bm_ok == T, 
                pangea_bm_obs_id := cumsum(pangea_bm_ok), 
                by = participant_id]

data_validation[pangea_no_bm_ok == T, 
                pangea_no_bm_obs_id := cumsum(pangea_no_bm_ok), 
                by = participant_id]

## Some needed time-related variables
data_validation[, time_between_obs := tstop - tstart]
data_validation[, last_date := max(tstop), by = participant_id]
data_validation[, ever_mm := any(prog_mm == 1), by = participant_id]

## Compute trajectories and add to dataset
data_validation <- computation_trajectories(data_validation)

## Compute risk scores and add to dataset
data_validation <- computation_risks(data_validation,
                                     models_folder = "models/",
                                     models = c("pangea_bm_traj", "pangea_no_bm_traj",
                                                "pangea_bm_traj_alt", "pangea_no_bm_traj_alt",
                                                "pangea_bm", "pangea_no_bm", "rolling_20_2_20"))
		# NOTE that the PANGEA 2.0 models are the pangea_bm_traj_alt and pangea_no_bm_traj_alt objects
		# their risk predictions (after calling the computation_risks functions) have analagous variable names

## Compute version of data where trajectories are set to 0 (used to assess trajectory models
## when no patient history is available)
data_validation_traj0 <- data.table(data_validation)
data_validation_traj0[, mspike_traj := 0]
data_validation_traj0[, creatinine_traj := 0]
data_validation_traj0[, iuratio_traj := 0]
data_validation_traj0[, hgb_traj := 0]
data_validation_traj0[, traj_hgb := 0]

data_validation_traj0 <- computation_risks(data_validation_traj0,
                                     models_folder = "models/",
                                     models = c("pangea_bm_traj", "pangea_no_bm_traj",
                                                "pangea_bm_traj_alt", "pangea_no_bm_traj_alt",
                                                "pangea_bm", "pangea_no_bm", "rolling_20_2_20"))

