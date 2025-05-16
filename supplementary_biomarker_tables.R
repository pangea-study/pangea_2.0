
# distribution of biomarkers in the cohorts at baseline (first pangea_bm obs)
biomarker_dists <- data_validation[pangea_bm_obs_id == 1, 
         rbind(hgb = quantile(hgb, c(0.5, 0.25, 0.75), na.rm = T),
               mspike = quantile(mspike, c(0.5, 0.25, 0.75), na.rm = T),
               iuratio = quantile(iuratio, c(0.5, 0.25, 0.75), na.rm = T),
               creatinine = quantile(creatinine, c(0.5, 0.25, 0.75), na.rm = T),
               bmpc = quantile(plasmacells, c(0.5, 0.25, 0.75), na.rm = T))]
			   

# frequency of BM measurements
bm_freq <- data_validation[, any(!is.na(plasmacells)), 
                           by = participant_id][, c(sum(V1), mean(V1))]
names(bm_freq) <- c("N not NA", "% not NA")


# frequency of the trajectories
traj_pt_freq <- sapply(data_validation[pangea_no_bm_ok == T, 
         list(hgb_traj = as.numeric(any(hgb_traj == 1)),
              mspike_traj = as.numeric(any(mspike_traj == 1)),
              iuratio_traj = as.numeric(any(iuratio_traj == 1)),
              creatinine_traj = as.numeric(any(creatinine_traj == 1))), 
         by = participant_id][,2:5], 
       mean, na.rm = T)
traj_obs_freq <- sapply(data_validation[pangea_no_bm_ok == T, 
                 list(hgb_traj, mspike_traj, iuratio_traj, creatinine_traj)], 
       mean, na.rm = T)
	   
	   
## by progression status
data_validation[, ever_mm := any(prog_mm == 1), by = participant_id]
traj_pt_prog_freq <- sapply(data_validation[pangea_no_bm_ok == T & ever_mm == 1, 
	                   list(hgb_traj = as.numeric(any(hgb_traj == 1)),
	                        mspike_traj = as.numeric(any(mspike_traj == 1)),
	                        iuratio_traj = as.numeric(any(iuratio_traj == 1)),
	                        creatinine_traj = as.numeric(any(creatinine_traj == 1))), 
	                   by = participant_id][,2:5], 
	          mean, na.rm = T)
traj_pt_nonprog_freq <- sapply(data_validation[pangea_no_bm_ok == T & ever_mm == 0, 
	                   list(hgb_traj = as.numeric(any(hgb_traj == 1)),
	                        mspike_traj = as.numeric(any(mspike_traj == 1)),
	                        iuratio_traj = as.numeric(any(iuratio_traj == 1)),
	                        creatinine_traj = as.numeric(any(creatinine_traj == 1))), 
	                   by = participant_id][,2:5], 
					   mean, na.rm = T)



# cohort-specific BM model
data_validation_pangea_bm_traj_mod <- coxph(model_traj_bm_alt$formula, data = data_validation)
summary(data_validation_pangea_bm_traj_mod)

# cohort-specific no-BM model
data_validation_pangea_no_bm_traj_mod <- coxph(model_traj_no_bm_alt$formula, data = data_validation)
summary(data_validation_pangea_no_bm_traj_mod)


## collect all results (other than models) into one list
supp_biomarker_list <- list(biomarker_dists = biomarker_dists,
                            bm_freq = bm_freq,
                            traj_pt_freq = traj_pt_freq,
                            traj_obs_freq = traj_obs_freq,
                            traj_pt_prog_freq = traj_pt_prog_freq,
                            traj_pt_nonprog_freq = traj_pt_nonprog_freq)

