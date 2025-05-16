#### Compute the risk categorization

# function to compute reclassification metrics for pangea vs. 20/2/20
# given specific thresholds to stratify the continuous pangea risk scores
pangea_reclassification <- function(cutoffs_pangea_bm_traj, cutoffs_pangea_no_bm_traj,
                                    cutoffs_pangea_bm, cutoffs_pangea_no_bm,
                                    pangea_labels,
                                    include_pangea_alt = T) {
  
  data_validation <- data.table(data_validation)
  
  ## compute risk categories for each of the models 
  
  # pangea trajectory models
  data_validation[, pangea_bm_traj_risk_cat := cut(pangea_bm_traj_2_year_risk,
                                                   breaks = cutoffs_pangea_bm_traj,
                                                   labels = pangea_labels)]
  data_validation[, pangea_bm_traj_risk_cat_num := as.numeric(pangea_bm_traj_risk_cat)]
  
  data_validation[, pangea_no_bm_traj_risk_cat := cut(pangea_no_bm_traj_2_year_risk,
                                                      breaks = cutoffs_pangea_no_bm_traj,
                                                      labels = pangea_labels)]
  data_validation[, pangea_no_bm_traj_risk_cat_num := as.numeric(pangea_no_bm_traj_risk_cat)]
  
  # pangea trajectory alt models
  if (include_pangea_alt == T) {
    data_validation[, pangea_bm_traj_alt_risk_cat := cut(pangea_bm_traj_alt_2_year_risk,
                                                         breaks = cutoffs_pangea_bm_traj,
                                                         labels = pangea_labels)]
    data_validation[, pangea_bm_traj_alt_risk_cat_num := as.numeric(pangea_bm_traj_alt_risk_cat)]
    
    data_validation[, pangea_no_bm_traj_alt_risk_cat := cut(pangea_no_bm_traj_alt_2_year_risk,
                                                            breaks = cutoffs_pangea_no_bm_traj,
                                                            labels = pangea_labels)]
    data_validation[, pangea_no_bm_traj_alt_risk_cat_num := as.numeric(pangea_no_bm_traj_alt_risk_cat)]
  }
  
  # pangea 1.0 models
  data_validation[, pangea_bm_risk_cat := cut(pangea_bm_2_year_risk,
                                              breaks = cutoffs_pangea_bm,
                                              labels = pangea_labels)]
  data_validation[, pangea_bm_risk_cat_num := as.numeric(pangea_bm_risk_cat)]
  
  data_validation[, pangea_no_bm_risk_cat := cut(pangea_no_bm_2_year_risk,
                                                 breaks = cutoffs_pangea_no_bm,
                                                 labels = pangea_labels)]
  data_validation[, pangea_no_bm_risk_cat_num := as.numeric(pangea_no_bm_risk_cat)]
  
  # mayo (20/2/20)
  labels_20_2_20 <- c("Low", "Intermediate", "High")
  data_validation[, rolling_20_2_20 := factor(rolling_20_2_20, 
                                              levels = labels_20_2_20)]
  data_validation[, table(pangea_bm_traj_risk_cat, rolling_20_2_20)]
  
  data_validation[, rolling_20_2_20_num := as.numeric(rolling_20_2_20)]
  
  
  ## compute reclassification tables
  
  # full sample
  reclass_bm_traj_all <- data_validation[, addmargins(table(rolling_20_2_20, 
                                                            pangea_bm_traj_risk_cat))]
  reclass_no_bm_traj_all <- data_validation[, addmargins(table(rolling_20_2_20, 
                                                               pangea_no_bm_traj_risk_cat))]
  
  if (include_pangea_alt == T) {
    reclass_bm_traj_alt_all <- data_validation[, addmargins(table(rolling_20_2_20, 
                                                              pangea_bm_traj_alt_risk_cat))]
    reclass_no_bm_traj_alt_all <- data_validation[, addmargins(table(rolling_20_2_20, 
                                                                 pangea_no_bm_traj_alt_risk_cat))]
  }
  
  reclass_bm_all <- data_validation[, addmargins(table(rolling_20_2_20, 
                                                       pangea_bm_risk_cat))]
  reclass_no_bm_all <- data_validation[, addmargins(table(rolling_20_2_20, 
                                                          pangea_no_bm_risk_cat))]
  
  # progressors (last obseration)
  reclass_bm_traj_prog_last <- data_validation[prog_mm == 1,
                                          addmargins(table(rolling_20_2_20, pangea_bm_traj_risk_cat))]
  reclass_no_bm_traj_prog_last <- data_validation[prog_mm == 1,
                                             addmargins(table(rolling_20_2_20, pangea_no_bm_traj_risk_cat))]
  
  if (include_pangea_alt == T) {
    reclass_bm_traj_alt_prog_last <- data_validation[prog_mm == 1,
                                            addmargins(table(rolling_20_2_20, pangea_bm_traj_alt_risk_cat))]
    reclass_no_bm_traj_alt_prog_last <- data_validation[prog_mm == 1,
                                               addmargins(table(rolling_20_2_20, pangea_no_bm_traj_alt_risk_cat))]
    
  }
  
  reclass_bm_prog_last <- data_validation[prog_mm == 1,
                                     addmargins(table(rolling_20_2_20, pangea_bm_risk_cat))]
  reclass_no_bm_prog_last <- data_validation[prog_mm == 1,
                                        addmargins(table(rolling_20_2_20, pangea_no_bm_risk_cat))]
  
  # progressors
  reclass_bm_traj_prog <- data_validation[ever_mm == 1,
                                          addmargins(table(rolling_20_2_20, pangea_bm_traj_risk_cat))]
  reclass_no_bm_traj_prog <- data_validation[ever_mm == 1,
                                             addmargins(table(rolling_20_2_20, pangea_no_bm_traj_risk_cat))]
  
  if (include_pangea_alt == T) {
    reclass_bm_traj_alt_prog <- data_validation[ever_mm == 1,
                                                addmargins(table(rolling_20_2_20, pangea_bm_traj_alt_risk_cat))]
    reclass_no_bm_traj_alt_prog <- data_validation[ever_mm == 1,
                                                   addmargins(table(rolling_20_2_20, pangea_no_bm_traj_alt_risk_cat))]
    
  }
  
  reclass_bm_prog <- data_validation[ever_mm == 1,
                                     addmargins(table(rolling_20_2_20, pangea_bm_risk_cat))]
  reclass_no_bm_prog <- data_validation[ever_mm == 1,
                                        addmargins(table(rolling_20_2_20, pangea_no_bm_risk_cat))]
  
  # non-progressors
  reclass_bm_traj_nonprog <- data_validation[ever_mm == 0,
                                          addmargins(table(rolling_20_2_20, pangea_bm_traj_risk_cat))]
  reclass_no_bm_traj_nonprog <- data_validation[ever_mm == 0,
                                             addmargins(table(rolling_20_2_20, pangea_no_bm_traj_risk_cat))]
  
  if (include_pangea_alt == T) {
    reclass_bm_traj_alt_nonprog <- data_validation[ever_mm == 0,
                                                addmargins(table(rolling_20_2_20, pangea_bm_traj_alt_risk_cat))]
    reclass_no_bm_traj_alt_nonprog <- data_validation[ever_mm == 0,
                                                   addmargins(table(rolling_20_2_20, pangea_no_bm_traj_alt_risk_cat))]
    
  }
  
  reclass_bm_nonprog <- data_validation[ever_mm == 1,
                                     addmargins(table(rolling_20_2_20, pangea_bm_risk_cat))]
  reclass_no_bm_nonprog <- data_validation[ever_mm == 1,
                                        addmargins(table(rolling_20_2_20, pangea_no_bm_risk_cat))]
  
  # put results into a list
  reclass_tables <- list(# all patients
                         reclass_bm_traj_all = reclass_bm_traj_all, 
                         reclass_no_bm_traj_all = reclass_no_bm_traj_all,
                         reclass_bm_traj_alt_all = if (include_pangea_alt) reclass_bm_traj_alt_all,
                         reclass_no_bm_traj_alt_all = if (include_pangea_alt) reclass_no_bm_traj_alt_all,
                         reclass_bm_all = reclass_bm_all,
                         reclass_no_bm_all = reclass_no_bm_all,
                         # progressors last obs
                         reclass_bm_traj_prog_last = reclass_bm_traj_prog_last, 
                         reclass_no_bm_traj_prog_last = reclass_no_bm_traj_prog_last,
                         reclass_bm_traj_alt_prog_last = if (include_pangea_alt) reclass_bm_traj_alt_prog_last,
                         reclass_no_bm_traj_alt_prog_last = if (include_pangea_alt) reclass_no_bm_traj_alt_prog_last,
                         reclass_bm_prog_last = reclass_bm_prog_last,
                         reclass_no_bm_prog_last = reclass_no_bm_prog_last,
                         # progressors
                         reclass_bm_traj_prog = reclass_bm_traj_prog, 
                         reclass_no_bm_traj_prog = reclass_no_bm_traj_prog,
                         reclass_bm_traj_alt_prog = if (include_pangea_alt) reclass_bm_traj_alt_prog,
                         reclass_no_bm_traj_alt_prog = if (include_pangea_alt) reclass_no_bm_traj_alt_prog,
                         reclass_bm_prog = reclass_bm_prog,
                         reclass_no_bm_prog = reclass_no_bm_prog,
                         # non-progressors
                         reclass_bm_traj_nonprog = reclass_bm_traj_nonprog, 
                         reclass_no_bm_traj_nonprog = reclass_no_bm_traj_nonprog,
                         reclass_bm_traj_alt_nonprog = if (include_pangea_alt) reclass_bm_traj_alt_nonprog,
                         reclass_no_bm_traj_alt_nonprog = if (include_pangea_alt) reclass_no_bm_traj_alt_nonprog,
                         reclass_bm_nonprog = reclass_bm_nonprog,
                         reclass_no_bm_nonprog = reclass_no_bm_nonprog)
  
  ## fraction of progressors ever classified high
  ever_obs <- data_validation[, list(ever_mayo_high = any(rolling_20_2_20 == "High"),
                                     ever_pangea_bm_traj_high = any(pangea_bm_traj_risk_cat == "High"),
                                     ever_pangea_bm_traj_alt_high = any(pangea_bm_traj_alt_risk_cat == "High"),
                                     ever_pangea_bm_high = any(pangea_bm_risk_cat == "High"),
                                     ever_pangea_no_bm_traj_high = any(pangea_no_bm_traj_risk_cat == "High"),
                                     ever_pangea_no_bm_traj_alt_high = any(pangea_no_bm_traj_alt_risk_cat == "High"),
                                     ever_pangea_no_bm_high = any(pangea_no_bm_risk_cat == "High"),
                                     max_pangea_bm_traj = max(pangea_bm_traj_2_year_risk),
                                     ever_mm = ever_mm[1]),
                              by = participant_id]
  
  # mayo
  ever_mayo_high <- ever_obs[, addmargins(table(ever_mm, ever_mayo_high))]
  ever_mayo_high_ss <- c(sens = ever_mayo_high[2,2]/ever_mayo_high[2,3],
                         spec = ever_mayo_high[1,1]/ever_mayo_high[1,3],
                         ppv = ever_mayo_high[2,2]/ever_mayo_high[3,2],
                         npv = ever_mayo_high[1,1]/ever_mayo_high[3,1])
  
  # PANGEA trajectory models
  ever_bm_traj_high <- ever_obs[, addmargins(table(ever_mm, ever_pangea_bm_traj_high))]
  ever_bm_traj_high_ss <- c(sens = ever_bm_traj_high[2,2]/ever_bm_traj_high[2,3],
                            spec = ever_bm_traj_high[1,1]/ever_bm_traj_high[1,3],
                            ppv = ever_bm_traj_high[2,2]/ever_bm_traj_high[3,2],
                            npv = ever_bm_traj_high[1,1]/ever_bm_traj_high[3,1])
  
  ever_no_bm_traj_high <- ever_obs[, addmargins(table(ever_mm, ever_pangea_no_bm_traj_high))]
  ever_no_bm_traj_high_ss <- c(sens = ever_no_bm_traj_high[2,2]/ever_no_bm_traj_high[2,3],
                            spec = ever_no_bm_traj_high[1,1]/ever_no_bm_traj_high[1,3],
                            ppv = ever_no_bm_traj_high[2,2]/ever_no_bm_traj_high[3,2],
                            npv = ever_no_bm_traj_high[1,1]/ever_no_bm_traj_high[3,1])
  
  # PANGEA trajectory alt models
  if (include_pangea_alt == T) {
    ever_bm_traj_alt_high <- ever_obs[, addmargins(table(ever_mm, ever_pangea_bm_traj_alt_high))]
    ever_bm_traj_alt_high_ss <- c(sens = ever_bm_traj_alt_high[2,2]/ever_bm_traj_alt_high[2,3],
                                  spec = ever_bm_traj_alt_high[1,1]/ever_bm_traj_alt_high[1,3],
                                  ppv = ever_bm_traj_alt_high[2,2]/ever_bm_traj_alt_high[3,2],
                                  npv = ever_bm_traj_alt_high[1,1]/ever_bm_traj_alt_high[3,1])
    
    ever_no_bm_traj_alt_high <- ever_obs[, addmargins(table(ever_mm, ever_pangea_no_bm_traj_alt_high))]
    ever_no_bm_traj_alt_high_ss <- c(sens = ever_no_bm_traj_alt_high[2,2]/ever_no_bm_traj_alt_high[2,3],
                                     spec = ever_no_bm_traj_alt_high[1,1]/ever_no_bm_traj_alt_high[1,3],
                                     ppv = ever_no_bm_traj_alt_high[2,2]/ever_no_bm_traj_alt_high[3,2],
                                     npv = ever_no_bm_traj_alt_high[1,1]/ever_no_bm_traj_alt_high[3,1])
  }
  
  # PANGEA 1.0 models
  ever_bm_high <- ever_obs[, addmargins(table(ever_mm, ever_pangea_bm_high))]
  ever_bm_high_ss <- c(sens = ever_bm_high[2,2]/ever_bm_high[2,3],
                            spec = ever_bm_high[1,1]/ever_bm_high[1,3],
                            ppv = ever_bm_high[2,2]/ever_bm_high[3,2],
                            npv = ever_bm_high[1,1]/ever_bm_high[3,1])
  
  ever_no_bm_high <- ever_obs[, addmargins(table(ever_mm, ever_pangea_no_bm_high))]
  ever_no_bm_high_ss <- c(sens = ever_no_bm_high[2,2]/ever_no_bm_high[2,3],
                               spec = ever_no_bm_high[1,1]/ever_no_bm_high[1,3],
                               ppv = ever_no_bm_high[2,2]/ever_no_bm_high[3,2],
                               npv = ever_no_bm_high[1,1]/ever_no_bm_high[3,1])
  
  # put into list
  ever_high <- list(tables = list(ever_mayo_high = ever_mayo_high,
                                  ever_bm_traj_high = ever_bm_traj_high,
                                  ever_no_bm_traj_high = ever_no_bm_traj_high,
                                  ever_bm_traj_alt_high = if(include_pangea_alt) ever_bm_traj_alt_high,
                                  ever_no_bm_traj_alt_high = if(include_pangea_alt) ever_no_bm_traj_alt_high,
                                  ever_bm_high = ever_bm_high,
                                  ever_no_bm_high = ever_no_bm_high),
                    ss = cbind(ever_mayo_high_ss = ever_mayo_high_ss,
                              ever_bm_traj_high_ss = ever_bm_traj_high_ss,
                              ever_no_bm_traj_high_ss = ever_no_bm_traj_high_ss,
                              ever_bm_traj_alt_high_ss = if(include_pangea_alt) ever_bm_traj_alt_high_ss,
                              ever_no_bm_traj_alt_high_ss = if(include_pangea_alt) ever_no_bm_traj_alt_high_ss,
                              ever_bm_high_ss = ever_bm_high_ss,
                              ever_no_bm_high_ss = ever_no_bm_high_ss))
  
  
  ## time-dependent predictive value
  top_cut <- tail(cutoffs_clinical, 2)[1]
  
  # pangea trajectory models
  sspn_pangea_bm_traj <- time_dep_sspn(dat = data_validation, 
                                       marker_var = "pangea_bm_traj_2_year_risk",
                                       obs_ok_var = "pangea_bm_ok",
                                       high_cut = top_cut)
  sspn_pangea_no_bm_traj <- time_dep_sspn(dat = data_validation, 
                                          marker_var = "pangea_no_bm_traj_2_year_risk",
                                          obs_ok_var = "pangea_no_bm_ok",
                                          high_cut = top_cut)
  
  # pangea trajectory alt models
  sspn_pangea_bm_traj_alt <- time_dep_sspn(dat = data_validation, 
                                  marker_var = "pangea_bm_traj_alt_2_year_risk",
                                  obs_ok_var = "pangea_bm_ok",
                                  high_cut = top_cut)
  sspn_pangea_no_bm_traj_alt <- time_dep_sspn(dat = data_validation, 
                                           marker_var = "pangea_no_bm_traj_alt_2_year_risk",
                                           obs_ok_var = "pangea_no_bm_ok",
                                           high_cut = top_cut)
  
  # pangea 1.0 models
  sspn_pangea_bm <- time_dep_sspn(dat = data_validation, 
                                  marker_var = "pangea_bm_2_year_risk",
                                  obs_ok_var = "pangea_bm_ok",
                                  high_cut = top_cut)
  sspn_pangea_no_bm <- time_dep_sspn(dat = data_validation, 
                                     marker_var = "pangea_no_bm_2_year_risk",
                                     obs_ok_var = "pangea_no_bm_ok",
                                     high_cut = top_cut)
  
  # mayo
  sspn_mayo <- time_dep_sspn(dat = data_validation, 
                             marker_var = "risk_rolling_20_2_20_prob_2_years",
                             obs_ok_var = "pangea_bm_ok",
                             high_cut = 0.4)
  
  # collect into a list
  time_sspn <- list(sspn_pangea_bm_traj = sspn_pangea_bm_traj,
                    sspn_pangea_no_bm_traj = sspn_pangea_no_bm_traj,
                    sspn_pangea_bm = sspn_pangea_bm,
                    sspn_pangea_no_bm = sspn_pangea_no_bm,
                    sspn_pangea_bm_traj_alt = sspn_pangea_bm_traj_alt,
                    sspn_pangea_no_bm_traj_alt = sspn_pangea_no_bm_traj_alt,
                    sspn_mayo = sspn_mayo)
  
  ## get rates of progression within 2 years for each risk strata by cohort_calib (mayo and pangea)
  
  # mayo
  mayo_high_calib <- cohort_calib(data_validation[rolling_20_2_20 == "High"])
  mayo_int_calib <- cohort_calib(data_validation[rolling_20_2_20 == "Intermediate"])
  mayo_low_calib <- cohort_calib(data_validation[rolling_20_2_20 == "Low"])
  
  # PANGEA trajectory models
  pangea_bm_traj_high_calib <- cohort_calib(data_validation[pangea_bm_traj_risk_cat == "High"])
  pangea_bm_traj_int_calib <- cohort_calib(data_validation[pangea_bm_traj_risk_cat == "Intermediate"])
  pangea_bm_traj_low_calib <- cohort_calib(data_validation[pangea_bm_traj_risk_cat == "Low"])
  
  pangea_no_bm_traj_high_calib <- cohort_calib(data_validation[pangea_no_bm_traj_risk_cat == "High"])
  pangea_no_bm_traj_int_calib <- cohort_calib(data_validation[pangea_no_bm_traj_risk_cat == "Intermediate"])
  pangea_no_bm_traj_low_calib <- cohort_calib(data_validation[pangea_no_bm_traj_risk_cat == "Low"])
  
  # PANGEA trajectory alt models
  if (include_pangea_alt == T) {
    pangea_bm_traj_alt_high_calib <- cohort_calib(data_validation[pangea_bm_traj_alt_risk_cat == "High"])
    pangea_bm_traj_alt_int_calib <- cohort_calib(data_validation[pangea_bm_traj_alt_risk_cat == "Intermediate"])
    pangea_bm_traj_alt_low_calib <- cohort_calib(data_validation[pangea_bm_traj_alt_risk_cat == "Low"])
    
    pangea_no_bm_traj_alt_high_calib <- cohort_calib(data_validation[pangea_no_bm_traj_alt_risk_cat == "High"])
    pangea_no_bm_traj_alt_int_calib <- cohort_calib(data_validation[pangea_no_bm_traj_alt_risk_cat == "Intermediate"])
    pangea_no_bm_traj_alt_low_calib <- cohort_calib(data_validation[pangea_no_bm_traj_alt_risk_cat == "Low"])
  }
  
  # PANGEA 1.0 models
  pangea_bm_high_calib <- cohort_calib(data_validation[pangea_bm_risk_cat == "High"])
  pangea_bm_int_calib <- cohort_calib(data_validation[pangea_bm_risk_cat == "Intermediate"])
  pangea_bm_low_calib <- cohort_calib(data_validation[pangea_bm_risk_cat == "Low"])
  
  pangea_no_bm_high_calib <- cohort_calib(data_validation[pangea_no_bm_risk_cat == "High"])
  pangea_no_bm_int_calib <- cohort_calib(data_validation[pangea_no_bm_risk_cat == "Intermediate"])
  pangea_no_bm_low_calib <- cohort_calib(data_validation[pangea_no_bm_risk_cat == "Low"])
  
  # put into a matrix
  risk_strat_calib <- cbind(mayo_high_calib = mayo_high_calib, 
                           mayo_int_calib = mayo_int_calib, 
                           mayo_low_calib = mayo_low_calib,
                           pangea_bm_traj_high_calib = pangea_bm_traj_high_calib,
                           pangea_bm_traj_int_calib = pangea_bm_traj_int_calib,
                           pangea_bm_traj_low_calib = pangea_bm_traj_low_calib,
                           pangea_no_bm_traj_high_calib = pangea_no_bm_traj_high_calib,
                           pangea_no_bm_traj_int_calib = pangea_no_bm_traj_int_calib,
                           pangea_no_bm_traj_low_calib = pangea_no_bm_traj_low_calib,
                           pangea_bm_traj_alt_high_calib = if(include_pangea_alt) pangea_bm_traj_alt_high_calib,
                           pangea_bm_traj_alt_int_calib = if(include_pangea_alt) pangea_bm_traj_alt_int_calib,
                           pangea_bm_traj_alt_low_calib = if(include_pangea_alt) pangea_bm_traj_alt_low_calib,
                           pangea_no_bm_traj_alt_high_calib = if(include_pangea_alt) pangea_no_bm_traj_alt_high_calib,
                           pangea_no_bm_traj_alt_int_calib = if(include_pangea_alt) pangea_no_bm_traj_alt_int_calib,
                           pangea_no_bm_traj_alt_low_calib = if(include_pangea_alt) pangea_no_bm_traj_alt_low_calib,
                           pangea_bm_high_calib = pangea_bm_high_calib,
                           pangea_bm_int_calib = pangea_bm_int_calib,
                           pangea_bm_low_calib = pangea_bm_low_calib,
                           pangea_no_bm_high_calib = pangea_no_bm_high_calib,
                           pangea_no_bm_int_calib = pangea_no_bm_int_calib,
                           pangea_no_bm_low_calib = pangea_no_bm_low_calib)
  
  
  ## get rates of progression within 2 years for each risk strata by cohort_calib (mayo and pangea)

  # pangea trajectory models
  pangea_bm_traj_cat_km <- risks_profile(data_validation, return_surv = save_deidentified,
                                         risk_cat_varname = "pangea_bm_traj_risk_cat",
                                         model_name = "PANGEA Trajectory (BM)")
  pangea_no_bm_traj_cat_km <- risks_profile(data_validation, return_surv = save_deidentified,
                                            risk_cat_varname = "pangea_no_bm_traj_risk_cat",
                                            model_name = "PANGEA Trajectory (No BM)")
  
  # pangea trajectory alt models
  if (include_pangea_alt == T) {
    pangea_bm_traj_alt_cat_km <- risks_profile(data_validation, return_surv = save_deidentified,
                                           risk_cat_varname = "pangea_bm_traj_alt_risk_cat",
                                           model_name = "PANGEA Trajectory (BM)")
    pangea_no_bm_traj_alt_cat_km <- risks_profile(data_validation, return_surv = save_deidentified,
                                              risk_cat_varname = "pangea_no_bm_traj_alt_risk_cat",
                                              model_name = "PANGEA Trajectory (No BM)")
  }
  
  # mayo
  mayo_cat_km <- risks_profile(data_validation, return_surv = save_deidentified,
                               risk_cat_varname = "rolling_20_2_20",
                               model_name = "20/2/20")
  
  # pangea 1.0 models
  pangea_bm_cat_km <- risks_profile(data_validation, return_surv = save_deidentified,
                                    risk_cat_varname = "pangea_bm_risk_cat",
                                    model_name = "PANGEA 1.0 (BM)")
  pangea_no_bm_cat_km <- risks_profile(data_validation, return_surv = save_deidentified,
                                       risk_cat_varname = "pangea_no_bm_risk_cat",
                                       model_name = "PANGEA 1.0 (No BM)")
  
  # put into list
  risk_strat_km <- list(pangea_bm_traj_cat_km = pangea_bm_traj_cat_km,
                        pangea_no_bm_traj_cat_km = pangea_no_bm_traj_cat_km,
                        pangea_bm_traj_alt_cat_km = if(include_pangea_alt) pangea_bm_traj_alt_cat_km,
                        pangea_no_bm_traj_alt_cat_km = if(include_pangea_alt) pangea_no_bm_traj_alt_cat_km,
                        mayo_cat_km = mayo_cat_km,
                        pangea_bm_cat_km = pangea_bm_cat_km,
                        pangea_no_bm_cat_km = pangea_no_bm_cat_km)
  
  ## average time at which progressors/nonprogressors are first classified as high risk
  data_validation[, first_pangea_bm_traj_high := (pangea_bm_traj_risk_cat == "High" & 
                                                    !duplicated(pangea_bm_traj_risk_cat)), 
                  by = participant_id]
  data_validation[, first_pangea_bm_high := (pangea_bm_risk_cat == "High" & 
                                               !duplicated(pangea_bm_risk_cat)), 
                  by = participant_id]
  data_validation[, first_mayo_high := (rolling_20_2_20 == "High" & 
                                          !duplicated(rolling_20_2_20)), 
                  by = participant_id]
  if (include_pangea_alt == T) {
    data_validation[, first_pangea_bm_traj_alt_high := (pangea_bm_traj_alt_risk_cat == "High" & 
                                                          !duplicated(pangea_bm_traj_alt_risk_cat)), 
                    by = participant_id]
  }
  
  data_validation[, first_pangea_no_bm_traj_high := (pangea_no_bm_traj_risk_cat == "High" & 
                                                    !duplicated(pangea_no_bm_traj_risk_cat)), 
                  by = participant_id]
  data_validation[, first_pangea_no_bm_high := (pangea_no_bm_risk_cat == "High" & 
                                               !duplicated(pangea_no_bm_risk_cat)), 
                  by = participant_id]
  if (include_pangea_alt == T) {
    data_validation[, first_pangea_no_bm_traj_alt_high := (pangea_no_bm_traj_alt_risk_cat == "High" & 
                                                          !duplicated(pangea_no_bm_traj_alt_risk_cat)), 
                    by = participant_id]
  }
  
  # progressors
  first_high_prog <- c(mayo = data_validation[ever_mm == 1 & first_mayo_high == T, 
                                              mean(last_date - tstart)],
                       pangea_bm_traj = data_validation[ever_mm == 1 & first_pangea_bm_traj_high == T, 
                                                        mean(last_date - tstart)],
                       pangea_bm_traj_alt = if(include_pangea_alt) data_validation[ever_mm == 1 & first_pangea_bm_traj_alt_high == T, 
                                                        mean(last_date - tstart)],
                       pangea_bm = data_validation[ever_mm == 1 & first_pangea_bm_high == T, 
                                                   mean(last_date - tstart)],
                       pangea_no_bm_traj = data_validation[ever_mm == 1 & first_pangea_no_bm_traj_high == T, 
                                                        mean(last_date - tstart)],
                       pangea_no_bm_traj_alt = if(include_pangea_alt) data_validation[ever_mm == 1 & first_pangea_no_bm_traj_alt_high == T, 
                                                                                   mean(last_date - tstart)],
                       pangea_no_bm = data_validation[ever_mm == 1 & first_pangea_no_bm_high == T, 
                                                   mean(last_date - tstart)])
  
  # non-progressors
  first_high_nonprog <- c(mayo = data_validation[ever_mm == 0 & first_mayo_high == T, 
                                              mean(last_date - tstart)],
                       pangea_bm_traj = data_validation[ever_mm == 0 & first_pangea_bm_traj_high == T, 
                                                        mean(last_date - tstart)],
                       pangea_bm_traj_alt = if(include_pangea_alt) data_validation[ever_mm == 0 & first_pangea_bm_traj_alt_high == T, 
                                                                                   mean(last_date - tstart)],
                       pangea_bm = data_validation[ever_mm == 0 & first_pangea_bm_high == T, 
                                                   mean(last_date - tstart)],
                       pangea_no_bm_traj = data_validation[ever_mm == 0 & first_pangea_no_bm_traj_high == T, 
                                                           mean(last_date - tstart)],
                       pangea_no_bm_traj_alt = if(include_pangea_alt) data_validation[ever_mm == 0 & first_pangea_no_bm_traj_alt_high == T, 
                                                                                      mean(last_date - tstart)],
                       pangea_no_bm = data_validation[ever_mm == 0 & first_pangea_no_bm_high == T, 
                                                      mean(last_date - tstart)])
  
  
  
  # put in matrix
  first_high_times <- cbind(first_high_prog = first_high_prog,
                            first_high_nonprog = first_high_nonprog)
  
  
  # return output
  list(reclass_tables = reclass_tables,
       ever_high = ever_high,
       time_sspn = time_sspn,
       risk_strat_calib = risk_strat_calib,
       risk_strat_km = risk_strat_km,
       first_high_times = first_high_times)
}

# reclassification tables where the pangea strata are based on clinical thresholds (low < 20%, high > 40%)
cutoffs_clinical <- c(0, 0.1, 0.4, 1)
pangea_labels <- c("Low", "Intermediate", "High")

reclass_results <- pangea_reclassification(cutoffs_pangea_bm_traj = cutoffs_clinical, 
                                           cutoffs_pangea_no_bm_traj = cutoffs_clinical,
                                           cutoffs_pangea_bm = cutoffs_clinical, 
                                           cutoffs_pangea_no_bm = cutoffs_clinical,
                                           pangea_labels = pangea_labels,
                                           include_pangea_alt = T)

