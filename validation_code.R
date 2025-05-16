#### Set pangea validation as working directory
path_data <- "validation_data.csv" # Put here the path to 
                                   # your data (see README.pdf for the format)

#### Set up data sharing variable
#### CHANGE THIS IF DESIRED
save_deidentified <- TRUE
   # if set to FALSE, no patient-level data will be saved into the results
   # if set to TRUE, some de-identified patient-level data will be saved.
   #   (ONLY: the measurement times, progression/censoring status, and the various
   #   risk scores. no patient biomarker or demographic data will be saved)
   # sharing the de-identified data is not necessary, but it would be helpful as 
   # it would give us more flexibility to modify some of the figures and tables.

#### Loading the packages and the functions
source("Packages.R")

source("Functions.R")

#### Loading the data
data_validation <- fread(path_data)
# ensure we only have smoldering patients
data_validation <- data_validation[current_diagnosis == "SMM"]

## load the models
for(f in list.files("models", full.names = T)) {
  load(f, verbose = T)
}


#### Compute the trajectories, risk scores, and other variables needed for validation
source("data_preparation.R")  # could take a few seconds
# this may produce warnings like "In summary.lm(curr_mod) : essentially perfect fit: summary may be unreliable"
   # that is expected and not an error


#### Compute descriptive statistics about the dataset
## summaries for Table 1
source("table_descriptives.R")

save(table_numeric_info, list_cat_info,
     file = "results/descriptive_tables.rda")

## summaries for Table A1 (median and IQR of biomarkers) 
## and A2 (frequencies and HR of trajectories, cohort-specific)
source("supplementary_biomarker_tables.R")

save(supp_biomarker_list, 
     data_validation_pangea_bm_traj_mod,
     data_validation_pangea_no_bm_traj_mod,
     file = "results/supplementary_biomarker_info.rda")

## Kaplan-Meier survival curves + 2-year progression rates (stratified by 20/2/20)
risk_cat_varname <- "rolling_20_2_20"
mayo_pfs_surv <- risks_profile(data_validation, return_surv = save_deidentified,
                               whole_cohort = F,
                               risk_cat_varname = "rolling_20_2_20",
                               model_name = "20/2/20")

save(mayo_pfs_surv, 
     file = "results/mayo_pfs_surv.rda")


## Kaplan-Meier survival curve + 2-year progression rate for whole cohort
data_validation[, cohort_dummy := "Overall"]
overall_pfs_surv <- risks_profile(data_validation, return_surv = save_deidentified,
                                  whole_cohort = T,
                                  risk_cat_varname = "cohort_dummy",
                                  model_name = "20/2/20")

save(overall_pfs_surv, 
     file = "results/overall_pfs_surv.rda")


#### Compute some descriptive statistics about the PANGEA risk predictions
## histograms of difference in risk predictions for various pairs of models
validation_risk_diff <- risk_diff_plot(data_validation)
validation_risk_diff_alt <- risk_diff_plot(data_validation,
                                           pangea_bm_traj_varname = "pangea_bm_traj_alt_2_year_risk",
                                           pangea_no_bm_traj_varname = "pangea_no_bm_traj_alt_2_year_risk")

save(validation_risk_diff,
     file = "results/validation_risk_diff.rda")
save(validation_risk_diff_alt,
     file = "results/validation_risk_diff_alt.rda")

#### Compute the concordance statistics (quality of ranking of risks)

## compute time-averaged C-statistics
# PANGEA trajectory models
c_stat_traj_bm <- concordance(model_traj_bm, newdata = data_validation)
c_stat_traj_no_bm <- concordance(model_traj_no_bm, newdata = data_validation)

# PANGEA trajectory models with no patient history
c_stat_traj_no_hist_bm <- concordance(model_traj_bm, newdata = data_validation_traj0)
c_stat_traj_no_hist_no_bm <- concordance(model_traj_no_bm, newdata = data_validation_traj0)

# PANGEA trajectory alt models
c_stat_traj_bm_alt <- concordance(model_traj_bm_alt, newdata = data_validation)
c_stat_traj_no_bm_alt <- concordance(model_traj_no_bm_alt, newdata = data_validation)

# PANGEA trajectory alt models with no patient history
c_stat_traj_no_hist_bm_alt <- concordance(model_traj_bm_alt, newdata = data_validation_traj0)
c_stat_traj_no_hist_no_bm_alt <- concordance(model_traj_no_bm_alt, newdata = data_validation_traj0)

# PANGEA 1.0 models
c_stat_bm <- concordance(mod_dep_time, newdata = data_validation)
c_stat_no_bm <- concordance(mod_dep_time_no_BM, newdata = data_validation)

# rolling 20/2/20
c_stat_20_2_20 <- concordance(mod_20_2_20, newdata = data_validation) 

## compute time-specific C-statistics
## (note this will produce NA results for any times at which the concordance
##  cannot be estimated since the dataset is too small)

# PANGEA trajectory models
c_stat_time_traj_bm <- c_stat_time_matrix(data_validation, model_traj_bm, 
                                          time_min = 0, time_max = 5, 
                                          time_increment = 0.1,
                                          mod_type = "pangea_bm")
c_stat_time_traj_no_bm <- c_stat_time_matrix(data_validation, model_traj_no_bm, 
                                          time_min = 0, time_max = 5, 
                                          time_increment = 0.1,
                                          mod_type = "pangea_no_bm")

# PANGEA trajectory models with no patient history
c_stat_time_traj_no_hist_bm <- c_stat_time_matrix(data_validation_traj0, 
                                                  model_traj_bm, 
                                                  time_min = 0, time_max = 5, 
                                                  time_increment = 0.1,
                                                  mod_type = "pangea_bm")
c_stat_time_traj_no_hist_no_bm <- c_stat_time_matrix(data_validation_traj0, 
                                                     model_traj_no_bm, 
                                                     time_min = 0, time_max = 5, 
                                                     time_increment = 0.1,
                                                     mod_type = "pangea_no_bm")

# PANGEA trajectory alt models
c_stat_time_traj_bm_alt <- c_stat_time_matrix(data_validation, model_traj_bm_alt, 
                                              time_min = 0, time_max = 5, 
                                              time_increment = 0.1,
                                              mod_type = "pangea_bm")
c_stat_time_traj_no_bm_alt <- c_stat_time_matrix(data_validation, model_traj_no_bm_alt, 
                                                 time_min = 0, time_max = 5, 
                                                 time_increment = 0.1,
                                                 mod_type = "pangea_no_bm")

# PANGEA trajectory alt models with no patient history
c_stat_time_traj_no_hist_bm_alt <- c_stat_time_matrix(data_validation_traj0, 
                                                      model_traj_bm_alt, 
                                                      time_min = 0, time_max = 5, 
                                                      time_increment = 0.1,
                                                      mod_type = "pangea_bm")
c_stat_time_traj_no_hist_no_bm_alt <- c_stat_time_matrix(data_validation_traj0, 
                                                         model_traj_no_bm_alt, 
                                                         time_min = 0, time_max = 5, 
                                                         time_increment = 0.1,
                                                         mod_type = "pangea_no_bm")

# PANGEA 1.0 models
c_stat_time_bm <- c_stat_time_matrix(data_validation, mod_dep_time, 
                                          time_min = 0, time_max = 5, 
                                          time_increment = 0.1,
                                          mod_type = "pangea_bm")
c_stat_time_no_bm <- c_stat_time_matrix(data_validation, mod_dep_time_no_BM, 
                                             time_min = 0, time_max = 5, 
                                             time_increment = 0.1,
                                             mod_type = "pangea_no_bm")

# rolling 20/2/20
c_stat_time_20_2_20 <- c_stat_time_matrix(data_validation, mod_20_2_20, 
                                          time_min = 0, time_max = 5, 
                                          time_increment = 0.1, 
                                          mod_type = "pangea_bm")

## save the c-statistics
c_stats_time_averaged <- list(c_stat_traj_bm = c_stat_traj_bm, 
                              c_stat_traj_no_bm = c_stat_traj_no_bm,
                              c_stat_traj_no_hist_bm = c_stat_traj_no_hist_bm, 
                              c_stat_traj_no_hist_no_bm = c_stat_traj_no_hist_no_bm,
                              c_stat_traj_bm_alt = c_stat_traj_bm_alt, 
                              c_stat_traj_no_bm_alt = c_stat_traj_no_bm_alt,
                              c_stat_traj_no_hist_bm_alt = c_stat_traj_no_hist_bm_alt, 
                              c_stat_traj_no_hist_no_bm_alt = c_stat_traj_no_hist_no_bm_alt,
                              c_stat_bm = c_stat_bm, 
                              c_stat_no_bm = c_stat_no_bm,
                              c_stat_20_2_20 = c_stat_20_2_20)

c_stats_time_specific <- list(c_stat_time_traj_bm = c_stat_time_traj_bm,
                              c_stat_time_traj_no_bm = c_stat_time_traj_no_bm,
                              c_stat_time_traj_no_hist_bm = c_stat_time_traj_no_hist_bm,
                              c_stat_time_traj_no_hist_no_bm = c_stat_time_traj_no_hist_no_bm,
                              c_stat_time_traj_bm_alt = c_stat_time_traj_bm_alt,
                              c_stat_time_traj_no_bm_alt = c_stat_time_traj_no_bm_alt,
                              c_stat_time_traj_no_hist_bm_alt = c_stat_time_traj_no_hist_bm_alt,
                              c_stat_time_traj_no_hist_no_bm_alt = c_stat_time_traj_no_hist_no_bm_alt,
                              c_stat_time_bm = c_stat_time_bm,
                              c_stat_time_no_bm = c_stat_time_no_bm,
                              c_stat_time_20_2_20 = c_stat_time_20_2_20)

save(c_stats_time_averaged, c_stats_time_specific, 
     file = "results/c_stats.rda")


#### Compute the risk categorization
source("risk_reclassification.R")

save(reclass_results,
     file = "results/reclass_results.rda")

#### Compute the calibration (comparison between predicted and actual risk)
# note that NAs will be produced for subcohorts which empty/too small
validation_calibration <- calibration_stats(data_validation,
                                            pred_varnames = c("pangea_bm_traj_2_year_risk",
                                                              "pangea_no_bm_traj_2_year_risk",
                                                              "pangea_bm_traj_alt_2_year_risk",
                                                              "pangea_no_bm_traj_alt_2_year_risk",
                                                              "pangea_bm_2_year_risk",
                                                              "pangea_no_bm_2_year_risk",
                                                              "risk_rolling_20_2_20_prob_2_years"))

save(validation_calibration, 
     file = "results/validation_calibration.rda")

#### If accepted above (save_deidentified == TRUE), make limited version of the 
#### dataset with just de-identified patient level survival and risk predictions 
#### (no biomarkers/demographics)
if (save_deidentified == TRUE) {
  remove_vars <- c("mspike", "iuratio", "plasmacells", "creatinine",
                   "sex", "race", "ethnicity", "age", "hgb", "immunofix2", "log_iuratio", "log_creatinine",
                   "mspike_traj", "creatinine_traj", "iuratio_traj", "hgb_traj",
                   "traj_hgb")
  
  data_validation_limited <- data.table(data_validation)
  data_validation_limited[, c(remove_vars) := NULL]
  
  write.csv(data_validation_limited, "results/data_validation_limited.csv")
}

