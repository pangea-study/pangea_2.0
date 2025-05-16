### code to analyze FISH data (only for cohorts where it's available)

## load packages/functions and set-up data
source("Packages.R")
source("Functions.R")

path_data <- "validation_data.csv" # Put here the path to 
                                   # your data (see README.pdf for the format)

data_validation <- fread(path_data)
data_validation <- data_validation[current_diagnosis == "SMM"]

source("data_preparation.R")

load("models/pangea_bm_traj_alt.RData", verbose = T)


data_validation[, pangea_bm_traj_ok := !(is.na(plasmacells) | 
                                           is.na(mspike) | 
                                           is.na(iuratio) | 
                                           is.na(creatinine) | 
                                           is.na(age)|
                                           is.na(hgb)|
                                           is.na(mspike_traj)|
                                           is.na(iuratio_traj)|
                                           is.na(hgb_traj)|
                                           is.na(creatinine_traj))]

data_validation[, pangea_risk_cat := ifelse(pangea_bm_traj_alt_2_year_risk < 0.1, "Low",
                                            ifelse(pangea_bm_traj_alt_2_year_risk < 0.4, "Intermediate", 
                                                   "High"))]



## List of probes available at Heidelberg:
# t(4;14)
# t(14;16)
# Gain(1q)
# Del17p
# 

### Combination of probes to test 
## (Modified to account for the fact that t_14_20 and del1p are unavailable for Heidelberg)

data_validation[,
                comp_fish_1 := fifelse(
                  any(is.na(t_4_14), is.na(t_14_16), is.na(gain_1q)) &        # any marker NA
                    sum(t_4_14, t_14_16,  gain_1q, na.rm = T) == 0,           # no markers positive
                  NA,
                  fifelse(
                    sum(t_4_14, t_14_16, na.rm = T)  > 0 & sum(gain_1q, na.rm = T) > 0,
                    1, 0
                  )
                ),
                by = 1:nrow(data_validation)                   
]

## any FISH positive
data_validation[, any_FISH := as.numeric((t_4_14 | t_14_16 | del17p | gain_1q ))]



### We can test (only partially) only the first combination of probe


# Function that, given a specific probe (or the name of the combination
# variable) return the info we need

probe_analysis <- function(probe, model = model_traj_bm_alt,
                           data = df_train){
  
  data = data.table(data)
  
  ## Reference concordance value
  data[,
       model_ok_probe := pangea_bm_traj_ok & !is.na(.SD),
       .SDcols = probe]
  data_red <- data[  model_ok_probe == 1,,]
  c_reference <- concordance(model, newdata = data_red)$concordance
  
  ## re-estimate model on data_red and get c-statistic
  model_red <- coxph(model$formula, data_red)
  c_model_red_reference <- concordance(model_red)$concordance
  
  # Number of patients for which the probe has been tested at least once
  n_patients <- data_red[, sum(!duplicated(participant_id))]
  
  # Number of patient presenting the alteration
  n_positive <- sum(data_red[, any(get(probe) == 1), 
                             by = participant_id]$V1)
  n_prog_positive <- sum(data_red[ever_mm == T, any(get(probe) == 1), 
                                  by = participant_id]$V1)
  
  # pangea risk strata vs. FISH probe status
  pangea_by_probe_table <- data[, table(pangea_risk_cat, get(probe), 
                                        useNA = "always")]
  
  # survival by probe/pangea status
  km_probe_only <- survfit(as.formula(paste0("Surv(tstart, tstop, prog_mm) ~ ", probe)),
                           data = data)
  km_probe_pangea <- survfit(as.formula(paste0("Surv(tstart, tstop, prog_mm) ~ pangea_risk_cat +", probe)),
                             data = data)
  km_pangea_only <- survfit(as.formula(paste0("Surv(tstart, tstop, prog_mm) ~ pangea_risk_cat")),
                            data = data)
  
  ## P-value and C-statistic of the model with the probe only
  probe_only_formula <- as.formula(paste0("Surv(tstart, tstop, prog_mm)~", probe))
  
  probe_only_mod <- coxph(probe_only_formula, data = data)
  probe_only_p_val <- summary(probe_only_mod)$coefficients[probe, "Pr(>|z|)"]
  probe_only_c <- concordance(probe_only_mod)$concordance
  
  info_probe_only_model <- summary(probe_only_mod)
  info_probe_only_coefficients <- summary(probe_only_mod)$coefficients
  hr_probe_only <- summary(probe_only_mod)$coefficients[probe, "exp(coef)"]
  
  ## P-value and C-statistic of the model with the probe and PANGEA 2.0 (BM)
  pangea_probe_formula <- as.formula(paste0("Surv(tstart, tstop, prog_mm)~mspike + 
  log_iuratio + plasmacells + log_creatinine + age +
        hgb_traj+ iuratio_traj + mspike_traj + creatinine_traj + ", probe))
  
  pangea_probe_mod <- coxph(pangea_probe_formula, data = data)
  pangea_probe_p_val <- summary(pangea_probe_mod)$coefficients[probe, "Pr(>|z|)"]
  pangea_probe_c <- concordance(pangea_probe_mod)$concordance
  
  info_pangea_probe_model <- summary(pangea_probe_mod)
  info_pangea_probe_coefficients <- summary(pangea_probe_mod)$coefficients
  hr_pangea_probe <- summary(pangea_probe_mod)$coefficients[probe, "exp(coef)"]
  
  results_row <- cbind(n_patients = n_patients, n_positive = n_positive, 
                       n_prog_positive = n_prog_positive, 
                       hr_probe_only = hr_probe_only,
                       probe_only_p_val = probe_only_p_val, 
                       hr_pangea_probe = hr_pangea_probe,
                       pangea_probe_p_val = pangea_probe_p_val,
                       pangea_probe_c = pangea_probe_c, 
                       c_model_red_reference = c_model_red_reference, 
                       c_reference = c_reference)
  
  
  return(list(results = results_row, 
              info_probe_only_model = info_probe_only_model,
              info_pangea_probe_model = info_pangea_probe_model,
              pangea_by_probe_table = pangea_by_probe_table,
              km_probe_only = km_probe_only,
              km_probe_pangea = km_probe_pangea,
              km_pangea_only = km_pangea_only))
}


## get results for each probe
list_fish_variables <- list(
  "t_4_14",
  "t_14_16",
  "del17p",
  "gain_1q",
  "comp_fish_1",
  "any_FISH"
)
names(list_fish_variables) <- list_fish_variables

list_fish_res <- lapply(list_fish_variables, probe_analysis,
                        data = data_validation,
                        model = model_traj_bm_alt)
names(list_fish_res) <- list_fish_variables

save(list_fish_res, file = "results/fish_results.RData")

