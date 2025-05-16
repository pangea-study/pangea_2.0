
# Function to find the first observation in the previous X months
find_index_start <- function(i, interval_range, follow_ups, labs_val, prev_obs = T, silent = T){
  if (silent == F) {
    print(i)
  }
  cumulative_times <- rev(cumsum(follow_ups[(i-1):1]))
  if (silent == F) {
    print(cumulative_times)
  }
  no_na_labsval <- which(!is.na(labs_val[1:(i-1)]))
  if (silent == F) {
    print(no_na_labsval)
  }
  cumulative_times <- cumulative_times[!is.na(labs_val[1:(i-1)])]
  if (silent == F) {
    print(cumulative_times)
  }
  if(any(cumulative_times <= interval_range)){
    index <- which.max(cumulative_times <= interval_range)
  }
  else{
    if(prev_obs == T){
      index <- i-1
    }
    else{
      return(NA)
    }
  }
  return(no_na_labsval[index])
}



##### Alternative evaluation strategy: we add the trajectories 
## Then, we adopt k-fold cross validation to select the optimal one
## This function divide the participant ids in k groups. Then, we will use
## k-1 sets to estimate the cox model and 1 to validate the estimated model

k_fold <- function(participant_id, k){
  participant_id <- unique(participant_id) # Just to be sure we call the function correctly
  x <- runif(length(participant_id))
  shuffled_participant <- participant_id[rank(x)]
  k_fold_sets <- split(shuffled_participant, cut(seq_along(shuffled_participant), 
                                                 k, labels = FALSE))
  return(k_fold_sets)
}

##### Arguments
### - training dataset
### - validation dataset
### - the name of the trajectory variable
c_statistic_traj <- function(data_train, data_val, traj = "traj", bm = T, tstart = "tstart_mm_or_early",
                             tstop = "tstop_mm_or_early", prog_mm = "prog_mm_or_early"){
  if(bm == T){
    formula <- as.formula(paste0("Surv(", tstart, ",", tstop, ",", prog_mm,") ~ mspike + log_iuratio + 
                   log_creatinine + plasmacells + age + noquote(traj)"))
  }
  else{
    formula <- as.formula(paste0("Surv(", tstart, ",", tstop, ",", prog_mm,") ~ mspike + log_iuratio + 
                   log_creatinine + age + noquote(traj)"))
    
  }
  
  mod <- tryCatch(coxph(formula, data = data_train), error=function(err) NULL)
  
  return(tryCatch(concordance(mod , newdata = data_val,
                              timewt = "n")$concordance, error=function(err) NA))
}


#### Combining k-fold and c_statistic_traj to return the mean c statistic
k_fold_validation <- function(data, traj = "traj", k = 5, bm = T, tstart = "tstart_mm_or_early",
                              tstop = "tstop_mm_or_early", prog_mm = "prog_mm_or_early"){
  
  folders <- k_fold(unique(data$participant_id), k = k)
  c_stat <- sapply(1:k, function(i){
    indeces_train <- (1:k)[-i]
    data_train <- data[participant_id %in% unlist(folders[indeces_train]),]
    data_val <- data[participant_id %in% folders[[i]], ]
    return(c_statistic_traj(data_train = data_train,
                            data_val = data_val,
                            traj = traj, bm = bm,
                            tstart = tstart, tstop = tstop, prog_mm = prog_mm))
  })
  return(mean(c_stat))
}

## The results are variable, I perform the k-fold R times and take
## the mean
repeated_k_fold_validation <- function(data, traj = "traj", k = 5, bm = T, tstart = "tstart_mm_or_early",
                                       tstop = "tstop_mm_or_early", prog_mm = "prog_mm_or_early", R = 50){
  
  c_stat_rep <- sapply(1:R, function(r){
    return(k_fold_validation(data = data, traj = traj, k = k, bm = bm, tstart = tstart,
                             tstop = tstop, prog_mm = prog_mm))
  })
  return(mean(c_stat_rep))
}



##### General CV C-statistics (I adapt the one for the traj to work in general)


c_statistic <- function(data_train, data_val, formula){
  
  mod <- tryCatch(coxph(formula, data = data_train), error=function(err) NULL)
  
  return(tryCatch(concordance(mod , newdata = data_val,
                              timewt = "n")$concordance, error=function(err) NA))
}


#### Combining k-fold and c_statistic_traj to return the mean c statistic
k_fold_validation_general <- function(data, k = 5, formula){
  
  folders <- k_fold(unique(data$participant_id), k = k)
  c_stat <- sapply(1:k, function(i){
    indeces_train <- (1:k)[-i]
    data_train <- data[participant_id %in% unlist(folders[indeces_train]),]
    data_val <- data[participant_id %in% folders[[i]], ]
    return(c_statistic(data_train = data_train,
                       data_val = data_val,
                       formula = formula
    ))
  })
  return(mean(c_stat))
}

## The results are variable, I perform the k-fold R times and take
## the mean
repeated_k_fold_validation_general <- function(data, formula, k = 5, R = 50){
  
  c_stat_rep <- sapply(1:R, function(r){
    return(k_fold_validation_general(data = data, formula = formula, k = k))
  })
  return(mean(c_stat_rep))
  
}



### Function to compute risks at X years given a model
get_risk_prob <- function(fit_prob, fit_time, year = 2){
  
  min_time <- min(abs(fit_time - year))
  ind_years <- which(abs(fit_time - year) == min_time)
  surv_years <- fit_prob[ind_years]
  
  return(round(1 - surv_years, 3))
  
}

# function_trajectories
traj_f <- function(values,ids,timediff, treshold_pvalue = 1){
  
  df <- cbind(values, ids, timediff) %>% as_tibble() %>% 
    arrange(ids, timediff)
  traj <- numeric(nrow(df))
  
  unique_ids <- unique(ids)
  for(i in unique_ids){
    
    ind_curr <- which(df$ids == i)
    df_curr <- df[ind_curr,]
    
    if(length(ind_curr) < 3){
      traj[ind_curr] = 0
    }else{
      
      traj[ind_curr[1:2]] = 0
      for(h in 3:length(ind_curr)){
        
        if(is.na(df_curr$values[h])){
          traj[ind_curr[h]] <- NA
        }else{
          curr_mod <- lm(values ~ timediff, data = df_curr[1:h,])
          
          if(!is.na(curr_mod$coefficients[2])){
            test_p <- summary(curr_mod)$coefficients[2,4] 
            if(!is.na(test_p)){
              if(test_p < treshold_pvalue){
                traj[ind_curr[h]] <- curr_mod$coefficients[2]
              }
            }
          }
          
          
        }
        
      }
    }
  }
  
  return(traj)
  
}


reduced.dataset <- function(data){
  
  data <- data[pangea_bm_ok == T,]
  
  data[, obs_id := order(tstart), by = participant_id]
  data[, last_obs := as.numeric(obs_id == max(obs_id)), by = participant_id]
  data[, tstop_mm_or_early := ifelse(last_obs == 0,tstop_mm_or_early, last_date),]
  
  data[, prog_mm_or_early := ifelse(last_obs == 0,prog_mm_or_early, ever_mm_or_early),
  ]
  return(data)
}


#### Which data do I want? Function to load the data of interest based on the
## following criteria
# - Cohort: train (DFCI) or val (UK + Greece)
# - Population: all, smm, mgus
# - Na carried forward: T/F
# - Progression criterion: redcap/slimcrab
# It is build so that the tstart, tstop, last_date and prog_mm variables
# contains the information we need

load_data <- function(cohort = "train", population = "all", 
                      na_carrying_forward = T, prog_criterion = "slimcrab",
                      data_directory = "../New_datasets",
                      dataset = "df_all_cohort.csv"){
  
  
  
  data_path <- paste(data_directory, "/", dataset, sep = "")
  dataset <- data.table::fread(data_path)
  
  
  if(prog_criterion == "slimcrab"){
    dataset <- dataset[is_lab_at_or_after_early_prog == 0,]
    dataset[,c("tstart","tstop","prog_mm","last_date","ever_mm") := 
              list(tstart_mm_or_early,
                   tstop_mm_or_early,
                   prog_mm_or_early,
                   last_date_mm_or_early,
                   ever_mm_or_early)]
  }
  
  if(population == "smm"){
    dataset <- dataset[current_diagnosis == "SMM",]
    dataset[,c("tstart", "tstop", "last_date") 
            := list(
              tstart_from_smm, tstop_from_smm, last_date_smm
            )]
  }
  
  if(population == "mgus"){
    dataset <-  dataset[current_diagnosis == "MGUS",]
  }
  
  if(cohort == "train"){
    dataset <- dataset[pi_cohort == "DFCI",]
  }
  else{
    dataset <- dataset[pi_cohort %in% c("UK", "Greece"),]
  }
  
  # Adjustment auxiliary variable
  dataset[, first_obs:= as.numeric(tstart == 0)]
  dataset[, last_obs:= as.numeric(tstop == last_date)]
  dataset[, time_to_mm_or_cens := last_date - tstart]
  
  return(dataset)
}

### Function that, given a dataset, computes the optimal trajectories
# The optimal trajectories are
# - M-spike: 0.2 increase with respect to any observation in the last 18 months
# - HgB: 1.5 decrease with respect to any observation in the last year
# - Iuratio: 20 increase with respect to any observation in the last 2 years
# - Creatinine: 25% increase with respect to any observation in the last year

computation_trajectories <- function(dataset, silent = T,
                                     traj_hgb_p_thresh = 0.1){
  dataset <- data.table(dataset)    # so changes are not made in place
  
  ## Trajectories of the Pangea Traj model
  dataset[!is.na(mspike),  mspike_traj:= sapply(1:length(mspike), 
                                                function(i, 
                                                         absolute_inc = 0.2, 
                                                         interval_range = 1.5){
                                                  if(i < 2){
                                                    return(0)
                                                  }
                                                  
                                                  i_start <- find_index_start(i = i, 
                                                                              interval_range = interval_range, 
                                                                              follow_ups = time_between_obs, 
                                                                              labs_val = mspike,
                                                                              prev_obs = F)
                                                  
                                                  if(is.na(i_start))
                                                    return(0)
                                                  
                                                  increase <- (mspike[i] - mspike[i_start:(i-1)])
                                                  
                                                  if(any(is.na(increase))){
                                                    increase <- increase[-which(is.na(increase))]
                                                  }
                                                  
                                                  if(any(increase >= absolute_inc)){
                                                    return(1)
                                                  }else{
                                                    return(0)
                                                  }
                                                }), 
          by = participant_id]
  
  
  
  dataset[!(is.na(creatinine)), 
          creatinine_traj :=  sapply(1:length(creatinine), 
                                     function(i, 
                                              percentage = 0.25, 
                                              interval_range = 1){
                                       if(i < 2){
                                         return(0)
                                       }
                                       
                                       i_start <- find_index_start(i = i, 
                                                                   interval_range = interval_range, 
                                                                   follow_ups = time_between_obs, 
                                                                   labs_val = creatinine,
                                                                   prev_obs = F)
                                       
                                       if(is.na(i_start))
                                         return(0)
                                       
                                       increase <- (creatinine[i] - creatinine[i_start:(i-1)])
                                       increase <- increase/creatinine[i_start:(i-1)]
                                       
                                       if(any(is.na(increase))){
                                         increase <- increase[-which(is.na(increase))]
                                       }
                                       
                                       if(any(increase >= percentage)){
                                         return(1)
                                       }else{
                                         return(0)
                                       }
                                     }), 
          by = participant_id]
  
  
  dataset[!is.na(hgb),  hgb_traj := sapply(1:length(hgb), 
                                           function(i, 
                                                    absolute_inc = 1.5 , 
                                                    interval_range = 1){
                                             if(i < 2){
                                               return(0)
                                             }
                                             
                                             i_start <- find_index_start(i = i, 
                                                                         interval_range = interval_range, 
                                                                         follow_ups = time_between_obs, 
                                                                         labs_val = hgb,
                                                                         prev_obs = F)
                                             
                                             if(is.na(i_start))
                                               return(0)
                                             
                                             increase <- (hgb[i] - hgb[i_start:(i-1)])
                                             
                                             if(any(is.na(increase))){
                                               increase <- increase[-which(is.na(increase))]
                                             }
                                             
                                             if(any(increase <= - absolute_inc)){
                                               return(1)
                                             }else{
                                               return(0)
                                             }
                                           }), 
          by = participant_id]
  
  
  dataset[!is.na(iuratio),  iuratio_traj := sapply(1:length(iuratio), 
                                                   function(i, 
                                                            absolute_inc = 20, 
                                                            interval_range = 2){
                                                     if(i < 2){
                                                       return(0)
                                                     }
                                                     
                                                     i_start <- find_index_start(i = i, 
                                                                                 interval_range = interval_range, 
                                                                                 follow_ups = time_between_obs, 
                                                                                 labs_val = iuratio,
                                                                                 prev_obs = F)
                                                     
                                                     if(is.na(i_start))
                                                       return(0)
                                                     
                                                     increase <- (iuratio[i] - iuratio[i_start:(i-1)])
                                                     
                                                     if(any(is.na(increase))){
                                                       increase <- increase[-which(is.na(increase))]
                                                     }
                                                     
                                                     if(any(increase >= absolute_inc)){
                                                       return(1)
                                                     }else{
                                                       return(0)
                                                     }
                                                   }), 
          by = participant_id]
  
  ###### Trajectory hgb (version for PANGEA 1.0), note different variable name traj_hgb
  dataset[!is.na(hgb), traj_hgb := {
    traj <- traj_f(- hgb, participant_id,
                   tstart,
                   treshold_pvalue = 0.5)
    as.numeric(traj > 0)
  }]
  dataset[is.na(hgb), traj_hgb := 0]
  # in the original paper traj_hgb was 0 even when hgb was missing (i.e. no evidence
  # of decreasing trajectory)
  
  
  
  return(dataset)
}

### Function that, given a dataset, return the dataset with various risk scores computed
computation_risks <- function(dataset, timepoint = 2,
                              models_folder = "../Models/",
                              models = c("pangea_bm_traj", "pangea_no_bm_traj", 
                                         "pangea_bm", "pangea_no_bm", "rolling_20_2_20")){
  
  dataset <- data.table(dataset)   # so changes will not be made in place
  
  # if the trajectory variables are not already in dataset, compute them
  if(!("mspike_traj" %in% names(dataset)) & any(grepl("traj", models))){
    dataset <- computation_trajectories(dataset)
  }
  
  if("pangea_bm_traj" %in% models){
    load(paste0(models_folder, "pangea_bm_traj.RData"), verbose = T)
    load(paste0(models_folder, "pangea_bm_traj_original_survfit.rda"), verbose = T)
    
    d_traj_bm <- dataset[pangea_bm_ok == T,]
    d_traj_bm$pangea_traj_bm_lin_pred <- predict(model_traj_bm, 
                                            d_traj_bm, 
                                            type = "lp")
    
    # manually construct a survfit for d_traj_bm_ok using the original one
    surv_traj_bm <- traj_bm_original_survfit
    n_traj_bm <- nrow(d_traj_bm)
    surv_traj_bm$surv0 <- replicate(n_traj_bm, 
                                           traj_bm_original_survfit$surv0)
    surv_traj_bm$surv <- t(t((surv_traj_bm$surv0))^exp(d_traj_bm$pangea_traj_bm_lin_pred))
    
    # use this surv and get_prob_pfs to calculate the risk score
    varname <- paste("pangea_bm_traj_", timepoint, "_year_risk", sep = "")
    dataset[pangea_bm_ok == T,
            (varname) := sapply(1:n_traj_bm, function(i) {
              get_risk_prob(surv_traj_bm$surv[,i], 
                            surv_traj_bm$time, 
                            year = timepoint)
              })]
  }
  
  if("pangea_no_bm_traj" %in% models){
    load(paste0(models_folder, "pangea_no_bm_traj.RData"), verbose = T)
    load(paste0(models_folder, "pangea_no_bm_traj_original_survfit.rda"), verbose = T)
    
    d_traj_no_bm <- dataset[pangea_no_bm_ok == T,]
    d_traj_no_bm$pangea_traj_no_bm_lin_pred <- predict(model_traj_no_bm, 
                                            d_traj_no_bm, 
                                            type = "lp")
    
    # manually construct a survfit for d_traj_bm_ok using the original one
    surv_traj_no_bm <- traj_no_bm_original_survfit
    n_traj_no_bm <- nrow(d_traj_no_bm)
    surv_traj_no_bm$surv0 <- replicate(n_traj_no_bm, 
                                       traj_no_bm_original_survfit$surv0)
    surv_traj_no_bm$surv <- t(t((surv_traj_no_bm$surv0))^exp(d_traj_no_bm$pangea_traj_no_bm_lin_pred))
    
    # use this surv and get_prob_pfs to calculate the risk score
    varname <- paste("pangea_no_bm_traj_", timepoint, "_year_risk", sep = "")
    dataset[pangea_no_bm_ok == T,
            (varname) := sapply(1:n_traj_no_bm, function(i) {
              get_risk_prob(surv_traj_no_bm$surv[,i], 
                            surv_traj_no_bm$time, 
                            year = timepoint)
            })]
    
  }
  
  if("pangea_bm_traj_alt" %in% models){
    load(paste0(models_folder, "pangea_bm_traj_alt.RData"), verbose = T)
    load(paste0(models_folder, "pangea_bm_traj_original_survfit_alt.rda"), verbose = T)
    
    d_traj_bm_alt <- dataset[pangea_bm_ok == T,]
    d_traj_bm_alt$pangea_traj_bm_alt_lin_pred <- predict(model_traj_bm_alt, 
                                                         d_traj_bm_alt, 
                                                         type = "lp")
    
    # manually construct a survfit for d_traj_bm_ok using the original one
    surv_traj_bm_alt <- traj_bm_original_survfit_alt
    n_traj_bm_alt <- nrow(d_traj_bm_alt)
    surv_traj_bm_alt$surv0 <- replicate(n_traj_bm_alt, 
                                    traj_bm_original_survfit_alt$surv0)
    surv_traj_bm_alt$surv <- t(t((surv_traj_bm_alt$surv0))^exp(d_traj_bm_alt$pangea_traj_bm_alt_lin_pred))
    
    # use this surv and get_prob_pfs to calculate the risk score
    varname <- paste("pangea_bm_traj_alt_", timepoint, "_year_risk", sep = "")
    dataset[pangea_bm_ok == T,
            (varname) := sapply(1:n_traj_bm_alt, function(i) {
              get_risk_prob(surv_traj_bm_alt$surv[,i], 
                            surv_traj_bm_alt$time, 
                            year = timepoint)
            })]
  }
  
  if("pangea_no_bm_traj_alt" %in% models){
    load(paste0(models_folder, "pangea_no_bm_traj_alt.RData"), verbose = T)
    load(paste0(models_folder, "pangea_no_bm_traj_original_survfit_alt.rda"), verbose = T)
    
    d_traj_no_bm_alt <- dataset[pangea_no_bm_ok == T,]
    d_traj_no_bm_alt$pangea_traj_no_bm_alt_lin_pred <- predict(model_traj_no_bm_alt, 
                                                               d_traj_no_bm_alt, 
                                                               type = "lp")
    
    # manually construct a survfit for d_traj_bm_ok using the original one
    surv_traj_no_bm_alt <- traj_no_bm_original_survfit_alt
    n_traj_no_bm_alt <- nrow(d_traj_no_bm_alt)
    surv_traj_no_bm_alt$surv0 <- replicate(n_traj_no_bm_alt, 
                                        traj_no_bm_original_survfit_alt$surv0)
    surv_traj_no_bm_alt$surv <- t(t((surv_traj_no_bm_alt$surv0))^exp(d_traj_no_bm_alt$pangea_traj_no_bm_alt_lin_pred))
    
    # use this surv and get_prob_pfs to calculate the risk score
    varname <- paste("pangea_no_bm_traj_alt_", timepoint, "_year_risk", sep = "")
    dataset[pangea_no_bm_ok == T,
            (varname) := sapply(1:n_traj_no_bm_alt, function(i) {
              get_risk_prob(surv_traj_no_bm_alt$surv[,i], 
                            surv_traj_no_bm_alt$time, 
                            year = timepoint)
            })]
  }
  
  if("pangea_bm" %in% models){
    load(paste0(models_folder, "mod_dep_time.rda"), verbose = T)
    load(paste0(models_folder, "pangea_bm_original_survfit.rda"), verbose = T)
    
    d_bm <- dataset[pangea_bm_ok == T,]
    d_bm$pangea_bm_lin_pred <- predict(mod_dep_time, 
                                       d_bm, 
                                       type = "lp")
    
    # manually construct a survfit for d_traj_bm_ok using the original one
    surv_bm <- bm_original_survfit
    n_bm <- nrow(d_bm)
    surv_bm$surv0 <- replicate(n_bm, bm_original_survfit$surv0)
    surv_bm$surv <- t(t((surv_bm$surv0))^exp(d_bm$pangea_bm_lin_pred))
    
    # use this surv and get_prob_pfs to calculate the risk score
    varname <- paste("pangea_bm_", timepoint, "_year_risk", sep = "")
    dataset[pangea_bm_ok == T,
            (varname) := sapply(1:n_bm, function(i) {
              get_risk_prob(surv_bm$surv[,i], 
                            surv_bm$time, 
                            year = timepoint)
            })]
    
  }
  
  if("pangea_no_bm" %in% models){
    load(paste0(models_folder, "mod_dep_time_no_BM.rda"), verbose = T)
    load(paste0(models_folder, "pangea_no_bm_original_survfit.rda"), verbose = T)
    
    d_no_bm <- dataset[pangea_no_bm_ok == T,]
    d_no_bm$pangea_no_bm_lin_pred <- predict(mod_dep_time_no_BM, 
                                             d_no_bm, 
                                             type = "lp")
    
    # manually construct a survfit for d_traj_bm_ok using the original one
    surv_no_bm <- no_bm_original_survfit
    n_no_bm <- nrow(d_no_bm)
    surv_no_bm$surv0 <- replicate(n_no_bm, no_bm_original_survfit$surv0)
    surv_no_bm$surv <- t(t((surv_no_bm$surv0))^exp(d_no_bm$pangea_no_bm_lin_pred))
    
    # use this surv and get_prob_pfs to calculate the risk score
    varname <- paste("pangea_no_bm_", timepoint, "_year_risk", sep = "")
    dataset[pangea_no_bm_ok == T,
            (varname) := sapply(1:n_no_bm, function(i) {
              get_risk_prob(surv_no_bm$surv[,i], 
                            surv_no_bm$time, 
                            year = timepoint)
            })]
    
  }
  
  if("rolling_20_2_20" %in% models){
    if (!("risk_rolling_20_2_20_prob_2_years" %in% names(dataset))) {
      dataset[, mspike_greater_2 := as.numeric(mspike > 2)]
      dataset[, iuratio_greater_20 := as.numeric(iuratio > 20)]
      dataset[, plasmacells_grater_20 := as.numeric(plasmacells > 20)]
      
      ## If at least 2 are NAs, we want the result to be NA
      dataset[, is_20_2_20_NA := apply(cbind(mspike_greater_2, 
                                             iuratio_greater_20, 
                                             plasmacells_grater_20),
                                       1, function(x) sum(is.na(x)) )]
      
      dataset[is_20_2_20_NA != 2, ind_20_2_20 := apply(cbind(mspike_greater_2, 
                                                             iuratio_greater_20, 
                                                             plasmacells_grater_20),
                                                       1, sum, na.rm = T)]
      
      
      dataset[, 
              rolling_20_2_20 := ifelse(ind_20_2_20 == 0 , "Low",
                                        ifelse(ind_20_2_20 == 1, "Intermediate", "High"))]
      
      dataset[, ':=' (mspike_greater_2 = NULL, iuratio_greater_20 = NULL, plasmacells_grater_20 = NULL,
                      is_20_2_20_NA = NULL,ind_20_2_20 = NULL )]
      
      dataset[!is.na(rolling_20_2_20), 
              risk_rolling_20_2_20_prob_2_years:= ifelse(rolling_20_2_20 == "High", 
                                                         0.442,
                                                         
                                                         ifelse(rolling_20_2_20 == "Intermediate", 0.179, 0.062))]
      
    }
    
    # also get baseline 20_2_20
    dataset[, baseline_20_2_20 := rolling_20_2_20[1], 
            by = participant_id]
    dataset[, risk_baseline_20_2_20_prob_2_years := risk_rolling_20_2_20_prob_2_years[1], 
            by = participant_id]
       # alternatively, we could base these on the first lab where 20/2/20 is not NA
  }
  
  return(dataset)
}


### function to evaluate model calibration for a given cohort
# compares the pangea model's expected # of events by the timepoint to the observed #
# timepoint is used to get the actual event rate; this should match the timepoint used
# in the risk prediction (given by pred_varname)
cohort_calib <- function(cohort_data, 
                         pred_varname = "pangea_bm_traj_2_year_risk",
                         timepoint = 2) {
  if(! pred_varname %in% names(cohort_data)) {
    stop(paste("The variable", pred_varname, "was not found in the supplied data set."))
  }
  
  cohort_data <- data.table(cohort_data)
  
  # predicted event rate (pre-existing in the dataset)
  pred_event_rate <- cohort_data[, mean(get(pred_varname), na.rm = T)]
  pred_event_rate_lower <- cohort_data[, quantile(get(pred_varname), 
                                                    probs = (0.025), 
                                                    na.rm = T)]
  pred_event_rate_upper <- cohort_data[, quantile(get(pred_varname), 
                                                    probs = (0.975), 
                                                    na.rm = T)]
  
  # observed event rate (based on simple KM curve)
  cohort_data[, time_to_mm_or_cens := last_date - tstart]
  cohort_surv <- survfit(Surv(time = time_to_mm_or_cens, 
                              event = ever_mm) 
                         ~ 1, 
                         data = cohort_data)
  
  time_index <- which.min(abs(cohort_surv$time - timepoint))
  actual_event_rate <- 1 - cohort_surv$surv[time_index]
  actual_event_rate_upper <- 1 - cohort_surv$lower[time_index]   # flip label b/c negative
  actual_event_rate_lower <- 1 - cohort_surv$upper[time_index]   # flip label b/c negative
  
  # output
  c(ratio = actual_event_rate/pred_event_rate,
    actual_event_rate = actual_event_rate, 
    pred_event_rate = pred_event_rate,
    pred_event_rate_lower = as.vector(pred_event_rate_lower),
    pred_event_rate_upper = as.vector(pred_event_rate_upper),
    actual_event_rate_lower = actual_event_rate_lower,
    actual_event_rate_upper = actual_event_rate_upper,
    n_patients = length(unique(cohort_data$participant_id)),
    n_observations = nrow(cohort_data))
}

### function to compute Kaplan Meier PFS curves for a cohort, stratified
### by baseline 20/2/20
risks_profile <- function(dataset, risk_cat_varname = "rolling_20_2_20", return_surv = F,
                          whole_cohort = F,
                          risk_categories = c("Low", "Intermediate", "High"),
                          model_name = "20/2/20"){
  #### Building the dataset
  ## For each participant, we create an additional variable for the first ob
  ## at which we can compute 20_2_20
  
  dataset <- data.table(dataset)
  
  dataset[, ok_risk_cat := !is.na(get(risk_cat_varname))]
  dataset[, first_risk_cat := cumsum(ok_risk_cat) == 1,
          by = participant_id]
  data_risks_profiles <- dataset[first_risk_cat == T, ]
  ## We update the tstop to be the final obs - the new tstart,
  # tstart = 0 
  # and prog_mm to be equal to ever_mm
  data_risks_profiles[, tstop := last_date - tstart]
  data_risks_profiles[, tstart := 0]
  data_risks_profiles[, prog_mm := ever_mm]
  
  ## Now, we create the survfit object to see the actual progression rates
  print(data_risks_profiles[, table(get(risk_cat_varname))])
  empirical_risks_classification <- survfit(as.formula(paste("Surv(tstart, tstop, prog_mm) ~ ",
                                              risk_cat_varname)), 
                                            data = data_risks_profiles)
  print(empirical_risks_classification)
  surv_2_year_summary <- summary(empirical_risks_classification, time = 2)
  surv_median <- quantile(empirical_risks_classification, probs = 0.5)
  
  # construct title
  plot_title <- if (whole_cohort == F) {
    paste("Progression free survival by ", model_name, "\nrisk category",
          sep = "")
  } else if (whole_cohort == T) {
    "Progression free survival"
  }
  
  legend_labs <- if (whole_cohort == F) {
    full_labels <- c("Low", "Intermediate", "High")
    actual_labels <- str_remove(names(empirical_risks_classification$strata), ".*=")
    full_labels[full_labels %in% actual_labels]
  } else if (whole_cohort == T) {
    "Overall"
  }
  
  ## Graphical representation
  p <- ggsurvplot(empirical_risks_classification, 
                  data = data_risks_profiles, 
                  conf.int = TRUE,
                  xlim = c(0, 10),
                  legend.labs = legend_labs)  %++% 
    scale_x_continuous(breaks = c(0, 2, 5, 10),
                       limits = c(0, 10)) %++% 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1)) %++%
    geom_vline(xintercept = 2,
               linetype = 2) %++%
    labs(title = plot_title,
         x = "Time (years since baseline)")
  
  
  # return information as appropriate
  list(plot = p, 
       surv_summary = surv_2_year_summary,
       surv_median = surv_median,
       surv_object = if(return_surv) empirical_risks_classification else NULL)
}


### function to plot the distribution of the difference between the 2.0 and 1.0 risk scores
risk_diff_plot <- function(dataset, 
                           pangea_bm_traj_varname = "pangea_bm_traj_2_year_risk",
                           pangea_no_bm_traj_varname = "pangea_no_bm_traj_2_year_risk") {
  dataset <- data.table(dataset)
  
  dataset[, any_traj := ifelse(mspike_traj + creatinine_traj + iuratio_traj + hgb_traj > 0,
                               "Yes",
                               "No")]
  
  risk_diff_dat <- rbind(dataset[, list(dif = get(pangea_bm_traj_varname) - pangea_bm_2_year_risk,
                                        traj = any_traj,
                                        model = "BM Model")],
                         dataset[, list(dif = get(pangea_no_bm_traj_varname) - pangea_no_bm_2_year_risk,
                                        traj = any_traj,
                                        model = "No BM Model")])
  
  ggplot(data = risk_diff_dat, 
         aes(x = dif,
             fill = traj, color = traj)) + 
    geom_density(adjust = 3, alpha = 0.5) + 
    theme_minimal() + 
    facet_wrap(facets = vars(model)) +
    xlim(c(-0.5, 0.5)) + 
    labs(fill = "Any trajectory?",
         color = "Any trajectory?",
         x = "PANGEA 2-year risk 2.0 risk score minus PANGEA 2-year risk 1.0 risk score")
}


### function to compute compute c stat for a given dataset, model, and timepoint
c_stat_time <- function(dataset, mod, given_time, 
                        mod_type = c("pangea_bm", "pangea_no_bm", "20_2_20")) 
{
  # grab mod_type value
  mod_type <- match.arg(mod_type)
  
  # get corresponding variable name to check if a visit has a complete set of labs
  complete_visit_varname <- ifelse(mod_type == "pangea_no_bm", "pangea_no_bm_ok", "pangea_bm_ok")
  
  #Collect the dataset based on participant id 
  #such that the last *complete* lab values chosen before given_time
  df_subset <- dataset[last_date >= given_time]
  df_subset[, last_complete_visit := which.max(tstart[(tstart <= given_time) &
                                                        (get(complete_visit_varname) == TRUE)]),
            by = participant_id]
  dt <- df_subset[tstart <= given_time & 
                    !is.na(last_complete_visit), 
                  .SD[last_complete_visit[1]], by = participant_id]
  
  if (nrow(dt) <= 1) {
    return(rep(NA, 5))
  }
  
  # update tstop and prog_mm to use the *full* followup for the patient, not just the current interval
  # this is necessary because the data are in the counting process format, but now we're doing an
  # analysis of just a fixed timepoint. the concordance function will look for these variables
  # since they're what are in the model.
  dt[, tstop := last_date]
  dt[, prog_mm := ever_mm]
  
  concord <- concordance(mod, newdata = dt)
  c_stat <- concord$concordance
  
  if(!is.na(c_stat) ){
    c_se <- sqrt(concord$var)
    c_lb <- c_stat - qnorm(0.975, mean = 0, sd = 1)*c_se
    c_ub <- c_stat + qnorm(0.975, mean = 0, sd = 1)*c_se
    return(c(c_stat = c_stat, 
             ci_lower = c_lb, 
             ci_upper = c_ub, 
             n = nrow(dt), 
             n_prog_mm = sum(dt$prog_mm)))
  }
  else
  {
    return(c(c_stat = c_stat, 
             ci_lower = NA, 
             ci_upper = NA,  
             n = nrow(dt), 
             n_prog_mm = sum(dt$prog_mm)))
  }
}


### function to compute c-statistics over a sequence of timepoints
c_stat_time_matrix <- function(dataset, mod, time_min, time_max, 
                               time_increment = 0.1,
                               mod_type = c("pangea_bm", "pangea_no_bm", "20_2_20")) {
  # grab mod_type value
  mod_type <- match.arg(mod_type)
  
  # set up vector of times 
  c_times <- seq(from = time_min, to = time_max, by = time_increment)
  
  # compute c stats
  c_matrix <- t(sapply(c_times, function(x) c_stat_time(dataset = dataset,
                                                        mod = mod,
                                                        given_time = x,
                                                        mod_type = mod_type)))
  c_matrix <- data.table(cbind(c_times, c_matrix))
  return(c_matrix)
}


#gives us the plot based on the given dataset, trained model, and time_rage.
#we also compare with 20/2/20
c_stat_plot <- function(dataset, mod_pangea,
                      time_min, time_max, time_increment = 1,
                      mod_20 = mod_20_2_20,
                      mod_type = c("pangea_bm", "pangea_no_bm", "20_2_20")){
  # grab mod_type value
  mod_type <- match.arg(mod_type)
  
  # Generate datasets
  c_matrix <- c_stat_time_matrix(dataset = dataset, mod = mod_pangea, 
                                 time_min = time_min, time_max = time_max, 
                                 time_increment = time_increment, 
                                 mod_type = mod_type)
  c_matrix_20 <- c_stat_time_matrix(dataset = dataset, mod = mod_20, 
                                    time_min = time_min, time_max = time_max, 
                                    time_increment = time_increment, 
                                    mod_type = mod_type)
  
  # Dynamically set model type based on presence of the "plasmacells" coefficient
  if(is.na(mod_pangea$coefficients["plasmacells"])){
    model_type_label <- "Model NO BM"
  } else {
    model_type_label <- "Model BM"
  }
  
  # Add model labels
  add_col_model <- rep(model_type_label, nrow(c_matrix))
  add_col_20_2_20 <- rep("Model 20/2/20", nrow(c_matrix_20))
  
  # Combine datasets
  c_matrix1 <- cbind(Model = add_col_model, c_matrix)
  c_matrix2 <- cbind(Model = add_col_20_2_20, c_matrix_20)
  cm <- rbind(c_matrix1, c_matrix2)
  
  # Define the dodge width for spreading out the points and error bars
  dodge <- position_dodge(width = 0.3)
  # Plot
  ggplot(cm, aes(x = c_times, y = c_stat, color = Model, alpha = n)) +
    geom_point(mapping = aes(size = n)) +
    geom_line(aes(color = Model)) +  # Add points to the lines
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, 
                  position = dodge,
                  linetype = "dashed") + 
    
    # Text for # of patients in blue and bold
    #geom_text(data = subset(cm, Model == model_type_label),
    geom_text(data = cm,
              aes(label = n), 
              vjust = 1.5,
              hjust = 1.2,
              show.legend = FALSE,
              color = "blue",
              fontface = "bold",
              position = dodge) +
    
    # Text for # of progressors in purple and bold
    geom_text(data = cm,
              aes(label = n_prog_mm), 
              vjust = -0.8, 
              hjust = -0.5,
              show.legend = FALSE, 
              color = "purple",
              fontface = "bold",
              position = dodge) +
    
    scale_size_continuous(range = c(1, 4)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    
    # Use setNames to create a dynamic color mapping
    scale_color_manual(values = 
                         setNames(c("darkolivegreen", "red3", "blue", "purple"),
                                  c(model_type_label, "Model 20/2/20", "# of patients", "# of progressors")),
                       name = " ") +
    
    scale_y_continuous(limits = c(0.3, 1)) +
    
    theme_minimal() +
    
    labs(x = "Years since diagnosis",
         y = "C-statistic",
         size = "# of patients",
         alpha = "# of patients",
         color = "Model Type") + 
    
    # Add fake points for manual legend
    geom_point(aes(x = Inf, y = Inf, color = "# of patients"), 
               show.legend = TRUE, size = 0) +  # Blue for "# of patients"
    geom_point(aes(x = Inf, y = Inf, color = "# of progressors"), 
               show.legend = TRUE, size = 0) +  # Purple for "# of progressors"
    
    # Customize the legend
    guides(color = guide_legend(override.aes = list(size = 5)),
           size = guide_legend(title = "# of patients"),
           alpha = guide_legend(title = "# of patients")) +
    
    theme(legend.text = element_text(size = 12),  # Adjust the size as needed
          legend.title = element_text(size = 14),  # Adjust the size as needed
          axis.title.x = element_text(size = 14),   # Increase x-axis label size
          axis.title.y = element_text(size = 14),   # Increase y-axis label size
          axis.text = element_text(size = 12))      # Increase axis tick label size 
}


### function to compute calibration statistics for a range of risk predictions
### and subcohorts
calibration_stats <- function(dataset, 
                              pred_varnames = c("pangea_bm_traj_2_year_risk",
                                                "pangea_no_bm_traj_2_year_risk",
                                                "pangea_bm_2_year_risk",
                                                "pangea_no_bm_2_year_risk",
                                                "risk_rolling_20_2_20_prob_2_years"),
                              conditions = c("TRUE",     # this gives the overall cohort (TRUE for all patients)
                                             "rolling_20_2_20 == \"Low\"",
                                             "rolling_20_2_20 == \"Intermediate\"",
                                             "rolling_20_2_20 == \"High\"",
                                             "mspike < 2",
                                             "mspike >= 2",
                                             "iuratio < 20",
                                             "iuratio >= 20",
                                             "plasmacells < 20",
                                             "plasmacells >= 20",
                                             "creatinine < 1",
                                             "creatinine >= 1",
                                             "age < 65",
                                             "age >= 65")) {
  # make pred_varnames a named vector so we get a named list from the lapply
  names(pred_varnames) <- pred_varnames
  
  # compute the calibration statistic for each combination of risk prediction and
  # cohort condition
  lapply(pred_varnames, function(pred_varname) {
    sapply(conditions, function(condition) {
      cohort <- dataset[!is.na(get(pred_varname)) & eval(str2expression(condition))]
      if (nrow(cohort) == 0) {
        rep(NA, 9)
      } else {
        cohort_calib(cohort,
                     pred_varname = pred_varname)
      }
    })
  })
  
}

### function to format for reporting the sensitivity/specificity/ppv/npv output 
### from the SeSpPPVNPV function in the timeROC package
process_sspn <- function(sspn) {
  sspn_pt <- sapply(sspn[c("TP", "FP", "PPV", "NPV")], function(x) x)[2,]
  sspn_pt[2] <- 1 - sspn_pt[2]            # convert FP to specificity
  names(sspn_pt) <- c("Sens", "Spec", "PPV", "NPV")
  
  sspn_ses <- c(sspn$inference$vect_se_Se[2], 
                sspn$inference$vect_se_Sp1[2],
                sspn$inference$vect_se_PPV[2],
                sspn$inference$vect_se_NPV1[2])
  names(sspn_ses) <- names(sspn_pt)
  
  sspn_lower <- sspn_pt + qnorm(0.025)*sspn_ses
  sspn_upper <- sspn_pt + qnorm(0.975)*sspn_ses
  
  sspn_report <- paste(round(sspn_pt, 3)*100, 
                       " (", round(sspn_lower, 3)*100,
                       " - ", round(sspn_upper, 3)*100,
                       ")", sep = "")
  names(sspn_report) <- names(sspn_pt)
  
  sspn_report
}

### function to compute the time-dependent sensitivity/specificity/ppv/npv of a 
### marker in a given dataset
time_dep_sspn <- function(dat, 
                          marker_var, 
                          high_cut = 0.4, 
                          eval_times = c(0, 0.1, 1, 2, 3, 4, 5),
                          obs_ok_var = "pangea_bm_ok") {
  data.table(t(sapply(eval_times, function(x) {
    dat[get(obs_ok_var) == T & 
          tstart <= x & tstop >= x, 
        {
          sspn <- SeSpPPVNPV(cutpoint = high_cut, 
                             T = last_date, 
                             delta = ever_mm,
                             marker = get(marker_var),
                             cause = 1, 
                             times = x + 2,
                             iid = T)
          sspn_report <- process_sspn(sspn)
          
          c(s = x,
            t = x + 2,
            sspn_report,
            n_risk = .N,
            n_prog = sum(ever_mm == 1 & last_date < x + 2),
            n_pred_high = sspn$Stats[2,"Positive (X>c)"],
            n_pred_low = sspn$Stats[2,"Negative (X<=c)"])
        }]
  })))
}


