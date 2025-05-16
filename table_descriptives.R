
### Create the column names (useful to cobine the results and keep track of what we do)
cohorts <- c("all", unique(data_validation$pi_cohort))
summaries_numeric_variables <- c("1_qrt", "med", "3_qrt", "NAs")
col_names_num <- do.call(paste, c(expand.grid(summaries_numeric_variables, cohorts)[,c(2,1)], sep = "_"))

summaries_cat_variables <- c("n", "prop")
col_names_cat <- do.call(paste, c(expand.grid(summaries_cat_variables, cohorts)[,c(2,1)], sep = "_"))

## We replicate table 1
data_validation[, first_obs := tstart == min(tstart), by = participant_id]
data_validation[, last_obs := tstart == max(tstart), by = participant_id]


# Number of patients
n_patients_all <- data_validation[first_obs == T, .N]


## Info age, number of visits, median follow up
quart_probs <- c(0.25, 0.5, 0.75)

# Age at diagnosis
age_info <- data_validation[first_obs == T, quantile(age, quart_probs, na.rm = T)]
age_info <- c(age_info, 
              NAs = data_validation[first_obs == T, sum(is.na(age))])

# Number of visits
data_validation[, n_visits:= .N, by = participant_id]

nvisits_info <- data_validation[first_obs == T, quantile(n_visits, quart_probs, na.rm = T)]
nvisits_info <- c(nvisits_info, 
                  NAs = data_validation[first_obs == T, sum(is.na(n_visits))])

# Length follow up (in years)
follow_up_info <- data_validation[last_obs == T, quantile(tstop, quart_probs, na.rm = T)]
follow_up_info <- c(follow_up_info, 
                    NAs = data_validation[last_obs == T, sum(is.na(tstop))])

# Length follow up (in years) using reverse Kaplan Meier method
follow_up_info <- data_validation[last_obs == T, quantile(tstop, quart_probs, na.rm = T)]
follow_up_info <- c(follow_up_info, 
                    NAs = data_validation[last_obs == T, sum(is.na(tstop))])

follow_up_reverse_km <- survfit(Surv(tstop, 1 - prog_mm) ~ 1, 
                                data = data_validation[last_obs == T])
follow_up_rkm_info <- quantile(follow_up_reverse_km, probs = c(0.25, 0.5, 0.75))$quantile
follow_up_rkm_info <- c(follow_up_rkm_info, NAs = 0)

# Interval between visits (in month)
data_validation[, time_between_obs := (tstop - tstart)*12]

time_between_obs_info <- data_validation[, quantile(time_between_obs, quart_probs, na.rm = T)]
time_between_obs_info <- c(time_between_obs_info, 
                           NAs = data_validation[, sum(is.na(time_between_obs))])


### organize the summaries of the numeric variables into a single table 
table_numeric_info <- rbind(age_info,
                            nvisits_info,
                            follow_up_info,
                            follow_up_rkm_info,
                            time_between_obs_info)


### organize the summaries of the categorical variables into a list of tables
### Baseline - Categorical anagraphic variables

# Sex
table_sex <- data_validation[first_obs == T,
                             table(sex, useNA = "always")]
table_sex_prop <- round(prop.table(table_sex), 2)

# Race
table_race <- data_validation[first_obs == T, 
                              table(race, useNA = "always")]
table_race_prop <- round(prop.table(table_race), 2)

# Race
table_ethnicity <- data_validation[first_obs == T, 
                              table(ethnicity, useNA = "always")]
table_ethnicity_prop <- round(prop.table(table_ethnicity), 2)

# Progression
table_progression <- data_validation[last_obs == T,
                                     table(prog_mm, useNA = "always")]
table_progression_prop <- round(prop.table(table_progression), 2)

# Immunofixation
table_immunofix <- data_validation[first_obs == T, 
                                   table(immunofix2, useNA = "ifany")]
table_immunofix_prop <- round(prop.table(table_immunofix), 2)

# Died or censored for treatment
table_end_date <- data_validation[last_obs == T, 
                                  table(end_date_type, useNA = "ifany")]
table_prop_end_date <- round(prop.table(table_end_date), 2)

## organize them into a list
list_cat_info <- list(counts = list(sex = table_sex, race = table_race,
                                    ethnicity = table_ethnicity,
                                    progression = table_progression, 
                                    immunofix = table_immunofix, 
                                    end_date = table_end_date),
                      rates = list(sex = table_sex_prop, race = table_race_prop,
                                   ethnicity = table_ethnicity_prop,
                                   progression = table_progression_prop, 
                                   immunofix = table_immunofix_prop,
                                   end_date = table_prop_end_date))


