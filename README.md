# Paper information

<b>Enhanced Dynamic Risk Stratification of Smoldering Multiple Myeloma</b>

Published in <i>Nature Medicine</i>: <a href="https://doi.org/10.1038/s41591-026-04304-x">10.1038/s41591-026-04304-x</a>

## Authors

Floris Chabrun*<sup>1,2,3</sup>, Daniel E Schwartz*<sup>2,4</sup>, Susanna Gentile*<sup>5</sup>, Elias K Mai*<sup>6</sup>, Tulika R Gupta<sup>1,7</sup>, Jacqueline Perry<sup>1</sup>, David M Cordas dos Santos<sup>1,2</sup>, Thomas Hielscher<sup>8</sup>, Annika Werly<sup>6</sup>, Sophia K Schmidt<sup>6</sup>, Foteini Theodorakakou<sup>9</sup>, Despina Fotiou<sup>9</sup>, Christine Ivy Liacos<sup>9</sup>, Nikolaos Kanellias<sup>9</sup>, Noelia Collado Gisbert<sup>10</sup>, Esperanza Martin-Sanchez<sup>10</sup>, Rosalinda Termini<sup>10</sup>, Johannes Waldschmidt<sup>11</sup>, Selina J Chavda<sup>12</sup>, Louise Ainley<sup>12</sup>, Matteo Claudio Da Vià<sup>13</sup>, Claudio De Magistris<sup>13</sup>, Loredana Pettine<sup>13</sup>, Michael A Timonian<sup>1,2</sup>, Jean-Baptiste Alberge<sup>1,2,14</sup>, Vidhi Patel<sup>1</sup>, Patrick Costello<sup>1</sup>, Catherine Tobia<sup>1</sup>, Sally Phan<sup>1</sup>, Jennifer Lamb<sup>1</sup>, Maria-Theresa Silverio<sup>1</sup>, Maya Davis<sup>1</sup>, Elizabeth K O’Donnell<sup>1,2</sup>, Catherine R Marinac<sup>1,2</sup>, Omar Nadeem<sup>1,2</sup>, Niccolo Bolli<sup>13,15</sup>, Kwee Yong<sup>12</sup>, K. Martin Kortüm<sup>11</sup>, Hermann Einsele<sup>11</sup>, María-Victoria Mateos<sup>16</sup>, Shaji Kumar<sup>17</sup>, Jesus San Miguel<sup>10</sup>, Bruno Paiva<sup>10</sup>, Efstathios Kastritis<sup>9</sup>, Meletios A Dimopoulos<sup>9,18</sup>, Marc S Raab<sup>6</sup>, Lorenzo Trippa<sup>+</sup><sup>19</sup>, Irene M Ghobrial<sup>+</sup><sup>1,2</sup>

*these authors contributed equally
<br><sup>+</sup>these authors jointly supervised this work

## Affiliations

<br><sup>1</sup> Center for Early Detection and Interception-Blood Cancers, Department of Medical Oncology, Dana-Farber Cancer Institute, Boston, MA
<br><sup>2</sup> Harvard Medical School, Boston, MA
<br><sup>3</sup> Laboratory of Biochemistry and Molecular Biology, University Hospital of Angers, University of Angers, Angers, France 
<br><sup>4</sup> Massachusetts General Hospital Biostatistics, Somerville, MA
<br><sup>5</sup> Sapienza University of Rome, Rome, Italy
<br><sup>6</sup> Heidelberg Myeloma Center, Internal Medicine V, Hematology, Oncology and Rheumatology, Heidelberg University Hospital and Medical Faculty Heidelberg, Heidelberg, Germany
<br><sup>7</sup> Yale School of Medicine, New Haven, CT
<br><sup>8</sup> German Cancer Research Center (DKFZ), Division of Biostatistics, Heidelberg, Germany
<br><sup>9</sup> National and Kapodistrian University of Athens, Athens, Greece
<br><sup>10</sup> Cancer Center Clinica Universidad de Navarra (CCUN), CIMA Universidad de Navarra, Instituto de Investigación Sanitaria de Navarra (IDISNA), CIBER-ONC number CB16/12/00369
<br><sup>11</sup> Department of Internal Medicine II, University Hospital of Würzburg, Würzburg, Germany
<br><sup>12</sup> Department of Haematology, Cancer Institute, University College London, UK
<br><sup>13</sup> Hematology section, Fondazione IRCCS Ca’ Granda Ospedale Maggiore Policlinico, Milan, Italy
<br><sup>14</sup> Cancer Program, Broad Institute of MIT and Harvard, Cambridge, MA, USA
<br><sup>15</sup> Department of Oncology and Hemato-Oncology, University of Milan, Milan, Italy
<br><sup>16</sup> Department of Hematology, University Hospital of Salamanca/IBSAL/CIC/CIBERONC, Salamanca, Spain
<br><sup>17</sup> Department of Hematology, Mayo Clinic Comprehensive Cancer Center, Rochester, Minnesota, USA
<br><sup>18</sup> Department of Medicine, Korea University, Seoul, South Korea
<br><sup>19</sup> Department of Data Science, Dana-Farber Cancer Institute, Boston, MA

## Abstract

Accurate prediction of risk of progression from smoldering multiple myeloma (SMM) to active multiple myeloma (MM) is paramount to individualized early therapeutic strategies with minimum risk of overtreatment. Current risk stratification models do not account for evolving biomarker trajectories. We assembled a cohort of 2,344 patients with SMM from seven international centers with longitudinal clinical and biological data to train and validate the Precursor Asymptomatic Neoplasms by Group Effort Analysis (PANGEA)-SMM risk models. Four evolving biomarkers were significantly associated with shorter time to progression: M-protein increase ≥0.2 g dl−1, involved/uninvolved serum free light chain ratio increase ≥20, creatinine increase >25% and hemoglobin decrease ≥1.5 g dl−1. PANGEA-SMM outperforms established models, including the 20/2/20 and IMWG models, by more accurately predicting progression (C-statistic = 0.79), even without biomarker history (C-statistic = 0.78) or recent bone marrow biopsy (C-statistic = 0.78). We present PANGEA-SMM to the community as an easy-to-use, open-access tool for risk stratification in SMM. Validation tools are available to compare PANGEA-SMM to established models.

# Code overview

This folder contains the necessary code to compute the various validation results for the updated PANGEA 2.0 trajectory models (BM and no BM) and alternatives (original PANGEA and 20/2/20).

The folder includes:
- `validation_code.R`: main script to run to perform the analysis
- `data_preparation.R`: script computing the trajectories, risk scores, and some other variables needed for the analysis
- `table_descriptives.R`: script to compute the basic descriptive statistics of the cohort (sourced in validation\_code.R)
- `risk_reclassification.R`: script to compute risk strata from the PANGEA risk scores (sourced in validation\_code.R)
- `Functions.R`: script defining the functions needed to run the various analyses
- `Packages.R`: script loading the needed packages
- `cutoffs_pangea_training.rda`: file storing cutoffs used to determine risk strata
- `models`: folder storing the PANGEA models and necessary information to compute the risk predictions
- `results`: folder to save the results at the end of the analysis (currently empty)

# How to perform the analysis

1. Set the `Validation - Code for External Collaborators` folder as working directory
2. Open the R script `validation_code.R`
3. Modify the path to load your data (at the beginning of the script)
4. If applicable, on line 7 set `save_deidentified` to TRUE (instead of FALSE, the default value)
5. Run the `validation_code.R` script

The results will be saved in the results folder.

The `validation_code.R` executes the following steps:
1. Source the `Functions.R` and `Packages.R` scripts
2. Load the data
3. Source the `data_preparation.R` script to compute some variables needed for the analysis, including the trajectories, the rolling 20/2/20 classification, and PANGEA risk scores (including for some previous developmental versions of the model)
    - **Note**: To just compute the PANGEA 2.0 risk scores, simply run `validation_code.R` up to line 32 where `data_preparation.R` is sourced.
    - The PANGEA 2.0 models are saved in the `pangea_bm_traj_alt.RData` and `pangea_no_bm_traj_alt.RData` files in the models folder, and their risk score variables created by `data_preparation.R` have similar names.
4. Compute descriptive statistics for the cohort
5. Compute basic statistics on progression rates for the cohort
6. Compute descriptive statistics on the risk scores for the cohort
7. Compute concordance statistics (ability of the risk scores to appropriately rank patients)
8. Compute classification statistics (based on risk strata)
9. Compute calibration statistics (ability of risk scores to match the actual progression risk for various subgroups)
10. Store the results in the results folder


# Format of the data

The dataset has to include the following variables:
- `tstart`: Time at which the current labs are measured (it should be 0 for the first observation) in years; used for validation.
- `tstop`: Time of the next lab or the end of the follow up (final observation) in years; used for validation.
- `prog_mm`: Binary progression status at the time tstop (1 = slim-crab myeloma, 0 otherwise); used for validation.
- `mspike`: M-spike in g/dL; used for the PANGEA risk prediction.
- `iuratio`: Involved/Uninvolved free light chain ratio; used for the PANGEA risk prediction.
- `plasmacells`: Numeric percentage of plasma cells (0-100); used for the PANGEA risk prediction.
- `creatinine`: Creatinine in mg/dL; used for the PANGEA risk prediction.
- `age`: Age of the patient at each visit in years; used for the PANGEA risk prediction.
- `hgb`: Hemoglobin in g/dL (needed to compute the trajectory); used for the PANGEA risk prediction.
- `sex`: Sex of the patient; used only for descriptive statistics.
- `current_diagnosis`: Categorical variable, is the patient smoldering, i.e. "SMM" (for this analysis we want only smoldering patients); used only for descriptive statistics.
- `end_date_type`: Description of the event that caused censoring or progression. possible values for censored patients are "Death", "Last appointment", and "Other treatment". For patients who progress the value should just be "Diagnosis". Used only for descriptive statistics.
- `race`: Race of patient. Possible values are "Asian", "Black or African American", "Declined", "Multiple", "Other", and "White". For patients where race is missing (it's ok if that's everyone in your cohort), the value should be NA. Used only for descriptive statistics.
- `ethnicity`: Ethnicity of the patient. Possible values are "Declined", "Hispanic or Latino", and "Not Hispanic or Latino". For patients where ethnicity is missing (it's ok if that's everyone in your cohort), the value should be NA. Used only for descriptive statistics.
- `immunofix2`: Immunofixation isotype. Possible values are "Biclonal", "IgA", "IgG", and "Light Chain Only". Missing values should be indicated as NA. Used only for descriptive statistics.

The data has to be in a long format, which is, each patient visit corresponds to
one row in the dataset. 

**Note**: PANGEA risk predictions can only be computed for the rows in which all predictor variables (`mspike`, `iuratio`, `creatinine`, `hgb`, and `age`) are available (not NA). Because of this, labs associated with the same patient visit but reported on slightly different dates (e.g. mspike on 5/9/24 and iuratio on 5/4/24) should be combined into a single row with a single date in order to get a risk prediction for that visit.

For a visual aide, please see the first rows of our dataset, 
obtained using the R function `tmerge` from the package `survival`. 
The data are loaded in `data.table` format.


``` r
## Dataset
library(data.table)
sample_data <- fread("fake_example_data.csv")
sample_data[participant_id == 1,]
```

```
##       V1 participant_id    tstart     tstop prog_mm mspike   iuratio plasmacells      age creatinine   hgb
##    <int>          <int>     <num>     <num>   <int>  <num>     <num>       <int>    <num>      <num> <num>
## 1:     1              1 0.0000000 0.9117243       0   0.37  7.903827          NA 42.96655       0.44  12.0
## 2:     2              1 0.9117243 1.4100240       0   1.93 17.548807          NA 43.87828       1.15  10.2
## 3:     3              1 1.4100240 2.7132696       0   1.58 19.952914          NA 44.37658       0.74  10.2
## 4:     4              1 2.7132696 3.2745413       0   1.39 67.277014          NA 45.67982       0.34  10.4
##    current_diagnosis    end_date_type                      race          ethnicity immunofix2    sex
##               <char>           <char>                    <char>             <char>     <char> <char>
## 1:               SMM Last appointment Black or African American Hispanic or Latino        IgA   Male
## 2:               SMM Last appointment Black or African American Hispanic or Latino        IgA   Male
## 3:               SMM Last appointment Black or African American Hispanic or Latino        IgA   Male
## 4:               SMM Last appointment Black or African American Hispanic or Latino        IgA   Male
```

``` r
sample_data[participant_id == 2,]
```

```
##       V1 participant_id    tstart     tstop prog_mm mspike  iuratio plasmacells      age creatinine   hgb
##    <int>          <int>     <num>     <num>   <int>  <num>    <num>       <int>    <num>      <num> <num>
## 1:     5              2 0.0000000 0.2956944       0   2.00       NA          NA 77.00214       0.56   9.8
## 2:     6              2 0.2956944 0.5256789       0   2.50 59.73815          45 77.29784       0.17   8.4
## 3:     7              2 0.5256789 0.8076837       0   4.11 78.20539          45 77.52782       0.27  10.1
## 4:     8              2 0.8076837 1.3032455       0   3.69 95.81232          45 77.80983       0.72  10.0
## 5:     9              2 1.3032455 1.4948993       1   3.85 91.67899          45 78.30539       1.19  10.9
##    current_diagnosis end_date_type   race ethnicity immunofix2    sex
##               <char>        <char> <char>    <char>     <char> <char>
## 1:               SMM     Diagnosis  White      <NA>        IgG Female
## 2:               SMM     Diagnosis  White      <NA>        IgG Female
## 3:               SMM     Diagnosis  White      <NA>        IgG Female
## 4:               SMM     Diagnosis  White      <NA>        IgG Female
## 5:               SMM     Diagnosis  White      <NA>        IgG Female
```

# R session
As a reference, we report the R version we used:


``` r
R.Version()$version.string
```

```
## [1] "R version 4.4.2 (2024-10-31)"
```

