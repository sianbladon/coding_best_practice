###
### Fit model over different imputed datasets
###

### Clear workspace
rm(list=ls())
### IF WANT TO MANUALLY SET WORKING DIRECTORY, DO IT HERE. RELATIVE FILEPATHS FROM THIS POINT ONWARDS.
#setwd("/path/to/root/directory/")
getwd()

### Read in research objects
dat_devel_impute <- readRDS("data/dat_devel_impute.rds")

### Load functions and libraries
source("R/functions.R")

### Load libraries
library(survival)
library(mice)
library(prodlim)
library(Cairo)

### Define t_eval
t_eval <- 5*365.25

### Write a function to
### 1) Extract one of the imputed datasets
### 2) Fit model
### 3) Run validation
### 4) Save calibration plots

run_analyses <- function(m){

  ### Extract one of the imputed datasets
  dat_devel_impute <- mice::complete(dat_devel_impute, action = m)

  ### Fit model
  fit_object <- coxph(Surv(time,status) ~ age + sex + chol + trt + platelet + bili + albumin  + protime + edema, data = dat_devel_impute)
  bhaz_object <- basehaz(fit_object, centered = TRUE)

  ### Save research objects
  saveRDS(fit_object, paste("data/fit_object", m, ".rds", sep = ""))
  saveRDS(bhaz_object, paste("data/bhaz_object", m, ".rds", sep = ""))

  ### Define validation dataset to be the same as development dataset
  ### This will lead to optimistic estimation of performance, this exact process should not be replicated in practice.
  dat_valid <- dat_devel_impute

  ### Evaluate calibration of this model using each of the functions that have been written
  ## Proportional hazards approach
  calib_ph <- est_calib_ph(data = dat_valid, fit  = fit_object, bhaz = bhaz_object, t = t_eval, nk = 3)
  ## IPCW approach
  calib_ipcw <- est_calib_ipcw(data = dat_valid, fit  = fit_object, bhaz = bhaz_object, t = t_eval, nk = 3,
                               cens_formula =
                                 as.formula("Surv(cens_time, cens_indicator) ~ age + sex + chol + trt + platelet + bili + albumin  + protime + edema"))
  ## Pseudo-value approach
  calib_pv <- est_calib_pv(data = dat_valid, fit  = fit_object, bhaz = bhaz_object, t = t_eval, nk = 3, split_n_groups = 3)

  ### Save Figure
  png(paste("figures/calibrate_pbc_m", m, ".png", sep = ""), res = 300, width = 6, height = 2, unit = "in")
  gridExtra::grid.arrange(calib_ph[["plot"]],
                          calib_ipcw[["plot"]],
                          calib_pv[["plot"]],
                          layout_matrix = base::matrix(seq_len(3), nrow = 1, ncol = 3, byrow = TRUE))
  dev.off()

}

### Run analyses for each imputed dataset (m = 1,2,3,4,5)
lapply(c(1,2,3,4,5), run_analyses)
