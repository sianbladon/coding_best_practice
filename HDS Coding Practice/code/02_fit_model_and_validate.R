###
### Fit model
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

### Extract one of the imputed datasets
dat_devel_impute_m1 <- mice::complete(dat_devel_impute, action = 1)
str(dat_devel_impute_m1)

### Check no missingness
var_miss <- colnames(dat_devel_impute_m1)[apply(dat_devel_impute_m1, 2, function(x) {any(is.na(x))})]
var_miss

### Fit model
fit_object <- coxph(Surv(time,status) ~ age + sex + chol + trt + platelet + bili + albumin  + protime + edema, data = dat_devel_impute_m1)
bhaz_object <- basehaz(fit_object, centered = TRUE)

### Save research objects
saveRDS(fit_object, "data/fit_object.rds")
saveRDS(bhaz_object, "data/bhaz_object.rds")

### Define validation dataset to be the same as development dataset
### This will lead to optimistic estimation of performance, this exact process should not be replicated in practice.
dat_valid <- dat_devel_impute_m1

### Evaluate calibration of this model using each of the functions that have been written
## Proportional hazards approach
calib_ph <- est_calib_ph(data = dat_valid, fit  = fit_object, bhaz = bhaz_object, t = t_eval, nk = 3)
## IPCW approach
calib_ipcw <- est_calib_ipcw(data = dat_valid, fit  = fit_object, bhaz = bhaz_object, t = t_eval, nk = 3,
                             cens_formula =
                               as.formula("Surv(cens_time, cens_indicator) ~ age + sex + chol + trt + platelet + bili + albumin  + protime + edema"))
## Pseudo-value approach
calib_pv <- est_calib_pv(data = dat_valid, fit  = fit_object, bhaz = bhaz_object, t = t_eval, nk = 3, split_n_groups = 3)

### Save research objects
saveRDS(calib_ph, "data/calib_ph.rds")
saveRDS(calib_ipcw, "data/calib_ipcw.rds")
saveRDS(calib_pv, "data/calib_pv.rds")

### Combine into single plot
plot(calib_ph[["plot"]])
plot(calib_ipcw[["plot"]])
plot(calib_pv[["plot"]])

### Save Figure
png("figures/calibrate_pbc_m1.png", res = 300, width = 6, height = 6, unit = "in")
gridExtra::grid.arrange(calib_ph[["plot"]],
                        calib_ipcw[["plot"]],
                        calib_pv[["plot"]],
                        layout_matrix = base::matrix(seq_len(4), nrow = 2, ncol = 2, byrow = TRUE))
dev.off()
