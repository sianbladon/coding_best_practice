###
### Data preperation and imputation
###

### Clear workspace
rm(list=ls())
### IF WANT TO MANUALLY SET WORKING DIRECTORY, DO IT HERE. RELATIVE FILEPATHS FROM THIS POINT ONWARDS.
#setwd("/path/to/root/directory/")
getwd()

###
### Preliminaries
###

### Load functions
source("R/functions.R")

### Load libraries
library(survival)
library(mice)

### Set seed
set.seed(101)

### Load PBC data
dat_devel <- pbc

### Add id variable required for pseudo-value functions
dat_devel$id <- 1:nrow(dat_devel)

### Change status value to be 0/1, not 0/1/2
dat_devel$status <- pmin(dat_devel$status, 1) # 0=censor, 1=transplant, 2=death

### Lets create a variable for censored or not censored
### (if an individual doesn't have an event, they are censored at end of follow up)
### This will be used in the IPCW method
dat_devel$cens_time <- dat_devel$time
dat_devel$cens_indicator <- 1 - dat_devel$status

### Identify missing variables
var_miss <- colnames(dat_devel)[apply(dat_devel, 2, function(x) {any(is.na(x))})]

### Get level of missingness
mice::md.pattern(dat_devel[, var_miss], plot = FALSE)

###
### Impute
###

### NB: THIS IS A BAD EXAMPLE OF HOW TO IMPUTE.
###     Input Parameters should be chosen with more care.
###     See Stef Van Buuren book for more: stefvanbuuren.name/fimd/

### Do a dry run to get predictor matrix
dry_run <- mice::mice(data = dat_devel, m = 1, maxit = 0)

### Save and edit predictor matrix
pred_matrix <- dry_run$predictorMatrix
pred_matrix[,"id"] <- 0
pred_matrix[,"time"] <- 0
pred_matrix[,"status"] <- 0
pred_matrix

### Impute using new predictor matrix
dat_devel_impute <- mice::mice(data = dat_devel, m = 5, maxit = 10, predictorMatrix = pred_matrix)

### Save research objects
saveRDS(dat_devel_impute, "data/dat_devel_impute.rds")
saveRDS(dat_devel, "data/dat_devel.rds")
