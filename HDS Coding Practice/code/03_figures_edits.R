###
### Fit model
###

### Clear workspace
rm(list=ls())
### IF WANT TO MANUALLY SET WORKING DIRECTORY, DO IT HERE. RELATIVE FILEPATHS FROM THIS POINT ONWARDS.
#setwd("/path/to/root/directory/")
getwd()

### Load libraries
library(survival)
library(mice)
library(prodlim)
library(Cairo)

### Read in research objects
calib_ph <- readRDS("data/calib_ph.rds")
calib_ipcw <- readRDS("data/calib_ipcw.rds")
calib_pv <- readRDS("data/calib_pv.rds")

### Save Figure
png("figures/calibrate_pbc_m1_edit.png", res = 300, width = 6, height = 2, unit = "in")
gridExtra::grid.arrange(calib_ph[["plot"]],
                        calib_ipcw[["plot"]],
                        calib_pv[["plot"]],
                        layout_matrix = base::matrix(seq_len(3), nrow = 1, ncol = 3, byrow = TRUE))
dev.off()
