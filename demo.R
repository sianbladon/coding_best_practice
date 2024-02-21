# demo

install.packages("NHSRdatasets", lib = "C:/Program Files/R/R-4.3.1/library/")

library(NHSRdatasets)
library(tidyverse)

data("LOS_model")

head(LOS_model)

write_csv(LOS_model, "data/trust_los.csv")

LOS_model <- as.data.frame(LOS_model)



library(here)
#ggsave(here::here("figs", paste0(org_name, ".png")))


# cleaning data ---------------------------------------------------------------


# demographic characteristics summary -------------------------------------


