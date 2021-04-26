# Clear workspace ---------------------------------------------------------
#rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
gravier_x <- read_tsv("data/01_gravier_data_x.tsv.gz")
gravier_y <- read_tsv("data/01_gravier_data_y.tsv.gz")

# Wrangle data ------------------------------------------------------------
gravier_data_clean <- 
  bind_cols(gravier_y,gravier_x) %>% 
  mutate(outcome = value, 
         .keep="unused",
         .before=1)


# Write data --------------------------------------------------------------
write_tsv(gravier_data_clean,
        path = "data/02_gravier_data_clean.tsv.gz")

# Remove Data -------------------------------------------------------------


rm(gravier_x,gravier_y)