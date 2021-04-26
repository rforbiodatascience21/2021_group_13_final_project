# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
gravier_data <- load("data/_raw/gravier.RData")


# Wrangle data ------------------------------------------------------------
gravier_data_x <- gravier %>% pluck("x") %>% as.tibble()
gravier_data_y <- gravier %>% pluck("y") %>% as.tibble()

# Write data --------------------------------------------------------------
write_tsv(gravier_data_x,
          path = "data/01_gravier_data_x.tsv.gz")
write_tsv(gravier_data_y,
          path = "data/01_gravier_data_y.tsv.gz")

# Remove Data -------------------------------------------------------------


rm(gravier_data_y,gravier_data_x,gravier_data,gravier)
