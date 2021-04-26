# Clear workspace ---------------------------------------------------------
#rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
gravier_data_clean <- read_tsv(file = "data/02_gravier_data_clean.tsv.gz")


# Wrangle data ------------------------------------------------------------
gravier_data_aug <- 
  gravier_data_clean %>% 
  mutate(outcome = case_when(outcome=="good"~0,
                             outcome=="poor"~1))


# Write data --------------------------------------------------------------
write_tsv(x = gravier_data_aug,
          path = "data/03_gravier_data_aug.tsv.gz")

# Remove Data -------------------------------------------------------------

rm(gravier_data_clean)

