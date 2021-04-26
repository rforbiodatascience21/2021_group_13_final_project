# Clear workspace ---------------------------------------------------------
rm(list = ls())

### Provided example code here

# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
dataframe <- load("some/data/example.RData")


# Wrangle data ------------------------------------------------------------
dataframe_wrangle <- dataframe %>% operations


# Write data --------------------------------------------------------------
write_tsv(dataframe_wrangle,
          path = "data/01_data.tsv.gz")


# Remove Data -------------------------------------------------------------


rm(dataframe, dataframe_wrangle)
