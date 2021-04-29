# Clear workspace ---------------------------------------------------------
#rm(list = ls())

### Provided example code here

# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
patients <-
  list("GSM4058963_025I","GSM4058942_8CO")

data <- map(patients,load_umicounts)

# Wrangle data ------------------------------------------------------------
data <- map(data,remove_zero_rows)

# Write data --------------------------------------------------------------


# Remove Data -------------------------------------------------------------


rm(variables created)