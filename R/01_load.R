# Clear workspace ---------------------------------------------------------
rm(list = ls())

### Provided example code here

# Load libraries ----------------------------------------------------------
library(tidyverse)

# tidyseurat is needed (at least for Karl's environment) for conversion of Matrix object to tibble
library(tidyseurat)

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
patients <-
  list("GSM4058963_025I","GSM4058942_8CO")

data <- 
  make_list_of_patient_data(patients)

# Wrangle data ------------------------------------------------------------
# We only want to look at the exons.

data <- 
  map(data,get_exon_umicounts)

# Write data --------------------------------------------------------------
walk2(data,patients,write_umicounts)

# Remove Data -------------------------------------------------------------
rm(patients,data)

