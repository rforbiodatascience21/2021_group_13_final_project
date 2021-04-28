# Clear workspace ---------------------------------------------------------
rm(list = ls())

### Provided example code here

# Load libraries ----------------------------------------------------------
library(tidyverse)

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
GSM4058963_025I <- read_rds(file = "Data/_raw/GSM4058963_025I.dgecounts.rds")


# Wrangle data ------------------------------------------------------------
# We only want to look at the exons.
umicount_exon_reads <- GSM4058963_025I %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()


# Write data --------------------------------------------------------------
write_tsv(umicount_exon_reads,
          path = "data/01_data_025I.tsv.gz")

# Remove Data -------------------------------------------------------------
rm(umicount_exon_reads, GSM4058963_025I)

