# Clear workspace ---------------------------------------------------------
rm(list = ls())

### Provided example code here

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(Seurat)

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
GSM4058963_025I <- read_rds(file = "Data/_raw/GSM4058963_025I.dgecounts.rds")


# Wrangle data ------------------------------------------------------------
# We only want to look at the exons.
umicount_exon_reads <- GSM4058963_025I %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all")

view(umicount_exon_reads)
umicount_exon_reads %>% 
  select("AAACCTGAGCAGCCTC")

# Write data --------------------------------------------------------------
write_tsv(dataframe_wrangle,
          path = "data/01_data.tsv.gz")


# Remove Data -------------------------------------------------------------
rm(dataframe, dataframe_wrangle)

