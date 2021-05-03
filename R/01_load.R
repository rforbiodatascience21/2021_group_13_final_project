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


# load first patient
patient_1 <- 
  read_rds(file = "/Data/_raw/GSM4058963_025I.dgecounts.rds")

patient_1 <- 
  data %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# sample

patient_1 <- 
  sample_n(data,100)

# ---------------------


patient_2 <- 
  read_rds(file = "/Data/_raw/GSM4058907_465C.dgecounts.rds")

patient_2 <- 
  data %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# sample

patient_2 <- 
  sample_n(data,100)

#--------------------------

patient_3 <- 
  read_rds(file = "/Data/_raw/GSM4058921_003C.dgecounts.rds")

patient_3 <- 
  data %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# sample

patient_3 <- 
  sample_n(data,100)

#------------------------


patient_4 <- 
  read_rds(file = "/Data/_raw/GSM4058936_207CO.dgecounts.rds")

patient_4 <- 
  data %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# sample

patient_4 <- 
  sample_n(data,100)

#-----------------------------


patient_5 <- 
  read_rds(file = "/Data/_raw/GSM4058944_056CO.dgecounts.rds")

patient_5 <- 
  data %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# sample

patient_5 <- 
  sample_n(data,100)

#---------------------

patient_6 <- 
  read_rds(file = "/Data/_raw/GSM4058977_177I.dgecounts.rds")

patient_6 <- 
  data %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# sample

patient_6 <- 
  sample_n(data,100)

# Wrangle data ------------------------------------------------------------
# We only want to look at the exons.

data <- 
  map(data,get_exon_umicounts)

metadata <-
  

data <- 
  map2(data,1000,sample_n)
  

# Write data --------------------------------------------------------------
walk2(data,patients,write_umicounts)

# Remove Data -------------------------------------------------------------
rm(patients,data)

