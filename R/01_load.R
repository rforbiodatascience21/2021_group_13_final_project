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
  read_rds(file = "Data/_raw/GSM4058963_025I.dgecounts.rds")

patient_1 <- 
  patient_1 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# downsample

patient_1 <- 
  sample_n(patient_1,100)

#-------------------------------------------------------------------------

patient_2 <- 
  read_rds(file = "Data/_raw/GSM4058907_465C.dgecounts.rds")

patient_2 <- 
  patient_2 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# downsample

patient_2 <- 
  sample_n(patient_2,100)

#-------------------------------------------------------------------------

patient_3 <- 
  read_rds(file = "Data/_raw/GSM4058921_003C.dgecounts.rds")

patient_3 <- 
  patient_3 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# downsample

patient_3 <- 
  sample_n(patient_3,100)

#-------------------------------------------------------------------------


patient_4 <- 
  read_rds(file = "Data/_raw/GSM4058936_207CO.dgecounts.rds")

patient_4 <- 
  patient_4 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# downsample

patient_4 <- 
  sample_n(patient_4,100)

#-------------------------------------------------------------------------


patient_5 <- 
  read_rds(file = "Data/_raw/GSM4058944_056CO.dgecounts.rds")

patient_5 <- 
  patient_5 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# downsample

patient_5 <- 
  sample_n(patient_5,100)

#-------------------------------------------------------------------------

patient_6 <- 
  read_rds(file = "Data/_raw/GSM4058977_177I.dgecounts.rds")

patient_6 <- 
  patient_6 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble()

# downsample

patient_6 <- 
  sample_n(patient_6,100)

# Write data --------------------------------------------------------------
# Write each downsampled patient into zipped csv
# Patient_1
write.csv(patient_1,
          path = "Data/01_patient_1.csv.gz")

# Patient_2
write.csv(patient_2,
          path = "Data/01_patient_2.csv.gz")

# Patient_3
write.csv(patient_3,
          path = "Data/01_patient_3.csv.gz")

# Patient_4
write.csv(patient_4,
          path = "Data/01_patient_4.csv.gz")

# Patient_5
write.csv(patient_5,
          path = "Data/01_patient_5.csv.gz")

# Patient_6
write.csv(patient_6,
          path = "Data/01_patient_6.csv.gz")

# Remove Data -------------------------------------------------------------
rm(patient_1,patient_2,patient_3,patient_4,patient_5,patient_6)

