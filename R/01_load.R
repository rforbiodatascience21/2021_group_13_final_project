# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(tidyverse)

# Tidyseurat is needed for conversion of Matrix object to tibble
library(tidyseurat)

# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data --------------------------------------------------------------------
# Patient_1
patient_1 <- read_rds(file = "Data/_raw/GSM4058963_025I.dgecounts.rds")

# Get gene identifiers from first patient
Genes <- patient_1 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>%   
  pluck("Dimnames") %>% 
  pluck(1) %>% 
  as_tibble()

write_csv(Genes,
          file = "Data/Genes_1")
rm(Genes)


# downsample in load as we cannot work with such big data on PCs
patient_1 <- patient_1 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble() %>% 
  select(c(1:2000))


# Repeat operations for remaining subjects -------------------------------------
# Patient_2                ---------------------------
patient_2 <- read_rds(file = "Data/_raw/GSM4058907_465C.dgecounts.rds")

Genes <- patient_2 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>%   
  pluck("Dimnames") %>% 
  pluck(1) %>% 
  as_tibble()

write_csv(Genes,
          file = "Data/Genes_2")
rm(Genes)

patient_2 <- patient_2 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble() %>% 
  select(c(1:2000))


# Patient_3                ---------------------------
patient_3 <- read_rds(file = "Data/_raw/GSM4058921_003C.dgecounts.rds")

Genes <- patient_3 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>%   
  pluck("Dimnames") %>% 
  pluck(1) %>% 
  as_tibble()

write_csv(Genes,
          file = "Data/Genes_3")
rm(Genes)

patient_3 <- patient_3 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble() %>% 
  select(c(1:2000))


# Patient_4                ---------------------------
patient_4 <- read_rds(file = "Data/_raw/GSM4058936_207CO.dgecounts.rds")

Genes <- patient_4 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>%   
  pluck("Dimnames") %>% 
  pluck(1) %>% 
  as_tibble()

write_csv(Genes,
          file = "Data/Genes_4")
rm(Genes)

patient_4 <- patient_4 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble() %>% 
  select(c(1:2000))


# Patient_5                ---------------------------
patient_5 <- read_rds(file = "Data/_raw/GSM4058944_056CO.dgecounts.rds")

Genes <- patient_5 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>%   
  pluck("Dimnames") %>% 
  pluck(1) %>% 
  as_tibble()

write_csv(Genes,
          file = "Data/Genes_5")
rm(Genes)

patient_5 <- patient_5 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble() %>% 
  select(c(1:2000))


# Patient_6                ---------------------------
patient_6 <- read_rds(file = "Data/_raw/GSM4058977_177I.dgecounts.rds")

Genes <- patient_6 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>%   
  pluck("Dimnames") %>% 
  pluck(1) %>% 
  as_tibble()

write_csv(Genes,
          file = "Data/Genes_6")
rm(Genes)

patient_6 <- patient_6 %>% 
  pluck("umicount") %>% 
  pluck("exon") %>% 
  pluck("all") %>% 
  as_tibble() %>% 
  select(c(1:2000))


# Write data -------------------------------------------------------------------

# Write each downsampled patient into zipped csv
# Patient_1
write.csv(patient_1,
          file = "Data/01_patient_1I.csv.gz")

# Patient_2
write.csv(patient_2,
          file = "Data/01_patient_2C.csv.gz")

# Patient_3
write.csv(patient_3,
          file = "Data/01_patient_3C.csv.gz")

# Patient_4
write.csv(patient_4,
          file = "Data/01_patient_4CO.csv.gz")

# Patient_5
write.csv(patient_5,
          file = "Data/01_patient_5CO.csv.gz")

# Patient_6
write.csv(patient_6,
          file = "Data/01_patient_6I.csv.gz")

# Remove Data -------------------------------------------------------------
rm(patient_1, 
   patient_2, 
   patient_3, 
   patient_4, 
   patient_5, 
   patient_6)