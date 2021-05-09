# Clear workspace ---------------------------------------------------------
#rm(list = ls())

# Load libraries ----------------------------------------------------------
library(tidyverse)

# Patchwork if we want to plot the boxplots below of more than one patient together
library(patchwork)

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------

# Patient 1
patient_1 <- 
  read_csv(file = "Data/01_patient_1I.csv.gz")


# Patient_2
patient_2 <- 
  read_csv(file = "Data/01_patient_2C.csv.gz")


# Patient_3
patient_3 <- 
  read_csv(file = "Data/01_patient_3C.csv.gz")


# Patient_4
patient_4 <- 
  read_csv(file = "Data/01_patient_4CO.csv.gz")


# Patient_5
patient_5 <- 
  read_csv(file = "Data/01_patient_5CO.csv.gz")


# Patient_6
patient_6 <-
  read_csv(file = "Data/01_patient_6I.csv.gz")


# Combining patients -----------------------------------------------------------

### Wrangle data ---------------------------------------------------------------

# ------------------- getting gene names and cell bar codes --------------------

# Loading metadata and getting cell types --------------------------------------
# Metadata supplied by the authors of the paper("ids.csv")
# This is just done for pretty gene names and does not provide any other benefit
gene_mapping <- 
  read_csv("Data/_raw/ids.csv",
           col_names=c("value","Gene")) %>%
  tibble()
  
  
# Patient_1
# transpose so that genes are columns and rows are cells
patient_1 <-
  patient_1 %>% 
  #select(-X1) %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname,
              values_from = value)

Genes <- read_csv("Data/Genes_1")
Genes <-
  Genes %>% 
  left_join(gene_mapping) %>% 
  pluck("Gene")

# rename genes and make nice name for column with barcodes
colnames(patient_1) <- 
  c("Cell_Barcode",Genes)


# -------------- now do the same for 5 other patients --------------

#Patient_2                ---------------------------
patient_2 <-
patient_2 %>%
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname,
              values_from = value)

Genes <- read_csv("Data/Genes_2")
Genes <-Genes %>% 
  left_join(gene_mapping) %>% 
  pluck("Gene")

colnames(patient_2) <- c("Cell_Barcode",Genes)


#Patient_3                ---------------------------
patient_3 <-patient_3 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname,
              values_from = value)

Genes <- read_csv("Data/Genes_3")
Genes <-Genes %>% 
  left_join(gene_mapping) %>% 
  pluck("Gene")

colnames(patient_3) <- c("Cell_Barcode",Genes)


#Patient_4                ---------------------------
patient_4 <-
  patient_4 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname,
              values_from = value)

Genes <- read_csv("Data/Genes_4")
Genes <-
  Genes %>% 
  left_join(gene_mapping) %>% 
  pluck("Gene")

colnames(patient_4) <- c("Cell_Barcode",Genes)


#Patient_5                ---------------------------
patient_5 <-
  patient_5 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname,
              values_from = value)

Genes <- read_csv("Data/Genes_5")
Genes <-
  Genes %>% 
  left_join(gene_mapping) %>% 
  pluck("Gene")

colnames(patient_5) <- c("Cell_Barcode",Genes)


#Patient_6                ---------------------------
patient_6 <-
  patient_6 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from = rowname, 
              values_from = value)

Genes <- read_csv("Data/Genes_6")
Genes <-
  Genes %>% 
  left_join(gene_mapping) %>% 
  pluck("Gene")

colnames(patient_6) <- c("Cell_Barcode",Genes)

# Write data -------------------------------------------------------------------

# Patient_1
write.csv(patient_1,
          file = "Data/02_patient_1I.csv.gz")

# Patient_2
write.csv(patient_2,
          file = "Data/02_patient_2C.csv.gz")

# Patient_3
write.csv(patient_3,
          file = "Data/02_patient_3C.csv.gz")

# Patient_4
write.csv(patient_4,
          file = "Data/02_patient_4CO.csv.gz")

# Patient_5
write.csv(patient_5,
          file = "Data/02_patient_5CO.csv.gz")

# Patient_6
write.csv(patient_6,
          file = "Data/02_patient_6I.csv.gz")

# Remove Data ------------------------------------------------------------------
rm(patient_1, patient_2, patient_3, patient_4, patient_5, patient_6, Genes, gene_mapping)
