# Clear workspace --------------------------------------------------------------
#rm(list = ls())

# Load libraries ---------------------------------------------------------------
library("tidyverse")

# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data --------------------------------------------------------------------
# Patient_1
patient_1 <- read_csv(file = "Data/02_patient_1I.csv.gz")

patient_1 <- patient_1 %>% 
  mutate(Patient_ID = "025I",
         .before = Cell_Barcode) %>% 
  select(Patient_ID:ncol(patient_1)) %>% 
  filter(Cell_Barcode != "X1")

# Patient_2                ---------------------------
patient_2 <- read_csv(file = "Data/02_patient_2C.csv.gz")

patient_2 <- patient_2 %>% 
  mutate(Patient_ID = "465C",
         .before = Cell_Barcode) %>% 
  select(Patient_ID:ncol(patient_2)) %>% 
  filter(Cell_Barcode != "X1")

# Patient_3                ---------------------------
patient_3 <- read_csv(file = "Data/02_patient_3C.csv.gz")

patient_3 <- patient_3 %>% 
  mutate(Patient_ID = "003C",
         .before = Cell_Barcode) %>% 
  select(Patient_ID:ncol(patient_3)) %>% 
  filter(Cell_Barcode != "X1")

# Patient_4                ---------------------------
patient_4 <- read_csv(file = "Data/02_patient_4CO.csv.gz")

patient_4 <- patient_4 %>% 
  mutate(Patient_ID = "207CO",
         .before = Cell_Barcode) %>% 
  select(Patient_ID:ncol(patient_4)) %>% 
  filter(Cell_Barcode != "X1")

# Patient_5                ---------------------------
patient_5 <- read_csv(file = "Data/02_patient_5CO.csv.gz")

patient_5 <- patient_5 %>% 
  mutate(Patient_ID = "056CO",
         .before = Cell_Barcode)%>% 
  select(Patient_ID:ncol(patient_5)) %>% 
  filter(Cell_Barcode != "X1")

# Patient_6                ---------------------------
patient_6 <- read_csv(file = "Data/02_patient_6I.csv.gz")

patient_6 <- patient_6 %>% 
  mutate(Patient_ID = "177I",
         .before = Cell_Barcode) %>% 
  select(Patient_ID:ncol(patient_6)) %>% 
  filter(Cell_Barcode != "X1")


# Combining patients------------------------------------------------------------
data <- bind_rows(
  patient_1,
  patient_2,
  patient_3,
  patient_4,
  patient_5,
  patient_6)


# Replacing NAs (gene not present in the patient but in another) by 0
data <- data %>%
  mutate(across(everything(),
           ~replace_na(.x,
                       0))) %>% 
  select(-starts_with("ENSG"))


# Clearing environment
rm(patient_1, 
   patient_2,
   patient_3,
   patient_4,
   patient_5,
   patient_6)


# Loading metadata -------------------------------------------------------------
metadata <- read_tsv(
  "Data/_raw/GSE136831_AllCells.Samples.CellType.MetadataTable.txt") %>% 
  tibble() %>% 
  select(CellBarcode_Identity,
         nUMI,
         nGene,
         CellType_Category,
         Subclass_Cell_Identity,
         Subject_Identity) %>% 
  rename(Patient_ID = Subject_Identity,
         Cell_Barcode = CellBarcode_Identity) %>% 
  filter(Patient_ID == "025I" |
           Patient_ID == "465C" |
           Patient_ID == "003C" |
           Patient_ID == "207CO"|
           Patient_ID == "056CO"|
           Patient_ID == "177I") %>% 
  mutate(Cell_Barcode = str_replace_all(Cell_Barcode,
                                        "(.+_)(.+)",
                                        "\\2"))


# Join metadata and patient data------------------------------------------------
data <- left_join(data,
                  metadata,
                  by = c("Patient_ID",
                         "Cell_Barcode")) %>% 
  relocate("nGene",
           "nUMI",
           "CellType_Category",
           "Subclass_Cell_Identity",
           .after = "Cell_Barcode") 


# Introduce group label and make that and Patient_ID factors
data <-
  data %>% 
  mutate(group = 
           factor(
             case_when(
             str_detect(Patient_ID,"I") ~ "IPF",
             str_detect(Patient_ID,"CO") ~ "COPD",
             TRUE ~ "Control")),         
         .after = Patient_ID,
         Patient_ID = factor(Patient_ID))

rm(metadata)


# Now there are quite some NAs as the metadata apparently contains only cells 
# filtered by the authors.
# We will filter further to second-check this filtering. Thus if it works well, 
# the NAs will mostly be gone

# Wrangle data -----------------------------------------------------------------

# Combine the patients' gene expression to the metadata data of the cells.
# Run patient_slicer, meta_slicer, combiner and binder

# If there are cells with fewer than 2000 transcripts recorded, these cells are
# filtered out from the 10000 starting cell count per patient
data <- data %>% 
  filter(nUMI > 2000)


# Filtering out the cells where the transcripts of the mitochondrial genes 
# represent more than 20% of the sum of the transcripts for a cell
# Run mito_filter
mt_selection <- select(data, 
                       starts_with("MT"))

mt_sum <- mt_selection %>% 
  mutate(mito_sum = rowSums(mt_selection))

data <- data %>%
  mutate(select(mt_sum,
                mito_sum),
         .after = Cell_Barcode) %>%
  filter(mito_sum / nUMI < 0.2)


# Now check if there are any NAs remaining
data %>% 
  filter(
  across(.cols = everything(),
         .fns = ~ is.na(.x)))

# --> no NAs left

#Breaking the complete data set into patient groups

Control <- data%>%
  filter(Patient_ID !="025I")%>%
  filter(Patient_ID !="207CO")%>%
  filter(Patient_ID !="056CO")%>%
  filter(Patient_ID !="177I")

COPD <- data%>%
  filter(Patient_ID !="025I")%>%
  filter(Patient_ID !="465C")%>%
  filter(Patient_ID !="003C")%>%
  filter(Patient_ID !="177I")

IPF <- data%>%
  filter(Patient_ID !="465C")%>%
  filter(Patient_ID !="207CO")%>%
  filter(Patient_ID !="056CO")%>%
  filter(Patient_ID !="003C")


# Write data -------------------------------------------------------------------
write_csv(data, 
          file = "Data/03_data.csv")

write_csv(Control, 
          file = "Data/03_Control.csv")
write_csv(COPD, 
          file = "Data/03_COPD.csv")
write_csv(IPF, 
          file = "Data/03_IPF.csv")



# Remove Data ------------------------------------------------------------------
rm(data,
   mt_selection,
   mt_sum)
