# Clear workspace ---------------------------------------------------------
#rm(list = ls())

# Load libraries ----------------------------------------------------------
library(tidyverse)
# patchwork if we want to plot the boxplots below of more than one patient together
library(patchwork)


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
#Patient 1
patient_1 <- read_csv(
          file = "Data/01_patient_1I.csv.gz")

# Patient_2
patient_2 <- read_csv(
          file = "Data/01_patient_2C.csv.gz")

# Patient_3
patient_3 <- read_csv(
          file = "Data/01_patient_3C.csv.gz")

# Patient_4
patient_4 <- read_csv(
          file = "Data/01_patient_4CO.csv.gz")

# Patient_5
patient_5 <- read_csv(
          file = "Data/01_patient_5CO.csv.gz")

# Patient_6
patient_6 <- read_csv(
          file = "Data/01_patient_6I.csv.gz")

#Loading metadata
meta <-  read_table2("Data/_raw/GSE136831_AllCells.Samples.CellType.MetadataTable.txt", col_names = TRUE)
### Wrangle data ------------------------------------------------------------


#--------getting gene names and cell barcodes----------------

gene_mapping <-
  read_csv("Data/_raw/ids.csv", col_names=c("value","Gene")) %>% 
  tibble()

Genes <- read_csv("Data/Genes")
Genes <-
  Genes %>% 
  left_join(gene_mapping) %>% 
  pluck("Gene")

rm(gene_mapping)

# Patient_1


# transpose so that genes are columnns and rows are cells

patient_1 <-
  patient_1 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)


# rename genes and make nice name for column with barcodes

colnames(patient_1) <- c("Cell_Barcode",Genes)

#-----now do the same for 5 other patients

#Patient_2

patient_2 <-
  patient_2 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

colnames(patient_2) <- c("Cell_Barcode",Genes)

#Patient_3

patient_3 <-
  patient_3 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

colnames(patient_3) <- c("Cell_Barcode",Genes)

#Patient_4

patient_4 <-
  patient_4 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

colnames(patient_4) <- c("Cell_Barcode",Genes)

#Patient_5

patient_5 <-
  patient_5 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

colnames(patient_5) <- c("Cell_Barcode",Genes)

#Patient_6

patient_6 <-
  patient_6 %>% 
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

colnames(patient_6) <- c("Cell_Barcode",Genes)


# check if there are rows with only zeros == cells with no expression, would correspond to empty droplets
data <- map(data,remove_zero_rows)

# no there arent. now check if there are any with particularly high umi count which could correspond to droplets with more than one cell

# get all the counts
ranges <- map(data,get_range)

# for the following first EDA work with only the counts as data is so big

# visualize to help understand what could be filtered out

patient1 <-
  ggplot(
    pluck(ranges,1),
    aes(row_sum)) + 
  geom_boxplot()

# there is a big range and a lot of ones

plot_data1 <-
  ranges %>% 
  pluck(1) %>% 
  filter(row_sum>10000)

patient1 <-
  ggplot(hist_data1,aes(row_sum)) + 
  geom_boxplot()

# still quite skewed
# lets next find good values to cut off (the small and the high ones)
# i would not trust cells with less than a few dozen transcripts...lets see above 50.
# lets also determine the upper value of where 90% of the data lie as there are certainly a few outliers coming from questionavle results


nineteeth <- ranges %>% pluck(1) %>% pull(1) %>% quantile(0.9)

patient1_data_refined <- 
  ranges %>% 
  pluck(1) %>% 
  filter(row_sum > 50 & row_sum < nineteeth) 

patient1 <-
  ggplot(patient1_data_refined,aes(row_sum)) + 
  geom_boxplot()

# data is still skewed but we know the variance in read counts in these experiments is high
# this EDA makes me write and use the following functions to get the cleaning streamlined

lower_cutoff <- 50
ranges_refined <- map2(ranges,lower_cutoff,refine_range)
rm(ranges)

### this function should work but my computer cannot handle it (it worked for one in the console) ###
data <- map2(data,ranges_refined,refine_data)


# we need to normalise the rest for our analysis 
#(it was written in the paper that they did that and scaled to 100000 umis per cell)
# i dont get why one should scale to 100000 if one can make it relative right away...
# i guess the question is here do they really compare absolute values
# which makes sense to me because that would certainly tell you something abouts a cell's acticity
# but given all of the uncertainties of the technology i cant imagine that it makes much sense

# I couldnt think of a tidyverse way to do it. as of now im thinking, loop through rows, divide by row_sum
# replace old row with new row, i dont think we want to go that way
# the normalise function does that but its very slow


# next up in the cleaning would be to join the genes expression data with some metadata

# Write data --------------------------------------------------------------


# Remove Data -------------------------------------------------------------


#rm(variables created)
