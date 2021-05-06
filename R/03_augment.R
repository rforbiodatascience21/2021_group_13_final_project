# Clear workspace ---------------------------------------------------------
#rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------

patient_1 <- read_csv(
  file = "Data/02_patient_1I.csv.gz")
patient_1 <-
  patient_1 %>% 
  mutate(Patient_ID = "1I",
         .before=Cell_Barcode) %>% 
  select(-1) %>% 
  slice(-1)

# Patient_2
patient_2 <- read_csv(
  file = "Data/02_patient_2C.csv.gz")
patient_2 <-
  patient_2 %>% 
  mutate(Patient_ID = "2C",
         .before=Cell_Barcode) %>% 
  select(-1) %>% 
  slice(-1)

# Patient_3
patient_3 <- read_csv(
  file = "Data/02_patient_3C.csv.gz")
patient_3 <-
  patient_3 %>% 
  mutate(Patient_ID = "3C",
         .before=Cell_Barcode) %>% 
  select(-1) %>% 
  slice(-1)

# Patient_4
patient_4 <- read_csv(
  file = "Data/02_patient_4CO.csv.gz")
patient_4 <-
  patient_4 %>% 
  mutate(Patient_ID = "4CO",
         .before=Cell_Barcode) %>% 
  select(-1) %>% 
  slice(-1)

# Patient_5
patient_5 <- read_csv(
  file = "Data/02_patient_5CO.csv.gz")
patient_5 <-
  patient_5 %>% 
  mutate(Patient_ID = "5CO",
         .before=Cell_Barcode) %>% 
  select(-1) %>% 
  slice(-1)

# Patient_6
patient_6 <- read_csv(
  file = "Data/02_patient_6I.csv.gz")
patient_6 <-
  patient_6 %>% 
  mutate(Patient_ID = "6I",
         .before=Cell_Barcode) %>% 
  select(-1) %>% 
  slice(-1)

#Combining patients--------------------------------------------------------

data <- bind_rows(
  patient_1,
  patient_2,
  patient_3,
  patient_4,
  patient_5,
  patient_6,
  .id = "patient_ID")

rm(patient_1,patient_2,patient_3,patient_4,patient_5,patient_6)

#Loading metadata and getting cell types--------------------------------------

metadata <-
  read_csv("Data/_raw/metatable.csv") %>% 
  tibble()

# Wrangle data ------------------------------------------------------------



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
#write_tsv(data/

# Remove Data -------------------------------------------------------------

rm(gravier_data_clean)

