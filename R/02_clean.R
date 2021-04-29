# Clear workspace ---------------------------------------------------------
#rm(list = ls())

### Provided example code here

# Load libraries ----------------------------------------------------------
library("tidyverse")
library(patchwork)


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
patients <-
  list("GSM4058963_025I","GSM4058942_8CO")

data <- map(patients,load_umicounts)

### Wrangle data ------------------------------------------------------------

# check if there are rows with only zeros == cells with no expression, would correspond to empty droplets
data <- map(data,remove_zero_rows)

# no there arent. now check if there are any with particularly high umi count which could correspond to droplets with more than one cell

# get all the counts
ranges <- map(data,get_range)

# visualize to help understand what could be filtered out

patient1 <-
  ggplot(
    pluck(ranges,1),
    aes(row_sum)) + 
  geom_histogram()

# there is a big range

hist_data1 <-
  ranges %>% 
  pluck(1) %>% 
  filter(row_sum>10000)

# the following graph needs improvement

patient1 <-
  ggplot(hist_data1,aes(row_sum)) + 
  geom_freqpoly() +
  xlim(min(hist_data1),max(hist_data1))

# lets next find good values to cut off (the small and the high ones)

# Write data --------------------------------------------------------------


# Remove Data -------------------------------------------------------------


rm(variables created)