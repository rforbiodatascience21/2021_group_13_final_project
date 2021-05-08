# Clear workspace ---------------------------------------------------------
#rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("viridis")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------

data <-
  read_csv("Data/03_data.csv")


# Model Data --------------------------------------------------------------

# maing some summary statistics

### first lets reorder the data according to groups

  

### how many cells are left of each patient

cells_per_group_and_patient <-
  data %>% 
  group_by(group,Patient_ID) %>% 
  summarise(total=n())


#levels=c("003C","465C","065CO","207CO","025I","177I")

#cells_plot <- 
  cells_per_group_and_patient %>%
  ggplot(aes(x=group,y=total)) +
  geom_col(aes(fill = fct_reorder2(as.factor(Patient_ID),total,group))) +
  scale_fill_viridis_d("Patients") +
  theme_minimal() +
  ylab("count") +
  labs(
    title = "Cell Counts per Group and Patient after Filtering"
  )

# Plot Data ---------------------------------------------------------------


# Save Plots --------------------------------------------------------------


# Remove Data -------------------------------------------------------------

