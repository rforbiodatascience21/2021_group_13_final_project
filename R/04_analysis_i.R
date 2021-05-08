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

#cells_count_plot <- 
  cells_per_group_and_patient %>%
  ggplot(aes(x=group,y=total)) +
  geom_col(aes(fill = fct_reorder2(as.factor(Patient_ID),total,group))) +
  scale_fill_viridis_d("Patients") +
  theme_minimal() +
  ylab("count") +
  labs(
    title = "Cell Counts per Group and Patient after Filtering"
  )

### How many cells are there of each type in each group 
  
sum_per_group <-
    data %>% 
    group_by(group) %>% 
    summarise(total_count=n())
  
cells_type_per_group <-
  data %>% 
  group_by(group,CellType_Category) %>% 
  summarise(cell_count=n()) %>% 
  pivot_wider(
    names_from = CellType_Category,
    values_from = cell_count
  ) %>% 
  ungroup() %>% 
  left_join(
    sum_per_group
  ) %>% 
  pivot_longer(
    Endothelial:Stromal,
    names_to = "CellType_Category",
    values_to = "cell_count"
  )

### How are cell types distributed

# cells_type_plot
  cells_type_per_group %>% 
    ggplot(aes(x=group,y=cell_count/total_count)) +
    geom_col(aes(fill = fct_reorder2(as.factor(CellType_Category),cell_count,group))) +
    scale_fill_viridis_d("Type",alpha=0.6) +
    theme_minimal() +
    coord_flip() +
    ylab("relative count") +
    labs(title = "Distribution of Cell Types per Group") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
#cells_type_plot
  

# Plot Data ---------------------------------------------------------------


# Save Plots --------------------------------------------------------------


# Remove Data -------------------------------------------------------------

