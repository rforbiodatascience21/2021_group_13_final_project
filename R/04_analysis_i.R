# Clear workspace ---------------------------------------------------------
#rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("viridis")
library("broom")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------

data <-
  read_csv("Data/03_data.csv")


# Model and Plot Data --------------------------------------------------------------

# making some summary statistics

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


  ### models ###
  
# Take only ciliated cells,
# pivot longer to get genes as rows for modelling
# normalise the Counts, !!!!!WE NEED TO DO THAT BEFORE!!!! because I think otherwise it will be normalised according to all levels, not the specific cells one

  
# Proposed solution, cant run this now
  
normalised <-
  data %>% 
  select(9:ncol(.)) %>% 
  map(normalise)
  
  
  
ciliated_COPD <-
    data %>% 
    filter(Subclass_Cell_Identity=="Ciliated") %>% 
    pivot_longer(
      cols = 9:ncol(data),
      names_to = "Gene",
      values_to = "Counts"
    ) %>% 
    select(group,Gene,Counts) %>%
#%>%
    #sample_n(250) %>% 
    mutate(
      Counts=normalise(Counts)
    ) %>% 
  filter(group != "IPF") %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()

# same for IPF

ciliated_IPF <-
  data %>% 
  filter(Subclass_Cell_Identity=="Ciliated") %>% 
  pivot_longer(
    cols = 8:ncol(data),
    names_to = "Gene",
    values_to = "Counts"
  ) %>% 
  select(group,Gene,Counts) %>%
  #%>%
  #sample_n(250) %>% 
  mutate(
    Counts=normalise(Counts)
  ) %>% 
  filter(group != "IPF") %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()
  

# now delete data as its not needed AND there is a column in the nested dataframes called data

rm(data)

# make two models

COPD_model <-
  ciliated_COPD %>% 
  mutate(mdl=map(data,
                 ~glm(group~Counts,
                      data=.x,
                      family=binomial(link="logit"))),
         tidy=map(mdl,
                  tidy,
                  conf.int=TRUE))

IPF_model <-
  ciliated_IPF %>% 
  mutate(mdl=map(data,
                 ~glm(group~Counts,
                      data=.x,
                      family=binomial(link="logit"))),
         tidy=map(mdl,
                  tidy,
                  conf.int=TRUE))
  
# Save Plots --------------------------------------------------------------


# Remove Data -------------------------------------------------------------

