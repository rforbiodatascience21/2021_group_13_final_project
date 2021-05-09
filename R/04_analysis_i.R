# Clear workspace ---------------------------------------------------------
#rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("viridis")
library("broom")
library("nnet")


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
    mutate(
      Counts=normalise(Counts),
      group=case_when(
        group=="Control"~0,
        group=="COPD"~1
      )
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
    cols = 9:ncol(data),
    names_to = "Gene",
    values_to = "Counts"
  ) %>% 
  select(group,Gene,Counts) %>%
  mutate(
    Counts=normalise(Counts),
    group=case_when(
      group=="Control"~0,
      group=="IPF"~1
    )
  ) %>% 
  filter(group != "IPF") %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()


COPD_against_IPF <-
  data %>% 
  filter(Subclass_Cell_Identity=="Ciliated") %>% 
  pivot_longer(
    cols = 9:ncol(data),
    names_to = "Gene",
    values_to = "Counts"
  ) %>% 
  select(group,Gene,Counts) %>%
  mutate(
    Counts=normalise(Counts),
    group=case_when(
      group=="COPD"~0,
      group=="IPF"~1
    )
  ) %>% 
  filter(group != "Control") %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()

multinomial <-
  data %>% 
  filter(Subclass_Cell_Identity=="Ciliated") %>% 
  pivot_longer(
    cols = 9:ncol(data),
    names_to = "Gene",
    values_to = "Counts"
  ) %>% 
  select(group,Gene,Counts) %>%
  mutate(
    Counts=normalise(Counts),
    group=case_when(
      group=="Control"~0,
      group=="COPD"~1,
      group=="IPF"~2
    )
  ) %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()

  

# now delete data as its not needed (memory) AND there is a column in the nested dataframes called data

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
                  conf.int=TRUE)) %>% 
  unnest(tidy) %>% 
  filter(term!="(Intercept)") %>% 
  mutate(adj_p.value=
           p.adjust(
             p.value,
             method="bonferroni"
           )) %>% 
  mutate(identified_as=
           case_when(
             p.value<0.05~"significant",
             TRUE~"unsignificant")) %>% 
  mutate(adj_identified_as=
           case_when(
             adj_p.value<0.05~"significant",
             TRUE~"unsignificant"))

# wihout correction 2,042 "significant"
# with bonferroni/holm/fdr correction: 0

IPF_model <-
  ciliated_IPF %>% 
  mutate(mdl=map(data,
                 ~glm(group~Counts,
                      data=.x,
                      family=binomial(link="logit"))),
         tidy=map(mdl,
                  tidy,
                  conf.int=TRUE)) %>%
  unnest(tidy) %>% 
  filter(term!="(Intercept)") %>% 
  mutate(adj_p.value=
           p.adjust(
             p.value,
             method="bonferroni"
           )) %>% 
  mutate(adj_identified_as=case_when(
    p.value<0.05~"significant",
    TRUE~"unsignificant")) %>% 
  mutate(adj_identified_as=
           case_when(
             adj_p.value<0.05~"significant",
             TRUE~"unsignificant"))

COPD_against_IPF_model <-
  COPD_against_IPF %>% 
  mutate(mdl=map(data,
                 ~glm(group~Counts,
                      data=.x,
                      family=binomial(link="logit"))),
         tidy=map(mdl,
                  tidy,
                  conf.int=TRUE)) %>%
  unnest(tidy) %>% 
  filter(term!="(Intercept)") %>% 
  mutate(adj_p.value=
           p.adjust(
             p.value,
             method="bonferroni"
           )) %>% 
  mutate(identified_as=case_when(
    p.value<0.05~"significant",
    TRUE~"unsignificant")) %>% 
  mutate(adj_identified_as=
           case_when(
             adj_p.value<0.05~"significant",
             TRUE~"unsignificant"))


multinomial_model <-  
  multinomial %>% 
  mutate(mdl=map(data,
                 ~multinom(group~Counts,
                      data=.x)),
         tidy=map(mdl,
                  tidy,
                  conf.int=TRUE)) %>%
  unnest(tidy) %>% 
  filter(term!="(Intercept)") %>% 
  mutate(adj_p.value=
           p.adjust(
             p.value,
             method="bonferroni"
           )) %>% 
  mutate(identified_as=case_when(
    p.value<0.05~"significant",
    TRUE~"unsignificant")) %>% 
  mutate(adj_identified_as=
           case_when(
             adj_p.value<0.05~"significant",
             TRUE~"unsignificant"))

p_without_correction <-
  bind_rows(
    COPD_model %>% 
      filter(p.value),
    IPF_model %>% 
      filter(p.value),
    COPD_against_IPF_model %>% 
      filter(p.value),
    multinomial_model %>% 
      filter(p.value),
  )

#significant_genes_COPD <-
  multinomial_model %>% 
  ggplot(aes(x=identified_as)) +
  geom_bar()
  
  
  
# Save Plots --------------------------------------------------------------


# Remove Data -------------------------------------------------------------

