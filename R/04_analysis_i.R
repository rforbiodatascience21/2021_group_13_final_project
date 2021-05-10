# Clear workspace --------------------------------------------------------------
#rm(list = ls())


# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("viridis")
library("broom")
library("nnet")
library("patchwork")


# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data --------------------------------------------------------------------
data <-
  read_csv("Data/03_data.csv")


# Model and Plot Data ----------------------------------------------------------
# Making some summary statistics

### First lets reorder the data according to groups
### How many cells are left of each patient
cells_per_group_and_patient <-
  data %>% 
  group_by(group,Patient_ID) %>% 
  summarise(total = n())

cells_count_plot <- 
  cells_per_group_and_patient %>%
  ggplot(aes(x = group,
             y = total)) +
  geom_col(aes(fill = fct_reorder2(as.factor(Patient_ID),
                                   total,group))) +
  scale_fill_viridis_d("Patients") +
  theme_minimal() +
  ylab("count") +
  labs(
    title = "Cell Counts per Group and Patient after Filtering")

  
### How many cells are there of each type in each group
sum_per_group <-
    data %>% 
    group_by(group) %>% 
    summarise(total_count = n())
  
cells_type_per_group <-
  data %>% 
  group_by(group,
           CellType_Category) %>% 
  summarise(cell_count = n()) %>% 
  pivot_wider(names_from = CellType_Category,
              values_from = cell_count) %>% 
  ungroup() %>% 
  left_join(sum_per_group) %>% 
  pivot_longer(Endothelial:Stromal,
               names_to = "CellType_Category",
               values_to = "cell_count")


### How are cell types distributed
cells_type_plot <-
  cells_type_per_group %>% 
    ggplot(aes(x = group,
               y = cell_count / total_count)) +
    geom_col(aes(fill = fct_reorder2(as.factor(CellType_Category),
                                     cell_count,group))) +
    scale_fill_viridis_d("Type", 
                         alpha = 0.6) +
    theme_minimal() +
    coord_flip() +
    ylab("relative count") +
    labs(title = "Distribution of Cell Types per Group") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())



# plot distribution of epithelial
 
subtypes_epithelial <-
  data %>% 
  filter(CellType_Category=="Epithelial") %>%  
  group_by(group,
           Subclass_Cell_Identity) %>% 
  summarise(cell_count = n()) %>% 
  pivot_wider(names_from = Subclass_Cell_Identity,
              values_from = cell_count) %>%
  ungroup() %>% 
  pivot_longer("ATII_High-Surfactants":"Mystery_Disease_Epithelial",
               names_to = "CellType_Subclass",
               values_to = "cell_count") %>% 
  mutate(
    cell_count=
      replace_na(
      cell_count,
      0)
  )

total_subtypes_per_group <-
  subtypes_epithelial %>% 
  group_by(group) %>% 
  summarise(total=
              sum(cell_count))

subtypes_epithelial <-
  left_join(
    subtypes_epithelial,
    total_subtypes_per_group)

plot_epithelial_subtypes <-
  subtypes_epithelial %>% 
    ggplot(aes(x = group,
               y = cell_count / total)) +
    geom_col(aes(fill = fct_reorder2(as.factor(CellType_Subclass),
                                     cell_count,group))) +
    scale_fill_viridis_d("Type", 
                         alpha = 0.6) +
    theme_minimal() +
    coord_flip() +
    labs(title = "Distribution of Subtypes of Ciliated Cells") +
    theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# remove to save memory

rm(cells_per_group_and_patient,
   sum_per_group,
   cells_type_per_group,
   subtypes_epithelial,
   total_subtypes_per_group)

### models ###

# Take only ciliated cells,
# Pivot longer to get genes as rows for modelling
# Normalise the Counts

## make datasets for models  

ciliated <-
  data %>% 
  filter(Subclass_Cell_Identity == "Ciliated")

  
ciliated_COPD <-
    ciliated %>% 
    pivot_longer(cols = TSPAN6:ncol(ciliated),
                 names_to = "Gene",
                 values_to = "Counts") %>% 
    select(group,
           Gene,
           Counts) %>%
  filter(!is.na(Counts)) %>% 
    mutate(Counts = normalise(Counts),
           group = case_when(group == "Control" ~ 0,
                             group == "COPD" ~ 1,
                             TRUE ~ 999)) %>% 
  filter(group != 999) %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()


# same for IPF
ciliated_IPF <-
  ciliated %>% 
  pivot_longer(cols = TSPAN6:ncol(ciliated),
               names_to = "Gene",
               values_to = "Counts") %>% 
  select(group,
         Gene,
         Counts) %>%
  filter(!is.na(Counts)) %>% 
  mutate(Counts = normalise(Counts),
         group = case_when(group == "Control" ~ 0,
                           group == "IPF" ~ 1,
                           TRUE ~ 999)) %>% 
  filter(group != 999) %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()


COPD_against_IPF <-
  ciliated %>% 
  pivot_longer(cols = TSPAN6:ncol(ciliated),
               names_to = "Gene",
               values_to = "Counts") %>% 
  select(group,
         Gene,
         Counts) %>%
  filter(!is.na(Counts)) %>% 
  mutate(Counts = normalise(Counts),
         group = case_when(group == "COPD" ~ 0,
                           group == "IPF" ~ 1,
                           TRUE ~ 999)) %>% 
  filter(group != 999) %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()

multinomial <-
  ciliated %>% 
  pivot_longer(cols = TSPAN6:ncol(ciliated),
               names_to = "Gene",
               values_to = "Counts") %>% 
  select(group,
         Gene,
         Counts) %>%
  filter(!is.na(Counts)) %>% 
  mutate(Counts = normalise(Counts),
         group = case_when(group == "Control" ~ 0,
                           group == "COPD" ~ 1,
                           group == "IPF" ~ 2)) %>% 
  group_by(Gene) %>% 
  nest() %>% 
  ungroup()


# now delete data as its not needed (memory) AND there is a column in the nested dataframes called data
rm(data)

# make four models
COPD_model <-
  ciliated_COPD %>% 
  mutate(mdl = map(data,
                 ~glm(group ~ Counts,
                      data = .x,
                      family = binomial(link = "logit"))),
         tidy = map(mdl,
                    tidy,
                    conf.int = TRUE)) %>% 
  unnest(tidy) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(adj_p.value=
           p.adjust(p.value,
                    method="bonferroni")) %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",
                                   TRUE ~ "unsignificant")) %>% 
  mutate(adj_identified_as = case_when(adj_p.value < 0.05 ~ "significant",
                                       TRUE ~ "unsignificant"))

# wihout correction 2,042 "significant"
# with bonferroni/holm/fdr correction: 0

IPF_model <-
  ciliated_IPF %>% 
  mutate(mdl = map(data,
                 ~glm(group ~ Counts,
                      data = .x,
                      family = binomial(link = "logit"))),
         tidy = map(mdl,
                    tidy,
                    conf.int = TRUE)) %>%
  unnest(tidy) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(adj_p.value = p.adjust(p.value,
                                method = "bonferroni")) %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",
                                   TRUE ~ "unsignificant")) %>% 
  mutate(adj_identified_as = case_when(adj_p.value < 0.05 ~ "significant",
                                       TRUE ~ "unsignificant"))

COPD_against_IPF_model <-
  COPD_against_IPF %>% 
  mutate(mdl = map(data,
                 ~glm(group ~ Counts,
                      data = .x,
                      family = binomial(link = "logit"))),
         tidy = map(mdl,
                    tidy,
                    conf.int = TRUE)) %>%
  unnest(tidy) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(adj_p.value = p.adjust(p.value,
                                method = "bonferroni")) %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",
                                   TRUE ~ "unsignificant")) %>% 
  mutate(adj_identified_as = case_when(adj_p.value < 0.05 ~ "significant",
                                       TRUE ~ "unsignificant"))


multinomial_model <-  
  multinomial %>% 
  mutate(mdl = map(data,
                 ~multinom(group ~ Counts,
                           data = .x)),
         tidy = map(mdl,
                    tidy,
                    conf.int = TRUE)) %>%
  unnest(tidy) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(adj_p.value = p.adjust(p.value,
                                method="bonferroni")) %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",
                                   TRUE ~ "unsignificant")) %>% 
  mutate(adj_identified_as = case_when(adj_p.value < 0.05 ~ "significant",
                                       TRUE ~ "unsignificant"))

sig_without_correction <-
  c(
    COPD_model %>% 
      filter(identified_as == "significant") %>% 
      nrow(),
    IPF_model %>% 
      filter(identified_as == "significant") %>% 
      nrow(),
    COPD_against_IPF_model %>% 
      filter(identified_as == "significant") %>% 
      nrow(),
    multinomial_model %>% 
      filter(identified_as == "significant") %>% 
      nrow()
  )

sig_with_bonferroni_correction <-
  c(
    COPD_model %>% 
      filter(adj_identified_as == "significant") %>% 
      nrow(),
    IPF_model %>% 
      filter(adj_identified_as == "significant") %>% 
      nrow(),
    COPD_against_IPF_model %>% 
      filter(adj_identified_as == "significant") %>% 
      nrow(),
    multinomial_model %>% 
      filter(adj_identified_as == "significant") %>% 
      nrow()
  )

significant_genes_without_correction <-
  tibble(sig_without_correction) %>% 
  ggplot() +
  geom_bar(aes(
    c("COPD_model",
      "IPF_model",
      "COPD_against_IPF",
      "multinomial_model"),
    sig_without_correction),
    stat="identity",
    fill = "lightblue") +
  xlab("Model") +
  ylab("Count") +
  theme_gray() +
  ggtitle("With Boniferri Correction") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2200))
  
significant_genes_with_correction <-
  tibble(sig_with_bonferroni_correction) %>% 
  ggplot() +
  geom_bar(aes(
    c("COPD_model",
      "IPF_model",
      "COPD_against_IPF",
      "multinomial_model"),
    sig_with_bonferroni_correction), 
    stat="identity",
    fill = "red") +
  xlab("Model") +
  ylab("Count") +
  theme_gray() +
  ggtitle("Without Boniferri Correction") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2200))
  
the_power_of_boniferri <- 
  significant_genes_without_correction + 
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  significant_genes_with_correction +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  plot_annotation(
    title = "Genes Identified As Significant"
    # ,subtitle = 
    # ,caption = 
    )
# Save Plots --------------------------------------------------------------

ggsave("Data/04_i_cells_count_plot.png", 
        plot = cells_count_plot)
ggsave("Data/04_i_cells_type_plot.png",
        plot = cells_type_plot)
ggsave("Data/04_i_plot_epithelial_subtypes.png",
        plot = plot_epithelial_subtypes)
ggsave("Data/04_i_the_power_of_boniferri.png",
        plot = the_power_of_boniferri)


# Remove Data -------------------------------------------------------------

rm(ciliated,
   cells_count_plot,
   cells_type_plot,
   plot_epithelial_subtypes,
   the_power_of_boniferri,
   ciliated_COPD,
   ciliated_IPF,
   COPD_against_IPF,
   multinomial,
   COPD_model,
   IPF_model,
   COPD_against_IPF_model,
   multinomial_model,
   sig_without_correction,
   sig_with_bonferroni_correction,
   significant_genes_without_correction,
   significant_genes_with_correction
)
