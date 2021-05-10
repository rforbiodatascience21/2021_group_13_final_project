# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("broom")
library("dendextend")
library("heatmapply")

# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data --------------------------------------------------------------------
data <- read_csv("Data/03_data.csv")

gois <- read_csv("Data/_raw/disease_genes.csv")


# Wrangling --------------------------------------------------------------------
ciliated <- data %>% 
  filter(Subclass_Cell_Identity == "Ciliated") %>% 
  pivot_longer(cols = TSPAN6:ncol(data),
               names_to = "gene",
               values_to = "Counts")

ciliated_gois <- right_join(ciliated,
                            gois, 
                            by=c("gene")) 


# Preparing the COPD model vs Control
ciliated_COPD <- ciliated_gois %>% 
  select(group,
         gene,
         Counts) %>%
  filter(!is.na(Counts)) %>% 
  filter(group != "IPF") %>% 
  mutate(Counts = normalise(Counts),
         group = case_when(group == "Control" ~ 0,
                           group == "COPD" ~ 1)) %>% 
  group_by(gene) %>% 
  nest() %>% 
  ungroup()


# Modelling the COPD vs Control 
COPD_model <- ciliated_COPD %>% 
  mutate(mdl = map(data,
                   ~glm(group ~ Counts,
                        data = .x,
                        family = binomial(link = "logit"))),
         tidy = map(mdl,
                    tidy,
                    conf.int = TRUE)) %>% 
  unnest(tidy) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(adj_p.value = 
           p.adjust(p.value,
                    method = "bonferroni")) %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",
                                   TRUE ~ "unsignificant")) %>% 
  mutate(adj_identified_as = case_when(adj_p.value < 0.05 ~ "significant",
                                       TRUE ~ "unsignificant"))


# Preparing the IPF model vs Control -------------------------------------------
ciliated_IPF <- ciliated_gois %>% 
  select(group,
         gene,
         Counts) %>%
  filter(!is.na(Counts)) %>% 
  mutate(Counts = normalise(Counts),
         group = case_when(group == "Control" ~ 0,
                           group == "IPF" ~ 1)) %>% 
  filter(group != "COPD") %>% 
  group_by(gene) %>% 
  nest() %>% 
  ungroup()


IPF_model <- ciliated_IPF %>% 
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
                    method="bonferroni")) %>% 
  mutate(identified_as = case_when(p.value < 0.05 ~ "significant",
                                   TRUE ~ "unsignificant")) %>% 
  mutate(adj_identified_as = case_when(adj_p.value < 0.05 ~ "significant",
                                       TRUE ~ "unsignificant"))


# Clustering, dendrograms and heatmaps -----------------------------------------
dendro <- COPD_model %>%
  select(gene,p.value) %>%
  filter(!is.na(p.value))

dendro <- column_to_rownames(dendro,
                           var = "gene")

d <- dist(sqrt(dendro))       
 
dend_row <- d %>% 
  hclust(method = "average") %>% 
  as.dendrogram 

dendro_plot <- dend_row %>% 
  highlight_branches %>% 
  plot


dend_k <- find_k(dend_row)

n_clusters <- plot(dend_k)


# Plotting the relationship of genes by significance ---------------------------
heat_dendro <- heatmaply(sqrt(dendro),
                         Colv = NULL, 
                         hclust_method = "average", 
                         fontsize_row = 8,fontsize_col = 6,
                         k_row = NA, margins = c(60,
                                                 50,
                                                 70,
                                                 90),
                         xlab = "p value", 
                         ylab = "genes", 
                         main = "Clustering of gene expression by significance",
                         plot_method = "plotly", 
                         row_dend_left = TRUE)


# PCA data preparation --- Dead end(sofus)
prefit <- ciliated_gois %>% 
  pivot_wider(names_from = "gene",
              values_from = "Counts")


prefit_genes <- prefit %>% 
  select(-mito_sum,
         -nGene,
         -nUMI) %>% 
  select(where(is.numeric))


prefit_genes %>% 
  var(na.rm = FALSE)


prefit_genes %>% 
  select(CP) %>% 
  var(na.rm = TRUE)


pca_fit <- prefit %>% 
  select(CP:ncol(prefit)) %>% # retain only numeric columns
  scale() %>% # scale data
  prcomp() # do PCA


# Save plots -------------------------------------------------------------------



# Removing data ----------------------------------------------------------------
rm(data,
   gois,
   ciliated,
   ciliated_gois,
   ciliated_COPD,
   COPD_model,
   ciliated_IPF,
   IPF_model,
   dendro,
   d,
   dend_row,
   dendro_plot,
   dend_k,
   n_clusters,
   heat_dendro,
   prefit,
   prefit_genes,
   pca_fit)