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


gois <-
  read_csv("Data/_raw/disease_genes.csv")


#Wrangling


ciliated <-data %>% 
  filter(Subclass_Cell_Identity == "Ciliated")%>% 
  pivot_longer(cols = TSPAN6:ncol(data),
               names_to = "gene",
               values_to = "Counts")

ciliated_gois<-right_join(ciliated, gois, by=c("gene")) 

ciliated_COPD<- ciliated_gois%>% 
  select(group,
         gene,
         Counts) %>%
  filter(!is.na(Counts)) %>% 
  mutate(Counts = normalise(Counts),
         group = case_when(group == "Control" ~ 0,
                           group == "COPD" ~ 1)) %>% 
  filter(group != "IPF") %>% 
  group_by(gene) %>% 
  nest() %>% 
  ungroup()


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




dendro <- COPD_model%>%
  select(gene,p.value)%>%
  filter(!is.na(p.value))
dendro<-column_to_rownames(dendro, var = "gene")

d<-dist(sqrt(dendro))       
 
dend_row <- d %>% hclust(method = "average") %>% as.dendrogram 
dend_row %>% highlight_branches %>% plot


dend_k <- find_k(dend_row)

plot(dend_k)

heatmaply(sqrt(dendro), Colv = NULL, hclust_method = "average", 
          fontsize_row = 8,fontsize_col = 6,
          k_row = NA, margins = c(60,170, 70,40),
          xlab = "p value", ylab = "genes", main = "A smoking good plot on KOL"
) 



heatmaply(sqrt(dendro), Colv = NULL, hclust_method = "average", 
          fontsize_row = 8,fontsize_col = 6,
          k_row = NA, margins = c(60,50, 70,90),
          xlab = "p value", ylab = "genes", main = "A smoking good plot on KOL",
          plot_method = "plotly", row_dend_left = TRUE
) 
