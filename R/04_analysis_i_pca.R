# Clear workspace ---------------------------------------------------------
#rm(list = ls())


# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("patchwork")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
gravier_data_aug <- read_tsv("data/03_gravier_data_aug.tsv.gz")


# Model Data --------------------------------------------------------------


#This code is a bit messy and parts are probably redundant

model1 <- gravier_data_aug %>% glm(outcome~g2E09,data=.,family=binomial(link="logit"))

gravier_data_long <- 
  gravier_data_aug %>% 
  pivot_longer(-outcome,names_to="gene",values_to = "log2_expr_level") %>% 
  mutate(log2_expr_level=round(log2_expr_level,5))

gravier_data_nested <- gravier_data_long %>% group_by(gene) %>% nest() %>% ungroup()

set.seed(42)
selected <- sample_n(gravier_data_nested,100)
gravier_data_nested_long <- 
  selected %>% 
  mutate(mdl=map(data,~glm(outcome~log2_expr_level,data=.x,family=binomial(link="logit"))))
rm(gravier_data_nested,selected)
gravier_data_nested_long <-
  gravier_data_nested_long %>% 
  mutate(tidy=map(mdl,tidy,conf.int=TRUE)) %>% 
  unnest(tidy) %>% 
  filter(term!="(Intercept)") %>% 
  mutate(identified_as=case_when(p.value<0.05~"significant",TRUE~"unsignificant"),neg_log_p=-log10(p.value))


#select outcome and gene expression levels for each patient

gravier_data_wide = gravier_data_aug %>%
  select(outcome, pull(gravier_data_nested_long, gene))

#do a PCA
pca_fit <- 
  gravier_data_wide %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale=TRUE)

#mutate to get better readable labels

gravier_data_wide <-
  gravier_data_wide %>% 
  mutate(outcome=as.factor(outcome)) %>% 
  mutate(outcome=case_when(outcome=="0"~"good",TRUE~"poor"))

rm(list=setdiff(ls(), c("gravier_data_wide","pca_fit")))


# Plot Data ---------------------------------------------------------------


plot_PCA1_vs_PCA2 <-
  pca_fit %>% 
  augment(gravier_data_wide) %>% 
  ggplot(aes(.fittedPC1,.fittedPC2,colour=outcome))+
  geom_point(size=1.5)+
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(filename = "results/04_plot_PCA1_vs_PCS2.png", plot=plot_PCA1_vs_PCA2, width = 16, height = 9, dpi = 72)


#tidy data to plot rotation matrix

rotation_matrix_data <-
  pca_fit %>% 
  tidy(matrix="rotation") %>% 
  pivot_wider(names_from = "PC",values_from = "value",names_prefix = "PC")

#arrow style for plot

arrow_style <- arrow(
  angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
)
#plot rotation matrix

plot_rotation_matrix <-
  rotation_matrix_data %>% 
  ggplot(aes(PC1,PC2)) +
  geom_segment(xend=0,yend=0,arrow = arrow_style)  +
  geom_text(
    aes(label = column),
    hjust = 1, nudge_x = -0.02, 
    color = "orchid3"
  ) +
  xlim(-0.4,0.2)

ggsave(filename = "results/04_plot_rotation_matrix.png", plot=plot_rotation_matrix, width = 16, height = 9, dpi = 72)


#plot of share of each PC for explaining variance

plot_PCA_shares <-
  pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "lightsalmon2", alpha = 0.8) +
  scale_x_continuous(labels=c(1,5,10,20,50,100),breaks=c(1,5,10,20,50,100))+
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_bw()

ggsave(filename = "results/04_plot_PCA_shares.png", plot=plot_PCA_shares, width = 16, height = 9, dpi = 72)

# cluster observations in PCA1 and PCA2

PCA1_2_data <-
  pca_fit %>% 
  augment(gravier_data_wide)

plot_PCA1_2_clusters <-
  pca_fit %>% 
  augment(gravier_data_wide) %>% 
  select(-.rownames,-outcome) %>% 
  kmeans(centers=2) %>% 
  augment(gravier_data_wide) %>% 
  select(outcome,.cluster) %>% 
  bind_cols(PCA1_2_data %>% select(.fittedPC1,.fittedPC2)) %>% 
  ggplot(aes(.fittedPC1,.fittedPC2,colour=.cluster))+
  geom_point(size=1.5) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  labs(title="Plot of k-means clustering vs. outcome for first two PCAs")

plot_PCA1_2_clusters <- plot_PCA1_2_clusters + plot_PCA1_vs_PCA2

ggsave(filename = "results/04_plot_kmeans_PCAs.png", plot=plot_PCA1_2_clusters, width = 16, height = 9, dpi = 72)


# Remove Data -------------------------------------------------------------

rm(PCA1_2_data,plot_PCA1_2_clusters,arrow_style,gravier_data_wide,pca_fit,plot_PCA_shares,plot_PCA1_vs_PCA2,plot_rotation_matrix,rotation_matrix_data)

