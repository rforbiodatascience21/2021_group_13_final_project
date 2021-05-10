# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("broom")
library("dendextend")
library("heatmapply")
library("ggdendro")
library("patchwork")
library("shiny")
library("shinythemes")

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

# k-means

cluster_data <- ciliated_gois %>% 
  pivot_wider(
    names_from = gene,
    values_from = Counts
  ) %>% 
  select(c("group",CP:AFP)) %>% 
  mutate(across(everything(),
                ~replace_na(.x,0)))

data_to_cluster <-
  cluster_data %>% 
  select(-group)

kclusts <- 
  tibble(k = 1:20) %>%
  mutate(
    kclust = map(k, ~kmeans(data_to_cluster, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, cluster_data)
  )

glanced <- kclusts %>%
  unnest(cols = glanced)

plot_k_clusters <- ggplot(glanced, 
                          aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(breaks=1:20,
                   limits=1:20) +
  ylab("Total Within Sum of Squares") +
  xlab("Number of Clusters") +
  ggplot(glanced,
         aes(k,
             (betweenss/totss))) +
  geom_line() +
  geom_point() +
  scale_x_discrete(breaks=1:20,
                    limits=1:20) +
  ylim(0:1) +
  ylab("Between/Tot Sum of Sq.") +
  xlab("Number of Clusters") +
  plot_annotation(title = 
                    "Evaluation of K-means with Different Number of Centroids")


### shiny ----------------------------------------------------------------------

significant_identification <- function(dataset,p){
  dataset <-
    dataset %>% 
    mutate(identified_as = 
             case_when(p.value<p~"significant",
                       TRUE~"unsignificant"))
}

manhatten_plot <- function(dataset,p){
  dataset %>% 
    mutate(gene = fct_reorder(as.factor(gene),
                              p.value,
                              .desc = TRUE)) %>% 
    ggplot(aes(gene,
               p.value,
               colour = identified_as)) + 
    geom_point(size = 2) + 
    geom_hline(yintercept = p,
               linetype = "dashed") + 
    labs(x="Gene",
         y="p-value") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle=45,
                                     size=3))
}

filter_sig_genes <- function(dataset){
  data <-
    dataset %>% 
    filter(identified_as=="significant") %>% 
    select("gene","p.value") %>% 
    mutate(p.value = as.character(p.value))
  return(data)
  
  
}

count_sig_genes <- function(dataset){
  data <-
    dataset %>% 
    filter(identified_as=="significant") %>% 
    select("gene","p.value") %>% 
    mutate(p.value = as.character(p.value))%>% 
    count()
  return(data)
  
  
}

ui <- fluidPage(
  theme = shinytheme("cyborg"),
  fluidRow(
    selectInput(
      "selected_data",
      label = "Select Model",
      choice = c("IPF_vs_Control","COPD_vs_Control"),
    )
  ),
  fluidRow(
    column(8,
           sliderInput("p","p-value",1e-3,0.1,value=0.05,step=0.01)
    ),
    column(4,
           checkboxInput("bon","Bonferroni Correction")
    )
  ),
  fluidRow(
    plotOutput("plot"),
    h3("Genes Identified as Significant:"),
    tableOutput("sig_genes"),
    
    
  )
)

server <- function(input,output,session){
  
  data <- eventReactive(
    {input$p
      input$bon
      input$selected_data},{
        if (input$selected_data=="COPD_vs_Control"){
          if (input$bon==FALSE){
            significant_identification(COPD_model,input$p)
          }
          else{
            significant_identification(COPD_model,input$p/63)
          }
        }
        else{
          if (input$bon==FALSE){
            significant_identification(IPF_model,input$p)
          }
          else{
            significant_identification(IPF_model,input$p/63)
          }
        }
      }
  )
  
  output$plot <- renderPlot(manhatten_plot(data(),input$p))
  
  output$sig_genes <- renderTable(filter_sig_genes(data()))
  
  output$gene_count <- renderTable(count_sig_genes(data()))
}

shinyApp(ui, server)



### end shiny ------------------------------------------------------------------

# dendrograms

#COPD model

dendroCOPD <- COPD_model %>%
  select(gene,p.value) %>%
  filter(!is.na(p.value))

dendroCOPD <- column_to_rownames(dendroCOPD,
                                 var = "gene")


hcCOPD<- hclust(dist(dendroCOPD), 
                "ave")

COPD_clusters <- dendro_data(hcCOPD,
                             type = "rectangle")

COPD_clusterplot <- ggplot() +
  geom_segment(data = segment(COPD_clusters), 
               aes(x = x,
                   y = y, 
                   xend = xend, 
                   yend = yend)) +
  geom_text(data = label(COPD_clusters), 
            aes(x = x,
                y = y, 
                label = label,
                hjust = 0), 
            size = 3) +
  xlab("genes") +
  ylab("distance") +
  labs(title = "COPD vs Control clustering by significance") +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 
                             0)) + 
  theme_minimal()


#IPF model

dendroIPF <- IPF_model %>%
  select(gene,p.value) %>%
  filter(!is.na(p.value))

dendroIPF <- column_to_rownames(dendroIPF,
                                var = "gene")

hcIPF<- hclust(dist(dendroIPF),
               "ave")

IPF_clusters <- dendro_data(hcIPF, 
                            type = "rectangle")

IPF_clusterplot <- ggplot() +
  geom_segment(data = segment(IPF_clusters), 
               aes(x = x, 
                   y = y, 
                   xend = xend, 
                   yend = yend)) +
  geom_text(data = label(IPF_clusters), 
            aes(x = x,
                y = y, 
                label = label,
                hjust = 0), 
            size = 3) +
  xlab("genes") +
  ylab("distance") +
  labs(title = "IPF vs Control clustering by significance") +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 
                             0)) + 
  theme_minimal()

### Start of none-tidyvers part ###

#A non-tidyverse compatible clustering method - NOT USED

d <- dist(sqrt(dendroIPF))     # change to dendroCOPD if wanted  
 
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

ggsave("Data/04_ii_plot_k_clusters.png",
       plot = plot_k_clusters)
ggsave("Data/04_ii_COPD_clusterplot.png",
       plot = COPD_clusterplot)
ggsave("Data/04_ii_IPF_clusterplot.png",
       plot = IPF_clusterplot)
ggsave("Data/04_ii_dendro_plot.png",
       plot = dendro_plot)
ggsave("Data/04_ii_heat/dendro_plot.png",
       plot = heat_dendro)

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
   pca_fit,
   kclusts,
   plot_k_clusters,
   ui)
