library(shiny)
library(shinythemes)
library(tidyverse)
library(ggdark)


load("Data/04_ii_COPD_model")
load("Data/04_ii_IPF_model")


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
