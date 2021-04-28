# Define project functions ------------------------------------------------
foo <- function(x){
  return(2*x)
}
bar <- function(x){
  return(x^2)
}

# Function for loading patients
load_patients <-function(patientfile){
  patientfile[21:25] <- read_rds(file = "Data/_raw/"+patientfile)
  umicount_exon_ <- patientfile[21:25] %>% 
    pluck("umicount") %>% 
    pluck("exon") %>% 
    pluck("all") %>% 
    as_tibble()
}