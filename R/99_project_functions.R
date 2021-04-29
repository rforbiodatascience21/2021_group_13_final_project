# Define project functions ------------------------------------------------

### loading functions ###-------------------------------------------------------
make_list_of_patient_data <- function(list_which_patients){
  data <- list()
  for (patient in list_which_patients){
    data <- c(data,list(read_rds(file = paste0("Data/_raw/",patient,".dgecounts.rds"))))
  }
  
  return(data)
}

get_exon_umicounts <-function(patient){
 
  umicount_exon_reads <- 
      as.tibble(patient) %>% 
      pluck("umicount") %>% 
      pluck("exon") %>% 
      pluck("all") %>% 
      as_tibble()
  
  return(umicount_exon_reads)
}

write_umicounts <- function(patient_data,patient_name){
  write_csv(patient_data,
            path = paste0("Data/01_data_",patient_name,".csv.gz"))
}

### cleaning functions ###-------------------------------------------------------

load_umicounts <- function(patient_name){
  read_csv(paste0("Data/01_data_",patient_name,".csv.gz"))
}

remove_zero_rows <-function(patient){
  patient <-
    patient %>% 
    mutate(row_sum=
             rowSums(patient))
    # this is apparently tidyverse but doesnt work
    #https://dplyr.tidyverse.org/articles/rowwise.html
    #rowwise()
    #mutate(row_sum=
    #         sum(c_across())
    #)
  patient <-
    patient %>% 
    filter(row_sum != 0)
  
  return(patient)
}

get_range <- function(patient){
  range <-
    patient %>% 
    select(row_sum) %>% 
    arrange(row_sum)
  
  return(range)
}

refine_range <- function(range, lower_cutoff){
  nineteeth <- 
    range %>% 
    pull(1) %>% 
    quantile(0.9)
  range <-
    range %>% 
    filter(row_sum > lower_cutoff & row_sum < nineteeth)
  
  return(range)
}

refine_data <- function(data,ranges_refined){
  data <-
    inner_join(
      data,ranges_refined,
    )
  
  return(data)
}