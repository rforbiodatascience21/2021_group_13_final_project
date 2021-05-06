# Define project functions ------------------------------------------------

library(tidyverse)

### loading functions ###-------------------------------------------------------


get_exon_umicounts <-function(patient){
 
  umicount_exon_reads <- 
      as.tibble(patient) %>% 
      pluck("umicount") %>% 
      pluck("exon") %>% 
      pluck("all") %>% 
      as_tibble()
  
  umicount_exon_reads <- 
    sample_n(data,100)
  
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
    filter(row_sum > 1000)
  
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


# this is a very slow function
normalise <- function(patient){
  first_row <-
    patient %>% 
    slice(1) %>% 
    as.numeric()
  first_row_sum <-
    first_row %>% 
    pluck(ncol(patient))
  result <- first_row/first_row_sum
  
  
  seq = 1:nrow(patient)
  for (row in seq){
    cur_row <-
      patient %>% 
      slice(row) %>% 
      as.numeric()
    cur_row_sum <-
      cur_row %>% 
      pluck(ncol(patient))
    cur_row <- cur_row/cur_row_sum
    result <- rbind(result,cur_row)
  }
    return(as.tibble(result))
}
