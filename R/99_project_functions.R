# Define project functions ------------------------------------------------
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

load_umicounts <- function(patient_name){
  read_csv(paste0("Data/01_data_",patient_name,".csv.gz"))
}

remove_zero_rows <-function(patient){
  patient <-
    patient %>% 
    mutate(row_sum=
             rowSums(patient),
                           .before=AAACCTGAGCAGCCTC)
  patient <-
    patient %>% 
    filter(total==0)
  
  return(patient)
}