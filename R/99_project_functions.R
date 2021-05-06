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


#Filters out the cells that have fewer than 1000 transcripts recorded
trancript_filter <-function(patient){
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



#Summarizes the mitochondrial transcripts of the cell
mitocompute <- function(patient){
  
  
    summarize(
      
      #mt_count = length(value[week<=2]),
      mt_sum = sum(MT-XXXX[str_match(patient, "MT-")])
  
    )
  
  "MT-"
  
  
  return(patient)
}


#Filters out the cells whose mitochondrial transcript represent >20% of the total
mito_filter  <-function(patient){
  patient %>% 
    mutate(row_sum=
             rowSums(patient))
  patient %>% 
    mutate(mito_sum=
             mitocompute(patient))
  
  patient <-
    patient %>%
    filter(mito_sum/row_sum<0.2)
  

  return(patient)
}


#Slices a patient and return the first x cells
patient_slicer <- function(data,up,down){
  
  patient <- filter(data, Patient_ID =="1I" )%>%
    slice(patient, up:down)
  return(patient)
}


#Slices the metadata of a patient and return the first x cells identities
meta_slicer <-function(meta, id){
  
  m_slice <-filter(meta, CellBarcode_Identity =="id")
  m_slice<- slice(m_slice, up:down)
  m_slice <-m_slice%>%
    select(nUMI,nGene,CellType_Category,Subclass_Cell_Identity)
  
  return(m_slice)
  
}


#Combining the patient's cells with their metadata
combiner <- function(patient,m_slice){
  
  augmented_patient <- patient%>%
    mutate(m_slice)
  
}
#Binding the patients back into an augmented data set
binder <- function (p1,p2,p3,p4,p5,p6){
  
  augmented_data <- bind_rows(
    patient_1,
    patient_2,
    patient_3,
    patient_4,
    patient_5,
    patient_6)
  
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
