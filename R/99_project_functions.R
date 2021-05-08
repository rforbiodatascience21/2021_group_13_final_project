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
transcript_filter <-function(data){
  
  data <-
    data %>% 
    filter(nUMI > 1000)
  
  return(data)
}


#Filters out the cells whose mitochondrial transcript represent >20% of the total
mito_filter  <-function(data,id,up,down){
  
  patient<- patient_slicer(data,id,up,down)
  
  mt_selection <-select(patient, starts_with("MT"))
  
  mt_selection<-mt_selection%>%
    discard(~all(is.na(.x))) %>%
    map_df(~.x)
  
  mt_sum<-mt_selection%>% 
    mutate(mito_sum=rowSums(mt_selection))

  patient <- patient %>%
    mutate(
      select(
        mt_sum,mito_sum))
  
  patient <- patient %>%
    filter(mito_sum/nUMI<0.2)
  

  return(patient)
}


#Slices a patient and return the first x cells
patient_slicer <- function(data,id,up,down){
  
  patient <- filter(data, Patient_ID ==id )%>%
    slice(patient, up:down)
  return(patient)
}


#Slices the metadata of a patient and return the first x cells identities
meta_slicer <-function(meta, id, up, down){
  
  m_slice <-filter(meta, Subject_Identity ==id)
  m_slice<- slice(m_slice, up:down)
  m_slice <-m_slice%>%
    select(CellBarcode_Identity,nUMI,nGene,CellType_Category,Subclass_Cell_Identity)
  
  return(m_slice)
  
}


#Matching the cells in the metadata file to the remaining cells in the patient data sets after the filtering
meta_matcher <-function(meta,patient,id){
  
  m_slice <- meta_slicer(meta,id)
  patient_barcodes<-select(patient,Cell_Barcode)
  meta_barcodes<-filter(metadata,Subject_Identity==id)
  meta_barcodes<-select(barcodes,Subject_Identity)
  #From that point on it doesn't work yet
  patient_barcodes<-meta_barcodes%>%mutate(Subject_Identity)
  patient_barcodes <- within(patient_barcodes, str_c(Subject_Identity,Cell_Barcode))
  patient <- patient_barcodes%>%mutate(main_column)
  patient<-patient%>%left_join(m_slice, by="CellBarcode_Identity")
}


#Combining the patient's cells with their metadata
#Change this to drop the second CellBarcode_Identity column that is appended in meta_slicer
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

normalise <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
