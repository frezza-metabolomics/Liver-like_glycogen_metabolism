


SampleRemoval <- function(Input_data, RemoveSamples = FALSE, SampleVector = NULL){
  
  if (RemoveSamples == T){
    if(is.null(SampleVector) == T){
      SampleVector <- rownames( Input_data)[ Input_data$Outliers != "no"]
    }
    if(length(SampleVector>0)){
      clean_dataset <- Input_data %>% filter(!rownames(Input_data) %in% SampleVector)
      message("Removing Outlier Samples: ", paste0(SampleVector, by=", "))
    }
  }else{
    message("Outlier Sample Removal == FALSE")
    clean_dataset <- Input_data
  }
  
  return(clean_dataset)
}