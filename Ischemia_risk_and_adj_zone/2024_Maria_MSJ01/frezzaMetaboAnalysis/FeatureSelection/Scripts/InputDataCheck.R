

InputDataCheck <- function(Input_data,
                           Input_SettingsFile,
                           Input_SettingsInfo
                           ){
  
  ## 1. ------------ Setup and installs ----------- ##
  suppressWarnings(suppressMessages(library(tidyverse))) # load tidyverse package
  
  #2. The Input_settings: Input_SettingInfo and Input_SettingFile
  if(is.vector(Input_SettingsInfo)==TRUE & is.null(Input_SettingsFile)==TRUE){
    stop("You have chosen Plot_SettingsInfo option that requires you to provide a DF Plot_SettingsFile.")
  }
  if(is.null(Input_SettingsInfo)==TRUE & is.null(Input_SettingsFile)==FALSE){
    warning("You have added a Plot_SettingsFile DF but the Plot_SettingsInfo option is empty. If you want to preprocess based some experimental condition please specify it in the Input_SettingsInfo.")
  }
  if (is.null(Input_SettingsInfo) == TRUE & is.null(Input_SettingsFile) == TRUE) {
    message("No Input_Settings have been added.")
  }
  
  if (is.vector(Input_SettingsInfo) == TRUE) {
    if ( "Conditions" %in% names(Input_SettingsInfo)){
      if(Input_SettingsInfo[["Conditions"]] %in% colnames(Input_SettingsFile) == FALSE){
        stop("The ",Input_SettingsInfo[["Conditions"]], " column selected as Conditions in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
      }else{# if true rename to Conditions
        Input_SettingsFile<- Input_SettingsFile %>%
          dplyr::rename("Conditions"= paste(Input_SettingsInfo[["Conditions"]]) )
      }
    }
    
    if ("Biological_Replicates" %in% names(Input_SettingsInfo)) {
      if(Input_SettingsInfo[["Biological_Replicates"]] %in% colnames(Input_SettingsFile) == FALSE){
        stop("The ",Input_SettingsInfo[["Biological_Replicates"]], " column selected as Biological_Replicates in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
      }else{# if true rename to Conditions
        Input_SettingsFile<- Input_SettingsFile %>%
          dplyr::rename("Biological_Replicates"= paste(Input_SettingsInfo[["Biological_Replicates"]]) )
      }
    }
    
    if ("Analytical_Replicates" %in% names(Input_SettingsInfo)) {
      if(Input_SettingsInfo[["Analytical_Replicates"]] %in% colnames(Input_SettingsFile)== FALSE){
        stop("The ",Input_SettingsInfo[["Analytical_Replicates"]], " column selected as Analytical_Replicates in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
      }else{# if true rename to Conditions
        Input_SettingsFile<- Input_SettingsFile %>%
          dplyr::rename("Analytical_Replicates"= paste(Input_SettingsInfo[["Analytical_Replicates"]]) )
      }
    }
    
    if ("PoolSamples" %in% names(Input_SettingsInfo)) {
      if(Input_SettingsInfo[["PoolSamples"]] %in% Input_SettingsFile[,Input_SettingsInfo[["Conditions"]]] == FALSE){
        stop("The ",Input_SettingsInfo[["PoolSamples"]], " column selected as PoolSamples in Input_SettingsInfo was not found in Input_SettingsFile. Please check your input.")
      }else{# if truethe conditions of Pools to Pool
        Input_SettingsFile[Input_SettingsFile[Input_SettingsInfo[["Conditions"]]]=="Pool","Conditions"]<- "Pool"
      }
    }
    
  }
  
  Input_data[] <- lapply(Input_data, as.numeric)
  Input_SettingsInfo[["Conditions"]]<- "Conditions"

  # test if all samples of data exist in SettingsFile
  if (all( rownames(Input_data) %in% rownames(Input_SettingsFile))){
    Input_SettingsFile <-   Input_SettingsFile[ rownames(Input_data), , drop = FALSE]
  }else{
    stop("Not all samples in the Input_data esixt in the SSF.")
  }
  
  
  assign("Input_data", Input_data, envir  = .GlobalEnv)
  assign("Input_SettingsFile", Input_SettingsFile, envir  = .GlobalEnv)
  assign("Input_SettingsInfo", Input_SettingsInfo, envir  = .GlobalEnv)
  
  #return()
}

