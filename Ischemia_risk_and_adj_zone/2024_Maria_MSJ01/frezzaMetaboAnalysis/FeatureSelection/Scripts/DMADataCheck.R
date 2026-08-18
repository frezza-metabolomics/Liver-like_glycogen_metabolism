

DMADataCheck <- function(RDS_folder,
                         DMA_data,
                         DMA_SettingsFile,
                         DMA_SettingsInfo,
                         ForceRun = TRUE,                         
                         FolderName = NULL,
                         OutputName = NULL
    ){
  

  
  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)){ OutputName <- "DMADataCheckRes.Rdata"  }else{ OutputName <- paste("DMADataCheck_" ,OutputName, ".Rdata" ,sep = "" )}
  if (file.exists(file.path(RDS_folder,  OutputName))==TRUE & ForceRun == FALSE) {
  message("Loading already existing ", paste(OutputName))
   # load(file.path(RDS_folder,  OutputName), .GlobalEnv)
    
    load(file.path(RDS_folder,  OutputName), envir = environment())
    
    assign("DMA_data", DMADataCheckRes$DMA_data, envir = .GlobalEnv)
    assign("DMA_SettingsFile", DMADataCheckRes$DMA_SettingsFile, envir = .GlobalEnv)
    assign("DMA_SettingsInfo", DMADataCheckRes$DMA_SettingsInfo, envir = .GlobalEnv)
    
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    
    
  ################################################################################################################################################################################################
  ## ------------ Check Input files ----------- ##
  #1. DMA_data and Conditions
  if(class(DMA_data) != "data.frame"){
    stop("DMA_data should be a data.frame. It's currently a ", paste(class(DMA_data), ".",sep = ""))
  }
  if(any(duplicated(row.names(DMA_data)))==TRUE){
    stop("Duplicated row.names of DMA_data, whilst row.names must be unique")
  } else{
    Test_num <- apply(DMA_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("DMA_data needs to be of class numeric")
    } else{
      Test_match <- merge(DMA_SettingsFile, DMA_data, by.x = "row.names",by.y = "row.names", all =  FALSE, sort = FALSE)
      if(nrow(Test_match) ==  0){
        stop("row.names DMA_data need to match row.names DMA_SettingsFile.")
      } else{
        DMA_data <- DMA_data
      }
    }
  }

  
  ## ------------ Check Input SettingsInfo ----------- ##
  #3. DMA_SettingsInfo
  if(DMA_SettingsInfo[["Conditions"]] %in% DMA_SettingsInfo==TRUE){
    if(DMA_SettingsInfo[["Conditions"]] %in% colnames(DMA_SettingsFile)== FALSE){
      stop("The ",DMA_SettingsInfo[["Conditions"]], " column selected as Conditions in DMA_SettingsInfo was not found in DMA_SettingsFile. Please check your input.")
    }else{# if true rename to Conditions
      DMA_SettingsFile<- DMA_SettingsFile%>%
        dplyr::rename("Conditions"= paste(DMA_SettingsInfo[["Conditions"]]) )
    }
  }else{
    stop("You have to provide a DMA_SettingsInfo for conditions.")
  }
  
  ##########################
  if("denominator" %in% names(DMA_SettingsInfo)==TRUE){
    if(DMA_SettingsInfo[["denominator"]] %in% DMA_SettingsFile$Conditions==FALSE){
      stop("The ",DMA_SettingsInfo[["denominator"]], " column selected as denominator in DMA_SettingsInfo was not found in DMA_SettingsFile. Please check your input.")
    }else{
      denominator <- DMA_SettingsInfo[["denominator"]]
    }
  }
  if("numerator" %in% names(DMA_SettingsInfo)==TRUE){
    if(DMA_SettingsInfo[["numerator"]] %in% DMA_SettingsFile$Conditions  == FALSE){
      stop("The ",DMA_SettingsInfo[["numerator"]], " column selected as numerator in DMA_SettingsInfo was not found in DMA_SettingsFile. Please check your input.")
    }else{
      numerator <- DMA_SettingsInfo[["numerator"]]
    }
  }
  if("denominator" %in% names(DMA_SettingsInfo)==FALSE  & "numerator" %in% names(DMA_SettingsInfo) ==TRUE){
    stop("Check input. The selected denominator option is empty while ",paste(DMA_SettingsInfo[["numerator"]])," has been selected as a numerator. Please add a denminator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison." )
  }
  
  ## ------------ Check Denominator/numerator ----------- ##
  #4.  Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
  if("denominator" %in% names(DMA_SettingsInfo)==FALSE  & "numerator" %in% names(DMA_SettingsInfo) ==FALSE){
    # all-vs-all: Generate all pairwise combinations
    conditions = DMA_SettingsFile$Conditions
    denominator <-unique(DMA_SettingsFile$Conditions)
    numerator <-unique(DMA_SettingsFile$Conditions)
    comparisons <- combn(unique(conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("denominator" %in% names(DMA_SettingsInfo)==TRUE  & "numerator" %in% names(DMA_SettingsInfo)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    conditions = DMA_SettingsFile$Conditions
    denominator <- DMA_SettingsInfo[["denominator"]]
    numerator <-unique(DMA_SettingsFile$Conditions)
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }else if("denominator" %in% names(DMA_SettingsInfo)==TRUE  & "numerator" %in% names(DMA_SettingsInfo)==TRUE){
    # one-vs-one: Generate the comparisons
    comparisons <- matrix(c(denominator, numerator))
    #Settings:
    MultipleComparison = FALSE
    all_vs_all = FALSE
  }
  
  #7. Are sample numbers enough?
  Num <- DMA_data %>%
    filter(DMA_SettingsFile$Conditions %in% numerator) %>%
    select_if(is.numeric)#only keep numeric columns with metabolite values
  Denom <- DMA_data %>%
    filter(DMA_SettingsFile$Conditions %in% denominator) %>%
    select_if(is.numeric)
  
  if(nrow(Num)==1){
    stop("There is only one sample available for ", numerator, ", so no statistical test can be performed.")
  } else if(nrow(Denom)==1){
    stop("There is only one sample available for ", denominator, ", so no statistical test can be performed.")
  }else if(nrow(Num)==0){
    stop("There is no sample available for ", numerator, ".")
  }else if(nrow(Denom)==0){
    stop("There is no sample available for ", denominator, ".")
  }
  
  ## ------------ Check Missingness ------------- ##
  #7.
  # If missing value imputation has not been performed the input data will most likely contain NA or 0 values for some metabolites, which will lead to Log2FC = NA.
  # Here we will check how many metabolites this affects in Num and Denom, and weather all replicates of a metabolite are affected.
  Num_Miss <- replace(Num, Num==0, NA)
  Num_Miss <- Num_Miss[, (colSums(is.na(Num_Miss)) > 0), drop = FALSE]
  
  Denom_Miss <- replace(Denom, Denom==0, NA)
  Denom_Miss <- Denom_Miss[, (colSums(is.na(Denom_Miss)) > 0), drop = FALSE]
  
  if((ncol(Num_Miss)>0 & ncol(Denom_Miss)==0)){
    message("In `numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- colnames(Num_Miss)
  } else if(ncol(Num_Miss)==0 & ncol(Denom_Miss)>0){
    message("In `denominator` ",paste0(toString(denominator)), ", NA/0 values exist in ", ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- colnames(Denom_Miss)
  } else if(ncol(Num_Miss)>0 & ncol(Denom_Miss)>0){
    message("In `numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), " and in `denominator`",paste0(toString(denominator)), " ",ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "),". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
  } else{
    message("There are no NA/0 values")
    Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
    Metabolites_Miss <- unique(Metabolites_Miss)
  }
  
  assign("DMA_data", DMA_data, envir = .GlobalEnv)
  assign("DMA_SettingsFile", DMA_SettingsFile, envir = .GlobalEnv)
  assign("DMA_SettingsInfo", DMA_SettingsInfo, envir = .GlobalEnv)
  
  DMADataCheckRes = list("DMA_data" = DMA_data, "DMA_SettingsFile" =  DMA_SettingsFile, "DMA_SettingsInfo" = DMA_SettingsInfo)
  # Save output RDS
  save(DMADataCheckRes, file = file.path(RDS_folder,  OutputName))
  }

}