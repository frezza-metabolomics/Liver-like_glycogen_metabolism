
FeatureFiltering <- function(Result_Folder,
                             RDS_folder,
                             OutputName = NULL,
                             ForceRun = T,
                             
                             Input_data,
                             Input_SettingsFile,
                             Input_SettingsInfo,
                             
                             Mode="Modified", # Modified, Standard or none
                             Feature_Filt_Value = 0.8
                             ){
  
 
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)) {OutputName = ""} else{OutputName = paste("_",OutputName, sep="")}
  if (file.exists(file.path(RDS_folder, paste("FeatureFilteringRes", OutputName,".Rdata", sep = ""))) == TRUE & ForceRun == FALSE) { # Check if RDS exists and forceRun is false then skip script and load result from RDS, otherwise run it
    message("Loading already existing ", paste("FeatureFilteringRes", OutputName,".Rdata", sep = ""))
    load(file.path(RDS_folder, paste("FeatureFilteringRes", OutputName,".Rdata", sep = "")), .GlobalEnv)
  
  ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running Feature filtering script main body.")
    Input_data <- replace(Input_data, Input_data==0, NA)
    
    # ignore pool samples
    if ("PoolSamples" %in% names(Input_SettingsInfo)) {
      feat_filt_data <- Input_data %>% filter(!Input_SettingsFile[,"Conditions"] == Input_SettingsInfo[["PoolSamples"]])
      feat_filt_Conditions <- Input_SettingsFile[!c(Input_SettingsFile[,"Conditions"] == Input_SettingsInfo[["PoolSamples"]]),"Conditions" ]
      }else{
      feat_filt_data <- Input_data
      feat_filt_Conditions <-  Input_SettingsFile$Conditions
    }
    NAMES=0
    if (Mode ==  "Modified"){
    message("Here we apply the modified 80%-filtering rule that takes the class information (Column `Conditions`) into account, which additionally reduces the effect of missing values. REF: Yang et. al., (2015), doi: 10.3389/fmolb.2015.00004)")
    
    # Input checks
    if(is.null(unique(feat_filt_Conditions)) ==  TRUE){
      stop("Condition information is missing from the Experimental design.")
    }
    if(length(unique(feat_filt_Conditions)) ==  1){
      stop("To perform the Modified feature filtering there have to be at least 2 different Conditions in the `Condition` column in the Experimental design. Consider using the Standard feature filtering option.")
    }
    
    miss <- c()
    split_Input <- split(feat_filt_data, feat_filt_Conditions) # split data frame into a list of data frames by condition
    for (m in split_Input){ # Select dataframe  # m <- split_Input[[1]]
      for(i in 1:ncol(m)) { # for each metabolute # i=1
        if(length(which(is.na(m[,i]))) >= (1-Feature_Filt_Value)*nrow(m)) # if we have more NAs than the percentage selected
          miss <- append(miss,i) # add metabolite to the prefiltering list
      }
    }
    
    filt_metabs <- c()
    if(length(miss)>0){ # if there are metabolites in the prefiltering list
      count <- table(miss) 
      for (k in 1:length(count)){ # for each metabolite
        if (count[[k]]== length(split_Input)){ # if the metab is missing for all dataframes 
          filt_metabs <- append(filt_metabs,as.numeric(names(count)[k])) # add metabolite in the filtering list
        }
      }
    }
    
    if(length(filt_metabs) ==  0){ #remove metabolites if any are found
      message("There where no metabolites exluded")
      Input_data <- Input_data
      feat_file_res <- "There where no metabolites exluded"
      write.table(feat_file_res,row.names =  FALSE, file = paste(Result_Folder,"/Filtered_metabolites","_",Feature_Filt_Value*100,"%_",Mode,OutputName,".csv",sep =  ""))
    } else {
      
      NAMES<-colnames(Input_data)[filt_metabs]
      message(length(unique(filt_metabs)) ," metabolites where removed")
      #message(length(unique(filt_metabs)) ," metabolites where removed: ", paste0(names, collapse = ", "))
      
      Input_data <- Input_data[,-filt_metabs]
      write.table(unique(colnames(Input_data)[filt_metabs]),row.names = FALSE, file =  paste(Result_Folder,"/Filtered_metabolites","_",Feature_Filt_Value*100,"%_",Mode,OutputName,".csv",sep =  ""))
    }
  }
    if (Mode ==  "Standard"){
    message("Here we apply the so-called 80%-filtering rule, which removes metabolites with missing values in more than 80% of samples. REF: Smilde et. al. (2005), Anal. Chem. 77, 6729–6736., doi:10.1021/ac051080y")
    message(paste("filtering value selected:", Feature_Filt_Value))
    
    split_Input <- feat_filt_data
    # split_Input[split_Input==0] <- NA
    miss <- c()
    message("***Performing standard feature filtering***")
    for(i in 1:ncol(split_Input)) { # Select metabolites to be filtered for one condition
      if(length(which(is.na(split_Input[,i]))) > (1-Feature_Filt_Value)*nrow(split_Input))
        miss <- append(miss,i)
    }
    
    if(length(miss) ==  0){ #remove metabolites if any are found
      message("There where no metabolites exluded")
      Input_data <- Input_data
      feat_file_res <- "There where no metabolites exluded"
      write.table(feat_file_res,row.names =  FALSE, file = paste(Result_Folder,"/Filtered_metabolites","_",Feature_Filt_Value*100,"%_",Mode,OutputName,".csv",sep =  ""))
    } else {
      NAMES<-unique(colnames(Input_data)[miss])
      message(length(unique(miss)) ," metabolites where removed: ", paste0(NAMES, collapse = ", "))
      Input_data <- Input_data[,-miss]
      write.table(unique(colnames(Input_data)[miss]),row.names =  FALSE, file = paste(Result_Folder,"/Filtered_metabolites","_",Feature_Filt_Value*100,"%_",Mode,OutputName,".csv",sep =  ""))
    }
  }
    if (Mode ==  "None"){
      warning("No feature filtering is selected.")
      Input_data <- as.data.frame(Input_data)
    }
    
    if (length(NAMES) !=0){
    features_filtered <- NAMES %>% as.vector()
    } else {
      features_filtered <- "no metabolites filtered"
    }
    Input_data <- as.data.frame(mutate_all(as.data.frame(Input_data), function(x) as.numeric(as.character(x))))
  
    # Save results
    writexl::write_xlsx(Input_data, file.path(Result_Folder, paste0("Feature_Input_data_filtered_by_",Feature_Filt_Value*100,"%_",Mode,OutputName,".xlsx")))#,showNA = TRUE)
    
    # Make Result list
    FeatureFilteringRes = list("Filtered_matrix"=Input_data, "Features_filtered"=features_filtered)
    
    # Save output RDS
    save(FeatureFilteringRes, file = file.path(RDS_folder, paste("FeatureFilteringRes", OutputName,".Rdata", sep = "")))
    }
  return(FeatureFilteringRes)
  }
                          

