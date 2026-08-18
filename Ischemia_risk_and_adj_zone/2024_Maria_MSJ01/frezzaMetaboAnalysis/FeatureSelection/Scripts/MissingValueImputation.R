
MissingValueImputation <- function(Result_Folder,
                                   RDS_folder,
                                   OutputName = NULL,
                                   ForceRun = TRUE,
                                   CoRe= FALSE,
                                   
                                   
                                   Filtered_matrix,
                                   Input_SettingsFile,
                                   Input_SettingsInfo,
                                   Algorithm = "Mean" # or "MAI"
                                   ){
  

    
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)) {OutputName = ""} else{OutputName = paste("_",OutputName, sep="")}
  if (file.exists(file.path(RDS_folder, paste("MissingValueImputationRes", OutputName,".Rdata", sep = ""))) == TRUE & ForceRun == FALSE) { # Check if RDS exists and forceRun is false then skip script and load result from RDS, otherwise run it
    message("Loading already existing ", paste("MissingValueImputationRes", OutputName,".Rdata", sep = ""))
    load(file.path(RDS_folder, paste("MissingValueImputationRes", OutputName,".Rdata", sep = "")), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running Missing Value Imputation script main body.")
    message("Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `mean of values`.")
    
      MVI_data <- Filtered_matrix 
      MVI_SettingsFile <-  Input_SettingsFile 
      MVI_SettingsFile <- MVI_SettingsFile[match(rownames(MVI_data), rownames(MVI_SettingsFile)), ]
      halfmin <- min(MVI_data, na.rm = T)/2
    
      
      if(CoRe==TRUE){
        for(mediaCond in  unique(Input_SettingsFile[ grepl(Input_SettingsInfo["CoRe_media"],Input_SettingsFile$Conditions),"Conditions"])){ # mediaCond = unique(Input_SettingsFile[ grepl("BasalMedia",Input_SettingsFile$Conditions),"Conditions"])[1]
          replaceNAdf <- MVI_data[ grepl(mediaCond,Input_SettingsFile$Conditions),]

          # find metabolites with NA
          na_percentage <- colMeans(is.na(replaceNAdf)) * 100
          highNA_metabs <- na_percentage[na_percentage>20]
          
          # report metabolites with NA
          if(sum(na_percentage)>0){
            message("NA values were found in Control_media samples for metabolites.")
            # if(sum(na_percentage>20)>0){
            #   message("Metabolites with high NA load in Control_media samples are: ",paste(names(highNA_metabs), collapse = ", "), ".")
            # }
          }
          # if all values are NA set to 0
          replaceNAdf[,which(sapply(replaceNAdf, function(x)all(is.na(x))))]=0
          # If there is at least 1 value use the half minimum per feature
          replaceNAdf <- apply(replaceNAdf, 2,  function(x) {x[is.na(x)] <-  min(x, na.rm = TRUE)/2
          return(x)
          }) %>% as.data.frame()
          
          # replace the samples in the original dataframe
          MVI_data[rownames(MVI_data) %in% rownames(replaceNAdf), ] <- replaceNAdf
        }
      }
    
      
      
      
      if(Algorithm=="Mean"){
      for (feature  in colnames(MVI_data)){ #feature =  colnames(MVI_data)[17]
        feature_data <- merge(MVI_data[feature] , MVI_SettingsFile %>% select(Input_SettingsInfo[["Conditions"]]), by= 0, sort = FALSE)
        feature_data <- column_to_rownames(feature_data, "Row.names")
        feature_data[feature_data==0] <- NA
        
        split_Input <- split(feature_data, MVI_SettingsFile[,"Conditions"]) # split data frame into a list of data frames by condition
        
        for (m in split_Input){ # Select dataframe  # m <- split_Input[[2]]
          if( sum(is.na(m[,1]))== (length(m[,1])) |  sum(is.na(m[,1]))== (length(m[,1])-1) ){ #  if  all is NA or we have only 1 value then put half min or 0 in whole group 
            MVI_data[rownames(MVI_data) %in% rownames(m),feature] <- halfmin # or 0

          }else{ # if we have more than 1 value
              MVI_data[rownames(MVI_data) %in% rownames(m),feature][is.na(m[,1])]<- mean(m[,1],na.rm = TRUE )
          }
        }
        }   
      } else if( Algorithm=="MAI"){
      # Use mechanism aware imputation
      set.seed(1234)
      Results = MAI::MAI(data_miss = MVI_data, # The data with missing values
                         MCAR_algorithm = "BPCA", # The MCAR algorithm to use
                         MNAR_algorithm = "Single", # The MNAR algorithm to use
                         assay_ix = 1, # If SE, designates the assay to impute
                         forest_list_args = list( # random forest arguments for training
                           ntree = 300,
                           proximity = FALSE),
                         verbose = TRUE # allows console message output
      )
      Results <- Results$Imputed_data %>% as.data.frame()
      rownames(Results) <- rownames(MVI_data)
      colnames(Results) <- colnames(MVI_data)
      MVI_data <- Results
    }
 
      set.seed(123) # Set seed for reproducibility
      MVI_data[MVI_data == 0] <- 1 + runif(sum(MVI_data == 0), min = 0.01, max = 0.1)
      
      
      
    # Make Result list
    MissingValueImputationRes = list("NA_removed_matrix"= MVI_data)
    
    # Save output RDS
    save(MissingValueImputationRes, file = file.path(RDS_folder, paste("MissingValueImputationRes",OutputName,".Rdata", sep = "")))
  }
  return(MissingValueImputationRes)
}
  
