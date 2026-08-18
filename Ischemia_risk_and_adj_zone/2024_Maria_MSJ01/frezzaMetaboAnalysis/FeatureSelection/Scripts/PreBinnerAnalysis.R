
#' @param ProjectFolder String with path to a project folder. (ie. "Z:\1. Mass Spec Projects\RUGARLI\Rrejusha PARAYIL\2023_Rrejusha_RP01")
#' @param FolderName The output folder name inside the 5. XXXX Analysed Data folder, as this is the analysis done before Binner the name is standardized to preBinnerAnalysis
#' @param CDInputFolder  The standard CD Analysis output should be in the folder '3. XXXX CD Data/XXXX'. In case we have another CD Analysis in the folder (ie. XXXX_test) and would like to use it as input we set the CDInputFolder to the additional name (ie. CDInputFolder = "test") 


PreBinnerAnalysis <- function(Folder_paths,
                              scriptpath = scriptpath,
                              CDInputFolder = NULL, #  CDInputFolder = "test"
                              FolderName = NULL,
                              InputData = RawDataPrepared,
                              RemoveOutlierSamples = TRUE,
                              SampleVector = NULL,
                              
                              Input_SettingsInfo = c(Conditions = "Conditions",PoolSamples = "Pool"),
                              OutputName = NULL,
                              ForceRun = FALSE
){
  
  
  Raw_Data_Annotated_SplitRemoved <- InputData
  Raw_Data_Annotated_SplitRemoved <- column_to_rownames(Raw_Data_Annotated_SplitRemoved, "ID")
  
####-------------------- Deal with adducts aka Binner ---------------------------####
## Select only 1 from all adducts for a feature. Preferably M+-H 
# Binner section. TO run Binner we have to preprocess the data

pos_binner_data <- Raw_Data_Annotated_SplitRemoved[-grep("]-", Raw_Data_Annotated_SplitRemoved$`Reference Ion`),] %>% select(-name, -`Reference Ion`)
pos_binner_data$Compound <- rownames(pos_binner_data)#paste0(pos_binner_data$mz, "_", pos_binner_data$RT)
neg_binner_data <- Raw_Data_Annotated_SplitRemoved[grep("]-", Raw_Data_Annotated_SplitRemoved$`Reference Ion`),] %>% select(-name, -`Reference Ion`)
neg_binner_data$Compound<- rownames(neg_binner_data)#paste0(neg_binner_data$mz, "_", neg_binner_data$RT)
data_list <- list(data_negative = neg_binner_data, data_positive = pos_binner_data)

####---------------------- Load Metadata -------------------------####
# usethis::edit_file(file.path(scriptpath,"LoadMetadata.R"))
source(file.path(scriptpath,"LoadMetadata.R"))
LoadMetadataRes <- LoadMetadata(Folder_paths = Folder_paths,
                                scriptpath = scriptpath,
                                data = neg_binner_data) # here we use one dataste. it is used to add the pool samples on the SSF

# Filter the metadata based on the data 
LoadMetadataRes$Metadata <- LoadMetadataRes$Metadata[rownames(LoadMetadataRes$Metadata) %in% colnames(data_list$data_negative),]

#####################################################
####----------- Make Changes to SSF Here --------####
library(dplyr)
library(stringr)
ssf <- LoadMetadataRes$Metadata
ssf <- ssf %>%
  mutate(
    First_Part = str_split_fixed(Conditions, "_", 2)[, 1],
    First_Part = if_else(
      str_detect(First_Part, "\\d$"),
      First_Part,
      str_c(str_split_fixed(Conditions, "_", 3)[, 1], "_", str_split_fixed(Conditions, "_", 3)[, 2])
    ),
    Rest = str_remove(Conditions, str_c(First_Part, "_?"))
  )

ssf <- ssf %>%
  mutate(First_Part = str_sub(First_Part, 1, -2))
ssf$Conditions <- paste0(ssf$First_Part, "_", ssf$Rest)
ssf[grep("pool", tolower(ssf$Conditions)),1] <- "Pool"
LoadMetadataRes$Metadata <- ssf%>% select(1:5)




### ------------s Preprocessing of Pos and Neg for Binner ---------------------####
Processed_Datasets <- list()
# For each mode in the list do the full preprocessing
for (i in seq_along(data_list)) { # i= 1
  dataset <- data_list[[i]]
  OutputName <- names(data_list)[i]
  
  # Process dataset
  dataset <- as.data.frame(dataset)
  rownames(dataset) <- dataset$Compound
  dataset <- dataset %>% select(-Compound, -mz, -RT)
  dataset  <- t(dataset) %>% as.data.frame()
  
  # Check input data
  # usethis::edit_file(file.path(scriptpath,"InputDataCheck.R"))   
  source(file.path(scriptpath,"InputDataCheck.R"))
  InputDataCheck(Input_data=dataset,
                 Input_SettingsFile = LoadMetadataRes$Metadata,
                 Input_SettingsInfo = Input_SettingsInfo)
  
  message("The dimentions of Input Data: ",dim(dataset)[1], " x ", dim(dataset)[2])
  message("The dimentions of Meta Data: ",dim(LoadMetadataRes$Metadata)[1], " x ", dim(LoadMetadataRes$Metadata)[2])
  message("The dimentions of order samples Data: ",length(LoadMetadataRes$RunOrder))
  
  ####----------------- Zero variance check for features --------------####
  # usethis::edit_file(file.path(scriptpath,"ZeroVarianceCheck.R"))   
  source(file.path(scriptpath,"ZeroVarianceCheck.R"))
  ZeroVarCheckRes <- ZeroVarCheck(Result_Folder = Folder_paths$PreBinnerFolder,
                                  Input_data= Input_data, 
                                  OutputName = OutputName,
                                  RemoveFeatures = TRUE)
  
  
  # Pool Estimation
  # usethis::edit_file(file.path(scriptpath,"PoolEstimation.R"))   
  source(file.path(scriptpath,"PoolEstimation.R"))
  Pool_Estimation_res <- PoolEstimation(Result_Folder = Folder_paths$PreBinnerFolder,
                                        RDS_folder = Folder_paths$PreBinnerRDS,
                                        scriptpath = scriptpath,
                                        OutputName = OutputName,
                                        ForceRun = ForceRun,
                                        #COLORS = "Conditions",
                                        Input_data = ZeroVarCheckRes$Input_data_filtered,
                                        Input_SettingsFile = Input_SettingsFile,
                                        Input_SettingsInfo = Input_SettingsInfo,
                                        Sample_Order = LoadMetadataRes$RunOrder,
                                        #common entries
                                        Remove_unstable_features = TRUE,
                                        Threshold_cv = 30,
                                        Save_as_Results = "xlsx", 
                                        Save_as_Plot = "svg")
  
  # Feature filtering
  # usethis::edit_file(file.path(scriptpath,"FeatureFiltering.R"))   
  source(file.path(scriptpath,"FeatureFiltering.R"))
  FeatureFilteringRes <- FeatureFiltering(Result_Folder = Folder_paths$PreBinnerFolder,
                                          RDS_folder = Folder_paths$PreBinnerRDS,
                                          OutputName = OutputName,
                                          ForceRun = ForceRun,
                                          Input_data=Pool_Estimation_res$DF$Filtered_Input_data,
                                          Input_SettingsFile=Input_SettingsFile,
                                          Input_SettingsInfo=Input_SettingsInfo,
                                          Mode="Modified", # Standard or Modified or None
                                          Feature_Filt_Value = 0.8)
  
  
  # Missing value imputation
  # usethis::edit_file(file.path(scriptpath,"MissingValueImputation.R"))   
  source(file.path(scriptpath,"MissingValueImputation.R"))
  MissingValueImputationRes <- MissingValueImputation(Result_Folder = Folder_paths$PreBinnerFolder,
                                                      RDS_folder = Folder_paths$PreBinnerRDS,
                                                      OutputName = OutputName,
                                                      ForceRun = ForceRun,
                                                      Filtered_matrix=FeatureFilteringRes$Filtered_matrix,
                                                      Input_SettingsFile=Input_SettingsFile,
                                                      Input_SettingsInfo=Input_SettingsInfo,
                                                      Algorithm = "Mean") # or "MAI"
  
  
  # Maybe add log2 normalisation here? 
  
  # Normalisation
  # usethis::edit_file(file.path(scriptpath,"Normalisation.R"))   
  source(file.path(scriptpath,"Normalisation.R"))
  NormalisationRes <- Normalisation(Result_Folder = Folder_paths$PreBinnerFolder,
                                    RDS_folder = Folder_paths$PreBinnerRDS,
                                    scriptpath = scriptpath,
                                    OutputName = OutputName,
                                    ForceRun = ForceRun,
                                    NA_removed_matrix = MissingValueImputationRes$NA_removed_matrix,
                                    Input_SettingsFile = Input_SettingsFile,
                                    Input_SettingsInfo = Input_SettingsInfo,
                                    Method = "Log2", # PQN or TIC
                                    PQN_selected_persentage=100) 
  
  
  
  # Outlier detection
  # usethis::edit_file(file.path(scriptpath,"OutlierDetection.R"))   
  source(file.path(scriptpath,"OutlierDetection.R"))
  OutlierDetectionRes <- OutlierDetection(Result_Folder = Folder_paths$PreBinnerFolder,
                                          RDS_folder = Folder_paths$PreBinnerRDS,
                                          scriptpath= scriptpath,
                                          FolderName = FolderName,
                                          OutputName = OutputName,
                                          ForceRun = ForceRun,
                                          data_norm = NormalisationRes$DF$NormalizedData,
                                          Input_SettingsFile = Input_SettingsFile,
                                          Input_SettingsInfo = Input_SettingsInfo,
                                          OutlierLoop = 2,
                                          npcs=NULL,
                                          HotellinsConfidence = 0.99,
                                          Save_as_Plot = "svg")
  
  
  # Make outlier sample vector
  outlier_sample <- rownames( OutlierDetectionRes$DFs$Processed_data)[ OutlierDetectionRes$DFs$Processed_data$Outliers != "no"]
  
  # Remove outlier samples
  RemoveOutlierSamples = F
  # usethis::edit_file(file.path(scriptpath,"SampleRemoval.R"))   
  source(file.path(scriptpath,"SampleRemoval.R"))
  SampleRemovalRes <- SampleRemoval(Input_data = OutlierDetectionRes$DFs$Processed_data, 
                                    RemoveSamples = RemoveOutlierSamples, 
                                    SampleVector = SampleVector) # or outlier_sample[1:3] for selected.
  
  processedData <- SampleRemovalRes
  # Save only the processed data
  output <-  processedData[,-c(1:grep("Outliers", colnames(processedData)))]
  output <- output %>% t() %>% as.data.frame()
  
  output <- output %>%  rownames_to_column(var = "townames")%>% separate(townames, into = c("mz", "RT"), sep = "_")
  output[1:2] <- lapply(output[1:2], as.numeric)
  output$Compound <- paste0(output$mz,"_", output$RT)
  
  write.table(output, row.names = FALSE, file =  file.path(Folder_paths$PreBinnerFolder,paste0("Compounds_",Folder_paths$ProjectCode,"_",OutputName, "_processed.csv")), sep = ",") # save zero var metabolite list
  
  # ####-------------- BinnerR ----------####
  # # usethis::edit_file(file.path(scriptpath,"Functions/BinnerUsed.R"))
  # source(file.path(scriptpath,"Functions/BinnerUsed.R"))
  # 
  # BinneR_Output <- BinneR(Result_Folder = Folder_paths$BinnerFolder,
  #                         RDS_folder = Folder_paths$RDSFolder,
  #                         ForceRun = ForceRun,
  #                         data=output, 
  #                         OutputName=OutputName,
  #                         CorrelationThreshold=0.9)
  # 
  # special_char_rows <- BinneR_Output[grepl("[+-]", BinneR_Output$Annotation), ]
  # special_char_rowsCompounds <- special_char_rows %>% select(Compound)
  # 
  # output$Compound %in% special_char_rowsCompounds$Compound
  # output <- output %>% filter(!output$Compound %in% special_char_rowsCompounds$Compound)
  
  Processed_Datasets[[OutputName]] <- output
  
  
}

# ####----------------- Merge datasets ----------------####
# message(paste("Merging positive and negative datasets")) # message back the loaded file
# #message("Output location: ", Folder_paths$RawDataFolder)
# # Merge positive and negative datasets
# # data <- merge(Processed_Datasets[[1]],Processed_Datasets[[2]], by=0)
# # colnames(data)[1] <- "Sample"
# # Save the cleaned data
# writexl::write_xlsx(data,file.path(ProjectFolder,paste("5. ", ProjectCode , " Analysed Data","/Results_preBinnerAnalysis/Data_clean_",ProjectCode,".xlsx" , sep = "")))#,showNA = TRUE)



##3. Select only one mode for each metabolite


# # combine datasets
# data <- merge(prepared_data[[1]],prepared_data[[2]], by=0)
# colnames(data)[1] <- "Sample"
# # Save the cleaned data
# writexl::write_xlsx(data,file.path(ProjectFolder,paste("5. ", ProjectCode , " Analysed Data","/Results_preBinnerAnalysis/Data_clean_",ProjectCode,".xlsx" , sep = "")))#,showNA = TRUE)


} # end

