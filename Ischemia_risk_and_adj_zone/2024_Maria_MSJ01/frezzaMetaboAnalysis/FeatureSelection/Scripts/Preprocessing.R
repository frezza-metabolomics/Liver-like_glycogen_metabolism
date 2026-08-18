Preprocessing <- function(Folder_paths = Folder_paths,
                          scriptpath = scriptpath,
                          FolderName = NULL,
                          OutputName = NULL,
                          ForceRun = FALSE,
                          Input_SettingsInfo = c(Conditions = "Conditions",
                                                 PoolSamples = "Pool", 
                                                 Biological_Replicates = "Biological Replicate",
                                                 Analytical_Replicates = "Technical Replicate"),
                          RemoveOutlierSamples = FALSE, 
                          SampleVector = NULL,
                          ReplicateSum = FALSE,
                          CoRe=FALSE
){
  

  

  # # usethis::edit_file("frezzaMetaboAnalysis/FeatureSelection/Scripts/Functions/LoadData.R")
  # source("frezzaMetaboAnalysis/FeatureSelection/Scripts/Functions/LoadData.R")
  # LoadDataRes <- LoadData(ProjectFolder = ProjectFolder, FolderName=FolderName,ForceRun = ForceRun, DataSource = NULL)
  
  data_location <- file.path (Folder_paths$ProjectFolder,"Results","PreBinner_folder",paste0("Data_clean_", Folder_paths$ProjectCode,".xlsx" ))
  message(paste("Input Data file:",data_location )) # message back the loaded file
  
  suppressMessages(data <- readxl::read_excel(data_location, sheet = 1, col_names = T) %>% as.data.frame()) # Load the file

  # Get only the intensity values
  data <- cbind(data %>% select("Name","m/z", "RT [min]"), data %>% select(grep("Area",colnames(data)))) 
  data <- data %>% select(-grep("(Max.)",colnames(data)))
  # Rename data col names
  if("Name" %in% colnames(data) & "m/z"%in% colnames(data) & "RT [min]"%in% colnames(data)){
    data <- data %>%
      rename("mz"= "m/z",
             "RT"= "RT [min]" ) 
    
    data$mz <- as.numeric(data$mz)
    data$RT <- as.numeric(data$RT)
  }else{
    "Error columns are missing from data"
  }
  
  if(sum(duplicated(paste0(data$`mz`,"_" ,data$RT)))>0){ # check for duplicated feature names mz_RT
    message("There are ", sum(duplicated(paste0(data$`mz`,"_" ,data$RT))) , " duplicated features in the data.")
  }
  rownames(data) <- make.unique(paste0(data$`mz`,"_" ,data$RT))
  
  data <- data %>% select(!c(Name,`mz`, RT ))  %>% as.data.frame()
  colnames(data) <-sub(paste0(".*", Folder_paths$ProjectCode), Folder_paths$ProjectCode, colnames(data))
  colnames(data) <- sub(".raw.*", "", colnames(data))
  data <- data[,!grepl("blank|blk|blnk", tolower(colnames(data)))] # remove any blank samples

  data <- as.data.frame(as.matrix(mutate_all(as.data.frame(data), function(x) as.numeric(as.character(x))))) # make all values numeric
  Output_data <- data
  
  ####---------------------- Load Metadata -------------------------####
  # usethis::edit_file(file.path(scriptpath,"LoadMetadata.R"))
  source(file.path(scriptpath,"LoadMetadata.R"))
  LoadMetadataRes <- LoadMetadata(Folder_paths = Folder_paths,
                                  scriptpath = scriptpath,
                                  data = Output_data) # here we use one dataset. it is used to add the pool samples on the SSF
  
  rownames(LoadMetadataRes$Metadata) <- gsub("Area: ", "", rownames(LoadMetadataRes$Metadata))
  colnames(Output_data) <- gsub("Area: ", "", colnames(Output_data))
  
  # Filter the metadata based on the data 
  LoadMetadataRes$Metadata <- LoadMetadataRes$Metadata[rownames(LoadMetadataRes$Metadata) %in% colnames(Output_data),]
  # Transpose the data df 
  Output_data <- as.data.frame(t(Output_data))
  
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
  
 ssf <- ssf%>%   mutate(First_Part = str_split_fixed(First_Part, "_", 2)[, 1])
  
  ssf <- ssf %>%
    mutate(First_Part = str_sub(First_Part, 1, -2))
  ssf$Conditions <- paste0(ssf$First_Part, "_", ssf$Rest)
  ssf[grep("pool", tolower(ssf$Conditions)),1] <- "Pool"
  
  # Add batch information
  ssf <- ssf %>% mutate(Batch = case_when(
    grepl("Pool", rownames(.)) ~ 3,
    grepl("_b", rownames(.)) ~ 2,
    TRUE ~ 1))
  
  LoadMetadataRes$Metadata <- ssf%>% select(c(1:5,13))
  
  
  # usethis::edit_file(file.path(scriptpath,"InputDataCheck.R"))
  source(file.path(scriptpath,"InputDataCheck.R"))
  InputDataCheck(Input_data=Output_data,
                 Input_SettingsFile = LoadMetadataRes$Metadata,
                 Input_SettingsInfo = Input_SettingsInfo)
  
  
  
  
  message("The dimentions of Input Data: ",dim(Input_data)[1], " x ", dim(Input_data)[2])
  message("The dimentions of Meta Data: ",dim(Input_SettingsFile)[1], " x ", dim(Input_SettingsFile)[2])
  message("The dimentions of order samples Data: ",length(LoadMetadataRes$RunOrder))

  # # # Get only the HG or the MM
   # Input_data_sel <- Input_data[grepl("HG", rownames(Input_SettingsFile)),]
   # Input_SettingsFile_sel <- Input_SettingsFile[grepl("HG", rownames(Input_SettingsFile)),]
   # OutputName = "_HG"
   Input_data_sel <- Input_data
   Input_SettingsFile_sel <- Input_SettingsFile
  
  ####----------------- Zero variance check for features --------------####
  # usethis::edit_file(file.path(scriptpath,"ZeroVarianceCheck.R"))   
  source(file.path(scriptpath,"ZeroVarianceCheck.R"))
  ZeroVarCheckRes <- ZeroVarCheck(Result_Folder = Folder_paths$PostBinnerFolder,
                                  Input_data= Input_data_sel, 
                                  OutputName = OutputName,
                                  RemoveFeatures = TRUE)
  

  
  # ####----------------- Check batch effects in raw Data --------------------------####
  # usethis::edit_file(file.path(scriptpath,"PCA.R"))
  source(file.path(scriptpath,"PCA.R"))
  pca_data <- ZeroVarCheckRes$Input_data_filtered %>% select(-c(as.vector(  which(is.na(colSums(ZeroVarCheckRes$Input_data_filtered[,-1])))  ) +1) )
  Batch_PCA <- suppressWarnings(PCA(Input_data=pca_data, Plot_SettingsInfo= c(color=  colnames(Input_SettingsFile_sel)[6] ),
                                    Plot_SettingsFile=Input_SettingsFile_sel,
                                    OutputPlotName = "PCA Raw Data Batch",
                                    Save_as_Plot =  NULL))
  ggsave(file=paste(Folder_paths$PostBinnerFolder,"/", "PCA_Batch_raw.svg", sep=""), plot=Batch_PCA, width=20, height=15, unit="cm")
  
  Cond_PCA <- suppressWarnings(PCA(Input_data=pca_data, Plot_SettingsInfo= c(color=  colnames(Input_SettingsFile_sel)[1] ),
                                    Plot_SettingsFile=Input_SettingsFile_sel,
                                    OutputPlotName = "PCA Raw Data Cond",
                                    Save_as_Plot =  NULL))
  ggsave(file=paste(Folder_paths$PostBinnerFolder,"/", "PCA_Cond_raw.svg", sep=""), plot=Cond_PCA, width=20, height=15, unit="cm")

#  # Batch effect correction
#   library(limma)
#   batch <- Input_SettingsFile$Batch
#   y2 <- removeBatchEffect(t(ZeroVarCheckRes$Input_data_filtered), batch)
#   par(mfrow=c(1,2))
#   boxplot(as.data.frame(t(ZeroVarCheckRes$Input_data_filtered)),main="Original")
#   boxplot(as.data.frame(y2),main="Batch corrected")
  
  # Pool Estimation
  # usethis::edit_file(file.path(scriptpath,"PoolEstimation.R"))   
  source(file.path(scriptpath,"PoolEstimation.R"))
  Pool_Estimation_res <- PoolEstimation(Result_Folder = Folder_paths$PostBinnerFolder,
                                        RDS_folder = Folder_paths$PostBinnerRDS,
                                        scriptpath = scriptpath,
                                        OutputName = OutputName,
                                        ForceRun = ForceRun,
                                        #COLORS = "Conditions",
                                        Input_data = ZeroVarCheckRes$Input_data_filtered, #Input_data = as.data.frame(t( y2)),
                                        Input_SettingsFile = Input_SettingsFile,
                                        Input_SettingsInfo = Input_SettingsInfo,
                                        Sample_Order = LoadMetadataRes$RunOrder,
                                        #common entries
                                        Remove_unstable_features = TRUE,
                                        Threshold_cv = 30,
                                        Save_as_Results = "xlsx", 
                                        Save_as_Plot = "svg")
  
  
  a<-Pool_Estimation_res$DF$CV_result
  remove <- c(remove,a[Pool_Estimation_res$DF$CV_result$HighVar,"Metabolite"])
  remove <- unique(remove)
  
  tmp <- ZeroVarCheckRes$Input_data_filtered[,!colnames(ZeroVarCheckRes$Input_data_filtered) %in% remove]
  Pool_Estimation_res$DF$Filtered_Input_data <- tmp
  

  
  # Feature filtering
  # usethis::edit_file(file.path(scriptpath,"FeatureFiltering.R"))   
  source(file.path(scriptpath,"FeatureFiltering.R"))
  FeatureFilteringRes <- FeatureFiltering(Result_Folder = Folder_paths$PostBinnerFolder,
                                          RDS_folder = Folder_paths$PostBinnerRDS,
                                          OutputName = OutputName,
                                          ForceRun = ForceRun,
                                          Input_data=Pool_Estimation_res$DF$Filtered_Input_data,
                                          Input_SettingsFile=Input_SettingsFile_sel,
                                          Input_SettingsInfo=Input_SettingsInfo,
                                          Mode="Modified", # Standard or Modified or None
                                          Feature_Filt_Value = 0.8)
  
  
  # Missing value imputation
  # usethis::edit_file(file.path(scriptpath,"MissingValueImputation.R"))   
  source(file.path(scriptpath,"MissingValueImputation.R"))
  MissingValueImputationRes <- MissingValueImputation(Result_Folder = Folder_paths$PostBinnerFolder,
                                                      RDS_folder = Folder_paths$PostBinnerRDS,
                                                      OutputName = OutputName,
                                                      ForceRun = ForceRun,
                                                      Filtered_matrix=FeatureFilteringRes$Filtered_matrix,
                                                      Input_SettingsFile=Input_SettingsFile_sel,
                                                      Input_SettingsInfo=Input_SettingsInfo,
                                                      Algorithm = "Mean",
                                                      CoRe=CoRe) # or "MAI"
  
  
  
  # Normalisation
  # usethis::edit_file(file.path(scriptpath,"Normalisation.R"))   
  source(file.path(scriptpath,"Normalisation.R"))
  NormalisationRes <- Normalisation(Result_Folder = Folder_paths$PostBinnerFolder,
                                    RDS_folder = Folder_paths$PostBinnerRDS,
                                    scriptpath = scriptpath,
                                    OutputName = OutputName,
                                    ForceRun = ForceRun,
                                    NA_removed_matrix = MissingValueImputationRes$NA_removed_matrix,
                                    Input_SettingsFile = Input_SettingsFile_sel,
                                    Input_SettingsInfo = Input_SettingsInfo,
                                    Method = "PQN", # PQN or TIC
                                    PQN_selected_persentage=40,
                                    CoRe=CoRe) 
  
  
  # Outlier detection
  # usethis::edit_file(file.path(scriptpath,"OutlierDetection.R"))   
  source(file.path(scriptpath,"OutlierDetection.R"))
  OutlierDetectionRes <- OutlierDetection(Result_Folder = Folder_paths$PostBinnerFolder,
                                          RDS_folder = Folder_paths$PostBinnerRDS,
                                          scriptpath= scriptpath,
                                          FolderName = FolderName,
                                          OutputName = OutputName,
                                          ForceRun = ForceRun,
                                          data_norm = NormalisationRes$DF$NormalizedData,
                                          Input_SettingsFile = Input_SettingsFile_sel,
                                          Input_SettingsInfo = Input_SettingsInfo,
                                          OutlierLoop = 2,
                                          npcs=NULL,
                                          HotellinsConfidence = 0.99,
                                          Save_as_Plot = "svg")
  
  
  # Make outlier sample vector
  outlier_sample <- rownames( OutlierDetectionRes$DFs$Processed_data)[ OutlierDetectionRes$DFs$Processed_data$Outliers != "no"]
  
  # Remove outlier samples
  # usethis::edit_file(file.path(scriptpath,"SampleRemoval.R"))   
  source(file.path(scriptpath,"SampleRemoval.R"))
  SampleRemovalRes <- SampleRemoval(Input_data = OutlierDetectionRes$DFs$Processed_data, 
                                    RemoveSamples = RemoveOutlierSamples, 
                                    SampleVector = SampleVector) # or outlier_sample[1:3] for selected.
  processedData <- SampleRemovalRes
  
  
  if(ReplicateSum==TRUE){
    # usethis::edit_file(file.path(scriptpath,"ReplicateSum.R"))   
    source(file.path(scriptpath,"ReplicateSum.R"))
    ReplicateSumres <- ReplicateSum(ProjectFolder = ProjectFolder,
                                    FolderName = FolderName,
                                    Input_data = processedData)
    
    processedData <- ReplicateSumres[["ProcessedDataSummed"]]
    
    # Make 5x5 PC for all conditions
    PCA_data <- processedData[,-c(1:which(names(processedData) == "n_Replicates_Summed"))]
    PCA_SettingsFile <- processedData %>% select("Conditions")
    
  }else{
    # Make 5x5 PC for all conditions
    PCA_data <- processedData[,-c(1:which(names(processedData) == "Outliers"))]
    PCA_SettingsFile <- processedData[,c(1:which(names(processedData) == "Outliers"))]
    
  }
  

  # PCA_data[-grep("hg", tolower(rownames(PCA_data))),]
  # PCA_SettingsFile[-grep("hg", tolower(rownames(PCA_SettingsFile))),]
  
  
  # usethis::edit_file(file.path(scriptpath,"PCA.R"))   
  source(file.path(scriptpath,"PCA.R"))
  
  
  Final_PCA <- suppressWarnings(PCA(Input_data=PCA_data, Plot_SettingsInfo= c(color=  colnames(PCA_SettingsFile)[1] ),
                                    Plot_SettingsFile=PCA_SettingsFile,
                                    OutputPlotName = "Sample PCA",
                                    Save_as_Plot =  NULL))

  
  ggsave(file=paste(Folder_paths$PostBinnerFolder,"/", "Final_PCA.svg", sep=""), plot=Final_PCA, width=20, height=15, unit="cm")
  
  
  # for colname in SSF make a 5x5 PCA and color using the selected column
  # for (colname in 1:length(colnames(PCA_SettingsFile[,c(1:which(names(processedData) == "Outliers"))]))){ # colname = 7
  for (colname in 1:(length(colnames(PCA_SettingsFile)))){ # colname = 1
    
    colnames(PCA_SettingsFile)[colname] <- gsub(" ", "_", colnames(PCA_SettingsFile)[colname])
    
    if(all(is.na(PCA_SettingsFile[colnames(PCA_SettingsFile)[colname]]))==FALSE){
      # usethis::edit_file(file.path(scriptpath,"QC_PCA.R"))   
      source(file.path(scriptpath,"QC_PCA.R"))
      QC_PCA <- suppressWarnings(PCA(Input_data=PCA_data, Plot_SettingsInfo= c(color=  colnames(PCA_SettingsFile)[colname]),
                                     Plot_SettingsFile=PCA_SettingsFile,
                                     Save_as_Plot =  NULL))
      

      # Results_folder_Preprocessing_folder_Quality_Control_PCA_folder = file.path(Folder_paths$PostBinnerFolder, "Quality_Control_PCA")
      # if (!dir.exists(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder)) {dir.create(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder)}
      
      ggsave(file=paste(Folder_paths$PostBinnerFolder,"/", "PCA_5x5_",colnames(PCA_SettingsFile)[colname],".svg", sep=""), plot=QC_PCA, width=60, height=40, unit="cm")
    }
  }
  
  Preprocessing_res <- list("processedData" = processedData)
  return(Preprocessing_res)
}





# #usethis::edit_file("frezzaMetaboAnalysis/FeatureSelection/Scripts/Functions/GroupVennDiagram.R")
#  source("frezzaMetaboAnalysis/FeatureSelection/Scripts/Functions/GroupVennDiagram.R")
# GroupVennDiagramRes <- GroupVennDiagram(ProjectFolder=ProjectFolder,
#                    Input_data = LoadDataRes$Data %>% select(!all_of(c(outlier_features, FeatureFilteringRes$Features_filtered))),
#                    Input_SettingsFile = LoadDataRes$Metadata,
#                    Input_SettingsInfo = c(PoolSamples = "Pool", Conditions="Conditions"),
#                    Ignore_group = NULL,
#                    OutputName = NULL,
#                    Processed_data = processedData,
#                    FolderName = FolderName,
#                    Boxplot = TRUE)


# colnames(processedData)[4:5] <- c("Analytical_Replicates", "Biological_Replicates")
# # Select metabolites

