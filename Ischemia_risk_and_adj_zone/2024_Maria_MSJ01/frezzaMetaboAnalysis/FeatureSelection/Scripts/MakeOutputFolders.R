# This script takes a ProjectFOlder path and makes output folders


MakeOutputFolders <- function(ProjectFolder,
                              OutputFolderName = NULL
                             #4 CDInputFolder=NULL
                              ){
  
  # Get project Code
  ProjectCode <- "MSj01" # Get the project code

  ### Make folders
  Analysed_data_folder <- file.path(ProjectFolder,paste("5.", ProjectCode, "Analysed Data", sep = " "))
  if (!dir.exists(Analysed_data_folder)) {dir.create(Analysed_data_folder)}
  
  if(is.null(OutputFolderName)==TRUE){
    ForCollaborators_folder <- file.path(Analysed_data_folder, "ForCollaborators") 
    Results_folder <- file.path(Analysed_data_folder, "Results")
  }else{
    ForCollaborators_folder <- file.path(Analysed_data_folder, paste0("ForCollaborators_",OutputFolderName) )
    Results_folder <- file.path(Analysed_data_folder, paste0("Results_",OutputFolderName) )
  }
  
  
  
  #Collaborators folders
  #ForCollaborators_folder <- file.path(Analysed_data_folder, "ForCollaborators") 
  if (!dir.exists(ForCollaborators_folder)) {dir.create(ForCollaborators_folder)} 
  
  # Data_ForCollaborators_folder <- file.path(ForCollaborators_folder, "Data") 
  # if (!dir.exists(Data_ForCollaborators_folder)) {dir.create(Data_ForCollaborators_folder)} 
  # 
  # Plots_ForCollaborators_folder <- file.path(ForCollaborators_folder, "Plots") 
  # if (!dir.exists(Plots_ForCollaborators_folder)) {dir.create(Plots_ForCollaborators_folder)} 
  #
  #Our folders
 # Results_folder <- file.path(Analysed_data_folder, "Results")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  
  PreBinner_folder <- file.path(Results_folder, "PreBinner_folder")
  if (!dir.exists(PreBinner_folder)) {dir.create(PreBinner_folder)}

  PostBinner_folder <- file.path(Results_folder, "PostBinner_folder")
  if (!dir.exists(PostBinner_folder)) {dir.create(PostBinner_folder)}
  
  DMA_folder <- file.path(Results_folder, "DMA_folder")
 # if (!dir.exists(DMA_folder)) {dir.create(DMA_folder)}

  PreBinner_RDS <- file.path(Results_folder, "PreBinner_RDS")
  if (!dir.exists(PreBinner_RDS)) {dir.create(PreBinner_RDS)}

  PostBinner_RDS <- file.path(Results_folder, "PostBinner_RDS")
  if (!dir.exists(PostBinner_RDS)) {dir.create(PostBinner_RDS)}
  
  DMA_RDS <- file.path(Results_folder, "DMA_RDS")
  #if (!dir.exists(DMA_RDS)) {dir.create(DMA_RDS)}

  # 
  # RawData_folder <- file.path(Results_folder, "1. RawData") 
  # if (!dir.exists(RawData_folder)) {dir.create(RawData_folder)}  
  # 
  FeatureSplits_folder <- file.path(PreBinner_folder, "Feature Splits")
  if (!dir.exists(FeatureSplits_folder)) {dir.create(FeatureSplits_folder)}
  # 
  # Binner_folder <- file.path(RawData_folder, "Binner") 
  # if (!dir.exists(Binner_folder)) {dir.create(Binner_folder)}  
  # 
  # # rawDataPoolEstimation_folder <- file.path(RawData_folder, "Pool Estimation") 
  # # if (!dir.exists(rawDataPoolEstimation_folder)) {dir.create(rawDataPoolEstimation_folder)}  
  # 
  # QualityControl_folder <- file.path(Results_folder, "2. QualityControl") 
  # if (!dir.exists(QualityControl_folder)) {dir.create(QualityControl_folder)}  
  # 
  # Preprocessing_folder <- file.path(Results_folder, "3. Preprocessing") 
  # if (!dir.exists(Preprocessing_folder)) {dir.create(Preprocessing_folder)}  
  # 
  # DifferentialMetaboliteAnalysis_folder <- file.path(Results_folder, "4. DifferentialMetaboliteAnalysis") 
  # if (!dir.exists(DifferentialMetaboliteAnalysis_folder)) {dir.create(DifferentialMetaboliteAnalysis_folder)} 
  # 
  # #RDS folders
  # RDS_folder <- file.path(Analysed_data_folder, "RDS") 
  # if (!dir.exists(RDS_folder)) {dir.create(RDS_folder)} 
  
  # Return paths of relevant folders
  return(list(
    ProjectFolder = ProjectFolder,
    ProjectCode = ProjectCode,
    AnalysedDataFolder = Analysed_data_folder,
    ForCollaborators_folder = ForCollaborators_folder,
    # Data_ForCollaborators_folder = Data_ForCollaborators_folder,
    # Plots_ForCollaborators_folder = Plots_ForCollaborators_folder,
    ResultsFolder = Results_folder,
    PreBinnerFolder = PreBinner_folder,
    PostBinnerFolder = PostBinner_folder,
    PreBinnerRDS = PreBinner_RDS,
    PostBinnerRDS = PostBinner_RDS,
    #RawDataFolder = RawData_folder,
    FeatureSplitsFolder = FeatureSplits_folder,
    #BinnerFolder = Binner_folder,
    #rawDataPoolEstimationFolder = rawDataPoolEstimation_folder,
    #QualityControlFolder = QualityControl_folder,
    #PreprocessingFolder = Preprocessing_folder,
    #DifferentialMetaboliteAnalysisFolder = DifferentialMetaboliteAnalysis_folder,
    #RDSFolder = RDS_folder
    DMAFolder = DMA_folder,
    DMARDS = DMA_RDS
  ))
}
