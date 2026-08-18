

DifferentialAbundantAnalysis<- function(Folder_paths,
                                        scriptpath,
                                        Result_Folder,
                                        RDS_folder,
                                        DMA_data,
                                        DMA_SettingsFile,
                                        DMA_SettingsInfo = c(Conditions = "Conditions"),
                                        ForceRun = ForceRun, #ForceRun = T
                                        Paired = FALSE,
                                        PairedData = NULL,
                                        OutputFolderName = NULL,
                                        OutputName = NULL,
                                        LogTransfold = FALSE,
                                        CoRe = FALSE
                                        ){
  
  
  # Create output folders for DMA
  if(is.null(OutputFolderName)!=TRUE){
    Result_Folder <- paste(Folder_paths$DMAFolder,OutputFolderName, sep="_" )
    RDS_folder <- paste(Folder_paths$DMARDS,OutputFolderName, sep="_" )
    if (!dir.exists(Result_Folder)) {dir.create(Result_Folder)}
    if (!dir.exists(RDS_folder)) {dir.create(RDS_folder)}
  }else{
    if (!dir.exists(Result_Folder)) {dir.create(Result_Folder)}
    if (!dir.exists(RDS_folder)) {dir.create(RDS_folder)}
  }
    
  
  # # SubSet data
  # seldatarownames <- DMA_SettingsFile %>% filter(DMA_SettingsFile$Conditions %in% DMA_SettingsInfo) %>% rownames() 
  # DMA_SettingsFile <- DMA_SettingsFile %>% filter( rownames(DMA_SettingsFile) %in%  seldatarownames)
  # DMA_data <- DMA_data %>% filter( rownames(DMA_data) %in%  seldatarownames)
  
    
  ####----------------- Zero variance check for features --------------####
  # usethis::edit_file(file.path(scriptpath,"ZeroVarianceCheck.R"))   
  source(file.path(scriptpath,"ZeroVarianceCheck.R"))
  ZeroVarCheckRes <- ZeroVarCheck(Result_Folder =Result_Folder,
                                  Input_data= DMA_data, 
                                  OutputName = OutputName,
                                  RemoveFeatures = TRUE)
  
  DMA_data <- ZeroVarCheckRes$Input_data_filtered
  
  if(LogTransfold==TRUE){
    originalData <- DMA_data
    DMA_data <- log2(DMA_data)}
  

  if(CoRe==T){
    set.seed(1234)
    DMA_data <- DMA_data %>%
      mutate_all(~ ifelse(. == 0, 1 + runif(1, min = 0.01, max = 0.1), .))
  }
  
  
  # usethis::edit_file(file.path(scriptpath,"DMADataCheck.R"))
  source((file.path(scriptpath,"DMADataCheck.R")))
  DMADataCheck(RDS_folder = RDS_folder,
               DMA_data =DMA_data,
               DMA_SettingsFile = DMA_SettingsFile,
               DMA_SettingsInfo = DMA_SettingsInfo,
               OutputName = OutputName)
  
  

  # usethis::edit_file(file.path(scriptpath,"Shapiro.R"))
  source((file.path(scriptpath,"Shapiro.R")))
  ShapiroTestRes <- ShapiroTest(Result_Folder = Result_Folder,
                                RDS_folder = RDS_folder,
                                scriptpath= scriptpath,
                                
                                DMA_data=DMA_data,
                                DMA_SettingsFile=DMA_SettingsFile,
                                DMA_SettingsInfo=DMA_SettingsInfo,
                                CoRe=FALSE,
                                ForceRun = ForceRun,
                                QQplots=FALSE,
                                Save_as_Plot="svg",
                                Save_as_Results="csv",
                                FolderName=FolderName,
                                OutputName=OutputName)
  

  # usethis::edit_file(file.path(scriptpath,"Bartlett.R"))
  source((file.path(scriptpath,"Bartlett.R")))
  BartlettTestRes <- BartlettTest(Result_Folder = Result_Folder,
                                  RDS_folder = RDS_folder,
                                  scriptpath= scriptpath,
                                  DMA_data=DMA_data,
                                  DMA_SettingsFile=DMA_SettingsFile,
                                  DMA_SettingsInfo=DMA_SettingsInfo,
                                  Save_as_Plot="svg",
                                  ForceRun = ForceRun,
                                  Save_as_Results="csv",
                                  Plot=TRUE,
                                  FolderName=FolderName,
                                  OutputName=OutputName)
  
  
  # usethis::edit_file(file.path(scriptpath,"Log2FC.R"))
  source((file.path(scriptpath,"Log2FC.R")))
  # if logtransform is true we use the original data to take the logfold change becaue otherwise we would take the log of log and that value would be very small
  if(LogTransfold==TRUE){
    Log2FCRes <- Log2FC(Result_Folder = Result_Folder,
                        RDS_folder =RDS_folder,
                        scriptpath= scriptpath,
                        DMA_data=originalData,
                        DMA_SettingsFile=DMA_SettingsFile,
                        DMA_SettingsInfo = DMA_SettingsInfo,
                        ForceRun = ForceRun,
                        FolderName=FolderName,
                        OutputName=OutputName,
                        CoRe=CoRe)
  }else{
    Log2FCRes <- Log2FC(Result_Folder =Result_Folder,
                        RDS_folder = RDS_folder,
                        scriptpath= scriptpath,
                        DMA_data=DMA_data,
                        DMA_SettingsFile=DMA_SettingsFile,
                        DMA_SettingsInfo = DMA_SettingsInfo,
                        ForceRun = ForceRun,
                        FolderName=FolderName,
                        OutputName=OutputName,
                        CoRe=CoRe)
  }

  
  
  
  ## IF we have 1 comparison
  if (length(DMA_SettingsInfo["numerator"]) & length(DMA_SettingsInfo["denominator"]) & is.na(DMA_SettingsInfo["numerator"])==FALSE) {
    # If we have parametric data
    if (mean(as.numeric(ShapiroTestRes$DF$Shapiro_result[,3]))>50){
       # If we have equal variances
      if (sum(BartlettTestRes$DF[,3], na.rm = T)/ (length(BartlettTestRes$DF[,3])- sum(is.na(BartlettTestRes$DF[,3]))) > 0.50){
        message("Using t-test with equal variances")
        # t- test with variance homogeneity
        # usethis::edit_file(file.path(scriptpath,"ttest.R"))
        source((file.path(scriptpath,"ttest.R")))
        DMARes <- tTest(Result_Folder = Result_Folder,
                        RDS_folder = RDS_folder,
                        scriptpath= scriptpath,
                        DMA_data=DMA_data,
                        DMA_SettingsFile=DMA_SettingsFile,
                        DMA_SettingsInfo = DMA_SettingsInfo,
                        Log2FC_table=Log2FCRes$Log2FC_table,
                        ForceRun = ForceRun,
                        STAT_padj = "fdr",
                        Paired = Paired,
                        PairedData = PairedData,
                        FolderName=FolderName,
                        OutputName=OutputName,
                        variance_equality=TRUE,
                        PCutoff= 0.05,
                        FCCutoff= 1,
                        CoRe=CoRe)
        DMARes <- DMARes[["DF"]][["DMA_result"]]
        
      }else{ 
        # If we do NOT have equal variances
        message("Using t-test with unequal variances")
        # t- test withOUT variance homogeneity
        # usethis::edit_file(file.path(scriptpath,"ttest.R"))
        source((file.path(scriptpath,"ttest.R")))
        DMARes <- tTest(Result_Folder = Result_Folder,
                        RDS_folder = RDS_folder,
                        scriptpath= scriptpath,
                        DMA_data=DMA_data,
                        DMA_SettingsFile=DMA_SettingsFile,
                        DMA_SettingsInfo = DMA_SettingsInfo,
                        Log2FC_table=Log2FCRes$Log2FC_table,
                        ForceRun = ForceRun,
                        STAT_padj = "fdr",
                        Paired = Paired,
                        PairedData = PairedData,
                        FolderName=FolderName,
                        OutputName=OutputName,
                        variance_equality=FALSE,
                        PCutoff= 0.05,
                        FCCutoff= 1,
                        CoRe=CoRe)
        DMARes <- DMARes[["DF"]][["DMA_result"]]
        
      }
    }else{     
      # IF we do not have parametric data
      message("Using w-test")
      # Wilcoxon test  = single comparison non parametric
      # usethis::edit_file(file.path(scriptpath,"Wilcoxtest.R"))
      source((file.path(scriptpath,"Wilcoxtest.R")))
      DMARes <- wTest(Result_Folder = Result_Folder,
                        RDS_folder = RDS_folder,
                        scriptpath= scriptpath,
                        DMA_data=DMA_data,
                        DMA_SettingsFile=DMA_SettingsFile,
                        DMA_SettingsInfo = DMA_SettingsInfo,
                        Log2FC_table=Log2FCRes$Log2FC_table,
                        STAT_padj = "fdr",
                        FolderName=FolderName,
                        OutputName=OutputName)
      DMARes <- DMARes[["DF"]][["DMA_result"]]
      
    }
  }else{
    # We have multiple comparisons
    
    # ANOVA  = multiple comparison  parametric
    # usethis::edit_file(file.path(scriptpath,"anova.R"))
    source((file.path(scriptpath,"anova.R")))
    AnovaTestRes <- AnovaTest(DMA_data=DMA_data,
                              DMA_SettingsFile=DMA_SettingsFile,
                              DMA_SettingsInfo = c(Conditions = "Conditions", numerator = NULL, denominator = NULL),
                              Log2FC_table=Log2FCRes$Log2FC_table,
                              STAT_padj = "fdr",
                              PCutoff= 0.05,
                              FCcutoff= 1,
                              CoRe=FALSE,
                              OutputName="",
                              ForceRun = ForceRun,
                              Save_as_Plot = "svg",
                              Labels = FALSE,
                              ProjectFolder = ProjectFolder,
                              FolderName="HG_All_vs_All")
    DMARes <- AnovaTestRes$DF$DMA_result
    
    
  }
    

  

  

   
## Welch ANOVA with  = multiple comparison parametricunequal variances
  # usethis::edit_file("frezzaMetaboAnalysis/FeatureSelection/Scripts/Functions/welchanova.R")
#   source("frezzaMetaboAnalysis/FeatureSelection/Scripts/Functions/welchanova.R")
# WelchAnovaTestRes <- WelchAnovaTest(DMA_data,
#                      DMA_SettingsFile=DMA_SettingsFile,
#                      DMA_SettingsInfo = c(Conditions="Conditions", numerator = NULL, denumerator = NULL),
#                      Log2FC_table=Log2FCRes$Log2FC_table,
#                      PCutoff= 0.05,
#                      FCcutoff= 1,
#                      CoRe=FALSE,                     
#                      STAT_padj = "fdr",
#                      Labels = TRUE,
#                      ForceRun = T,
#                      Save_as_Plot = "svg",
#                      ProjectFolder = ProjectFolder,
#                      FolderName=FolderName)
# 
# wTestRes <- WelchAnovaTestRes$DF$DMA_result
  ## Kruskal wallis = multiple comparison non parametric
  # usethis::edit_file("frezzaMetaboAnalysis/FeatureSelection/Scripts/Functions/kruskal.R")
  # source("frezzaMetaboAnalysis/FeatureSelection/Scripts/Functions/kruskal.R")
  # KruskalTestRes <- KruskalTest(DMA_data=DMA_data,
  #                           DMA_SettingsFile=DMA_SettingsFile,
  #                           DMA_SettingsInfo = c(Conditions = "Conditions", numerator = NULL, denominator = NULL),
  #                           Log2FC_table=Log2FCRes$Log2FC_table,
  #                           STAT_padj = "fdr",
  #                           PCutoff= 0.05,
  #                           FCcutoff= 1,
  #                           CoRe=FALSE,
  #                           OutputName="",
  #                           ForceRun = T,
  #                           Save_as_Plot = "svg",
  #                           Labels = TRUE,
  #                           ProjectFolder = ProjectFolder,
  #                           FolderName=FolderName)
  

  return(DMARes)
}







# ProjectCode <- strsplit(ProjectFolder, "_")[[1]][[ length(strsplit(ProjectFolder, "_")[[1]])]] # Get the project code
# Analysed_data_folder <- file.path(ProjectFolder,paste("5.", ProjectCode, "Analysed Data", sep = " ")) # Define Analysed data folder
# if (!dir.exists(Analysed_data_folder)) {dir.create(Analysed_data_folder)} # check if the Analysed data folder exists and Make it
# Results_folder <- file.path(Analysed_data_folder,  # Define results folder
#                             if (is.null(FolderName)) {paste("Results") }else{ paste("Results",FolderName,sep = "_" )}) # Add additional folder name to the output results folder name
# if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # check if results folder and if it does not make it
# RDS_folder <- file.path(Analysed_data_folder, # Define RDS folder 
#                         if (is.null(FolderName)) {paste("RDS")}else{ paste("RDS",FolderName,sep = "_" )})  # Add additional folder name to the output results folder name
# if (!dir.exists(RDS_folder)) {dir.create(RDS_folder)} # check if RDS folder and if it does not make it
# 
# df <- data.frame(name = character(),
#                  `m/z` = numeric(),
#                  RT = numeric(),
#                  stringsAsFactors = FALSE)
# 
# file_listData <- list.files(file.path(ProjectFolder,paste("3.", ProjectCode, "CD","Data"),ProjectCode), pattern = ".xlsx", full.names = TRUE) # Lists all the excel files in the tracefinder folder
# if(is_empty(file_listData) == TRUE){stop("The list of xlsxfiles neg or pos are empty. please check that the correct folder is used and that there are xlxs files")}
# for (file in file_listData){ # file = file_listData[1]
#   if(length(grep("neg", tolower(file)))==1){OutputName <- "neg"}else{ OutputName <- "pos"}
#   # Get data
#   data <- readxl::read_xlsx(file)
#   data <- as.data.frame(data)
#   data <- cbind(data %>% select("Name","m/z", "RT [min]"), data %>% select(grep("Area",colnames(data))))
#   colnames(data)[c(2,3)] <- c("m/z", "RT")
#   
#   # There might be duplicated rows the best would be to kinda merge them. Now I will do the folowwing
#   # if there are duplicated rows take the one with less NA
#   count_na <- function(x) {
#     sum(is.na(x))
#   }
#   data <- data %>%
#     mutate(na_count = apply(., 1, count_na)) %>%  # Apply count_na function to each row
#     group_by(`m/z`, RT) %>%
#     filter(na_count == min(na_count)) %>%
#     slice(1) %>%  # If there are ties, keep only the first row
#     ungroup() %>%
#     select(-na_count)  # Remove the temporary na_count column
#   
#   data <- as.data.frame(data)
#   
#   
#   rownames(data) <- paste0(data$`m/z`,"_" ,data$RT)
#   # add valine-d8
#   if(data[which.min(abs(data$`m/z` - 126.13645)),3] - 5.67< 0.2){
#     valine <- rownames(data)[which.min(abs(data$`m/z` - 126.13645))]
#     rownames(data)[which.min(abs(data$`m/z` - 126.13645))] <- "Valine-d8"
#   }
#   df <- rbind(df, data%>% select(Name, `m/z`, RT))
# }
# df$Metabolite <-rownames(df)
# 
# result <- merge(wdata, df, "Metabolite")