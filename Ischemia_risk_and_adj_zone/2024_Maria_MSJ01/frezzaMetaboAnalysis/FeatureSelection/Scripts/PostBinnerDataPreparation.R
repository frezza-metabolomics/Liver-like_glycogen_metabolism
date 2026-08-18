
#' @param ProjectFolder String with path to a project folder. (ie. "Z:\1. Mass Spec Projects\RUGARLI\Rrejusha PARAYIL\2023_Rrejusha_RP01")

PostBinnerDataPreparation <- function(Folder_paths = Folder_paths,
                                      scriptpath = scriptpath,
                                      CDInputFolder = NULL
                                      ){
  prepared_data <- list()
  filtering_feats_list <- list()
  
  if (is.null(CDInputFolder)==FALSE){
    CDInputFolder <- paste0("_",CDInputFolder)
  }

  # FInd the resundant features
  for(mode in c("positive", "negative")){ # mode="positive"
    message(paste("Preparing", mode, "dataset")) # message back the loaded file
    # Load binner output
    
    suppressMessages(filtering_feats <- readxl::read_excel(file.path (Folder_paths$ProjectFolder,paste("5.",  Folder_paths$ProjectCode, "Analysed Data"),"Results","PreBinner_folder",paste0("Compounds_", Folder_paths$ProjectCode ,"_data_",mode,"_processed_Report.xlsx")), sheet = 4, col_names = T))
    # Remove the rows with Mass error != NA and !=0
    filtering_feats <- filtering_feats[is.na(filtering_feats$`Mass Error`)==F & filtering_feats$`Mass Error`!=0, "Feature"] %>% unlist() %>% as.vector()
    # save list of redundant features
    write.csv(filtering_feats, file.path(Folder_paths$ProjectFolder,paste("5. ",  Folder_paths$ProjectCode , " Analysed Data","/Results","/PreBinner_folder/RedundantFeatures_",Folder_paths$ProjectCode,"_",mode,".csv" , sep = "")),row.names = F)
    filtering_feats_list[[mode]]<- filtering_feats
    }
  filtering_feats_list <- unlist(filtering_feats_list, recursive = FALSE)
  
  # For each mode find CD data
  file_listData <- list.files(file.path(Folder_paths$ProjectFolder,paste("3.", Folder_paths$ProjectCode, "CD","Data"),paste0(Folder_paths$ProjectCode, CDInputFolder)), pattern = ".xlsx", full.names = TRUE) # Lists all the excel files in the tracefinder folder
  if(is_empty(file_listData) == TRUE){stop("The list of xlsxfiles neg or pos are empty. please check that the correct folder is used and that there are xlxs files")}
  if(length(grep("hilic", tolower(file_listData)))>0 ){file_listData <- file_listData[grep("hilic", tolower(file_listData))[1]] }
  
  data <- readxl::read_xlsx(file_listData) %>% as.data.frame()
  rownames(data) <- make.unique(paste0(data$`m/z`,"_" ,data$`RT [min]`))
    
  # Remove redundant features from the data
  data <- data[!rownames(data) %in% filtering_feats_list,]

  message("Output location: ", file.path (Folder_paths$PreBinnerFolder,paste0("Data_clean_", Folder_paths$ProjectCode,".xlsx" )))
  # Save the cleaned data
  writexl::write_xlsx(data, file.path (Folder_paths$PreBinnerFolder,paste0("Data_clean_", Folder_paths$ProjectCode,".xlsx" )))#,showNA = TRUE)
}


  