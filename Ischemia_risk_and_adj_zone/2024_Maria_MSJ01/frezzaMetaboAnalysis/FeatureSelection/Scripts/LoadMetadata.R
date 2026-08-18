

LoadMetadata <- function(Folder_paths,
                         scriptpath,
                         data
                         ){
  

  # Get Metadata
  file_listSSF <- list.files(file.path(Folder_paths$ProjectFolder,paste("1.", Folder_paths$ProjectCode, "SSFs")), pattern = ".xlsx", full.names = TRUE) # get SSF excel file
  # remove open xls
  if(sum(grep("~", file_listSSF))>0){
    file_listSSF <- file_listSSF[-grep("~", file_listSSF)]
  }
  
  message(paste("Input SSF file:",  strsplit(file_listSSF, " SSFs/")[[1]][2])) # message the file that is loaded
  suppressMessages(SSF_data <- readxl::read_excel(file_listSSF, sheet = 1, col_names = T) %>% as.data.frame()) # load ssf file
  sel.col <- grep("Vial Identifier", SSF_data, ignore.case = FALSE, value = FALSE)  # find the column we want
  sel.row <- grep("Vial Identifier", t(SSF_data[,grep("Vial Identifier", SSF_data, ignore.case = FALSE, value = FALSE)]), ignore.case = FALSE, value = FALSE)  # find the row we want
  metadata <- SSF_data[sel.row:length(as.data.frame(t(SSF_data))),sel.col:ncol(SSF_data)] # subset the data
  if (sum(is.na(metadata[,1])) != 0){
    metadata <- metadata[1:min(which(is.na(metadata[,1]))) - 1,]} # Remove NA rows in the end
  colnames(metadata) <- metadata[1,] # make first row into col names
  metadata <- metadata[-c(1:4),] # remove the first 4 rows with information
  if(sum(grep("Example", metadata[1,]))>0 ){
    metadata <- metadata[-1,] # remove the first 4 rows with information
  }
  rownames(metadata) <- metadata[,1]  # make 1st column into row names
  metadata <- metadata[,-1] # delete fist column
  colnames(metadata)[1] <- "Conditions"
  
  # Remove any pools or bkanks from the SSF
  metadata <- metadata[!grepl("pool", tolower(rownames(metadata)),ignore.case=TRUE ) ,] 
  metadata <- metadata[!grepl("blank|blk|blnk", tolower(rownames(metadata)),ignore.case=TRUE ) ,] 
  
  # Add back inthe SSF the pool samples
  Pool_metadata <- data.frame(matrix(NA, nrow = sum(grepl("Pool", colnames(data),ignore.case=TRUE)), ncol = ncol(metadata))) # Mame a dataframe for pool metadata
  colnames(Pool_metadata) <- colnames(metadata) # add colnames from metadata
  rownames(Pool_metadata) <-  colnames(data)[grepl("Pool", colnames(data),ignore.case=TRUE)] # add rownames of the pool samples from the data
  Pool_metadata$Conditions <- "Pool" # add Pool as condition
  Output_metadata <- rbind(metadata, Pool_metadata) # merge the metadata with the pool metadata
  
  
  # Get sample run order
  # select XXXX folder in 2. CD Data
  tmp <- list.dirs(file.path(Folder_paths$ProjectFolder,paste("2.", Folder_paths$ProjectCode, "Raw Data")), recursive = F)
  #tmp <- tmp[grepl(paste0("_",Folder_paths$ProjectCode,"_"),tmp)]              
  
  file_listRunOrder <- list.files(file.path(tmp), pattern = ".csv",recursive = TRUE, full.names = TRUE) # get all csv files in the raw data folder and subfolders
  if (length(file_listRunOrder)>1){
    ordeR <- c(2,3)
    warning("more than 2 csv order lists were found. we will select the ", ordeR, "please change the entry ordeR accordingly to avoid an error" )
    RandomizedSamples <- NULL
    for (i in ordeR){
      RandomizedSample <- read.csv(file_listRunOrder[i]) %>% as.data.frame() 
      colnames(RandomizedSample) <- RandomizedSample[1,] # make first row into col names
      RandomizedSample <- RandomizedSample[-c(1),] # remove the first 4 rows with information
      RandomizedSample <- RandomizedSample[,1]
      RandomizedSamples <- c(RandomizedSamples, RandomizedSample)
    }
  }else{
    message(paste("Sample Run Order file:",  strsplit(file_listRunOrder, "E240/")[[1]][2]))
    RandomizedSamples <- read.csv(file_listRunOrder) %>% as.data.frame() 
    colnames(RandomizedSamples) <- RandomizedSamples[1,] # make first row into col names
    RandomizedSamples <- RandomizedSamples[-c(1),] # remove the first 4 rows with information
    RandomizedSamples <- RandomizedSamples[,1]
  }
  
  LoadMetadataRes = list("Metadata" =  Output_metadata, "RunOrder" = RandomizedSamples)
  
  return(LoadMetadataRes)
}
