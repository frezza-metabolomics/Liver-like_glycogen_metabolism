# Function to load the raw data from Ming
LoadData <- function(Input_data=data, Sheet = 2, SSF_DATA = NULL, ProjectName = NULL){
 # Input_data=file.path(pathtoFile,"20220620_Luisa_LB01_pHILIC_E240_updated.xlsx")
 # Sheet = 2
 # SSF_DATA = input_SSFXLX
 # ProjectName = ProjectName
 # Load libraries
 # Function to load the raw data from Ming
library(tidyverse)
library(readxl)
if (is.null(SSF_DATA)){
  Input_data <- read_excel(Input_data, sheet = Sheet)
  Input_data <- as.data.frame(Input_data)
  PooledData <- as.data.frame(subset(Input_data, grepl("Pool", Input_data[,1])))
  Input_data <- Input_data[1:min(which(is.na(Input_data[,1]))) - 1,] # take rows until you reach the first NA
  }else{
  Input_data <- read_excel(Input_data, sheet = Sheet, col_names = F)
  Input_data <- as.data.frame(Input_data)
  PooledData <- as.data.frame(subset(Input_data, grepl("Pool", Input_data[,1])))
}

#fix the pooled data
colnames(PooledData) <- Input_data[1,] # make first row into col names
PooledData <- PooledData[order(PooledData$Code),]
rownames(PooledData) <- PooledData[,1]# make 1st column into row names
PooledData <- PooledData[,-1] # delete fist row

#fix the input data
colnames(Input_data) <- Input_data[1,] # make first row into col names
Input_data <- Input_data[-1,] # delete first row
Input_data<- Input_data[order(Input_data[,1]),]
rownames(Input_data) <- Input_data[,1]# make 1st column into row names
Input_data <- Input_data[,-1] # delete fist row
Input_data <- subset(Input_data, !grepl("Pool", row.names(Input_data)))
#names(Input_data)[1] <- "Sample"
# Create the Experiments_design
if (is.null(SSF_DATA)) {
Design <- as.data.frame(Input_data$Sample)
names(Design) <- "Conditions"
rownames(Design) <- rownames(Input_data)
}else{
if (is.null(ProjectName)){
stop("impossible to having seperate SSF and empty Project Name, please assigne the 'ProjectName' value")
}else{
Design <- read_excel(SSF_DATA, sheet = 1)
Design <- as.data.frame(Design)

# find the column we want
sel.col <- grep("Vial Identifier", Design, ignore.case = FALSE, value = FALSE)
# find the row we want
sel.row <- grep("Vial Identifier", t(Design[,grep("Vial Identifier", Design, ignore.case = FALSE, value = FALSE)]), ignore.case = FALSE, value = FALSE)
# select the Samples 
# We will also take the Conditions as we will use them for stratified random ordering of samples
ncols = 3 # change this depending on the number of column we take from the SSF
Design <- Design[sel.row:length(as.data.frame(t(Design))),sel.col:(sel.col+ ncols)] # take the sample data
if(any(is.na(Design))==TRUE){
Design <- Design[1:min(which(is.na(Design[,1]))) - 1,] # take rows until you find the first NA
}
Design <- Design[-c(1:4),] # remove the first 4 rows with information
# make the design column name vector
COLNAMES <- c("Sample", "Conditions", "Tissue", "Species")
if (length(COLNAMES)!= (ncols+1)){
stop("The selected columns and the names in COLNAMES have different lengths")
}
names(Design) <- COLNAMES


Design <- Design[!is.na(names(Design))]
Design$Conditions <- gsub(" ", "_", Design$Conditions)
Design$Conditions <- gsub("_℃", "oC",Design$Conditions)
Design$Conditions <- gsub("℃", "oC",Design$Conditions)
row.names(Design) <- Design$Sample
Design$Sample <- NULL
}
}
Input_data$Sample <- NULL
Input_data$`total ion count` <- NULL
Input_data$TIC <- NULL

Input_data <- as.matrix(mutate_all(as.data.frame(Input_data), function(x) as.numeric(as.character(x)))) %>% as.data.frame()

# identify weird metabolite names
weirdMetabolites<- names(Input_data)[grepl("[*]",names(Input_data))]
#weirdMetabolites <- c(weirdMetabolites , names(Input_data)[grepl("-",names(Input_data))])
weirdMetabolites <- c(weirdMetabolites, names(Input_data)[grepl("/",names(Input_data))])

write(weirdMetabolites, file="metabolites_to_check.txt")
names(Input_data) <- gsub("[*]", "", names(Input_data))
#names(Input_data) <- gsub("-", "_", names(Input_data))
names(Input_data) <- gsub("/", ".", names(Input_data))
#remove column from data frame filled with NAs
PooledData <- PooledData[, unlist(lapply(PooledData, function(x) !all(is.na(x))))]
#Make result
#result <- list(Data = Input_data, Experimental_design = Design)
print("Done")
# pkg <-c("package:tidyverse","package:readxl")
# lapply(pkg, detach, character.only = TRUE, unload = TRUE)
return(list(Data = Input_data, Experimental_design = Design, weirdMetabolites = data.frame(weirdMetabolites), Pooled_data = PooledData))
}
