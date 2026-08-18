


LoadPackages <- function(){
  
  RequiredPackages <- c("tidyverse", # general scripting
                        "readxl", # to read xls files
                        "Hmisc", # correlation
                        "patchwork", # for output plot grid
                        "matrixStats", # to get median for signal drift
                        "MAI", # general scripting
                        "factoextra", # visualize PCA
                        "qcc", # for hotelling plots
                        "ggplot2", # For visualization PCA
                        "hash", # Dictionary in R for making column of outliers
                        "gridExtra",
                        "inflection",# For finding inflection point/ Elbow knee /PCA component selection # https://cran.r-project.org/web/packages/inflection/inflection.pdf # https://deliverypdf.ssrn.com/delivery.php?ID = 454026098004123081018105104090015093000085002012023032095093077109069092095000114006057018122039107109012089110120018031068078025094036037013095100070100076109026029024044005068010070117123085122016083112098002109001027028000024115096122101001083084026&EXT = pdf&INDEX = TRUE # https://arxiv.org/abs/1206.5478
                        "patchwork", # for output plot grid
                        "gtools",
                        "EnhancedVolcano"
  )
  
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])] #Check which packages are already installed
  if (length(new.packages)) install.packages(new.packages) # install missing packages
  suppressWarnings(suppressMessages(library(tidyverse))) # load tidyverse package
  
  
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install(new.packages)
  }
  
}
