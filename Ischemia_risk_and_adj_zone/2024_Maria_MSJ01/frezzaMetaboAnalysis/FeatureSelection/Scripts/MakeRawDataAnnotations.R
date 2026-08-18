# This script loads all excel files in CD folder and combines all annotations into one final excel.

MakeRawDataAnnotations <- function(Folder_paths,
                                   CDInputFolder = CDInputFolder,
                                   OutputName = NULL,
                                   ForceRun = FALSE
){
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)) {OutputName = ""} else{OutputName = paste("_",OutputName, sep="")}
  if (file.exists(file.path(Folder_paths$PreBinnerRDS, paste( "MakeRawDataAnnotations", OutputName,".Rdata", sep = ""))) == TRUE & ForceRun == FALSE) { # Check if RDS exists and forceRun is false then skip script and load result from RDS, otherwise run it
    message("Loading already existing ", paste( "MakeRawDataAnnotations", OutputName,".Rdata", sep = ""))
    load(file.path(Folder_paths$PreBinnerRDS, paste( "MakeRawDataAnnotations", OutputName,".Rdata", sep = "")), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running MakeRawDataAnnotations script main body.")
    
    if (is.null(CDInputFolder)==FALSE){
      CDInputFolder <- paste0("_",CDInputFolder)
    }
    
    ### ---------------------------- 1. Merge all files in CD Data folder -------------------------- ###
    # Lists all the excel files in the CD folder
    file_listData <- list.files(file.path(Folder_paths$ProjectFolder,paste("3.", Folder_paths$ProjectCode, "CD","Data"),paste0(Folder_paths$ProjectCode, CDInputFolder)), pattern = ".xlsx", full.names = TRUE) 
    if(is_empty(file_listData) == TRUE){stop("Unable to find CDInputFolder files. Check the 3.CD Data folder and the CDInputFolder parameter.")}
    # Reorder the list with "pHILIC" first
    if (length(which(grepl("pHILIC", file_listData))) > 0) {
      file_listData <- c(file_listData[which(grepl("pHILIC", file_listData))], file_listData[-which(grepl("pHILIC", file_listData))])
    }
    
    # Make list for processed input dfs
    input_df_list <- list()
    for (file in file_listData){ # file = file_listData[1]
      
      # Read the Excel file
      df <- readxl::read_excel(file)
      # Save file name
      file_name <- sub(".*Compounds_NR01_(.*?)\\.xlsx", "\\1", file)
      # Load file and parse the data
      if("Name" %in% colnames(df) & "m/z"%in% colnames(df) & "RT [min]"%in% colnames(df)){
        df <- df %>%
          select("Name", "m/z", "RT [min]") %>%
          rename(!!file_name := "Name",
                 "mz"= "m/z",
                 "RT"= "RT [min]" ) 
        
        df$mz <- as.numeric(df$mz)
        df$RT <- as.numeric(df$RT)
        
        input_df_list[[file_name]] <- df
      }else{
        warning(paste("Warning: Columns Name, mz and RT were NOT found in file:", file))
      }
      df$mz <- as.null(df$mz)
      df$RT <- as.null(df$RT)
    }# end of file loop
    
    # Merge all data frames in the list by 'mz' and 'RT'
    if (length(input_df_list) > 0) {
      merged_data <- reduce(input_df_list, full_join, by = c("mz", "RT"))
      merged_data <- merged_data[order(merged_data$mz),]
    }else{
      warning("No CD Input data were found")
    }
    # Make sure the first column is called name
    colnames(merged_data)[1] <- "name"
    
    ###############################################
    ### ### Here add the metaname converter ### ###
    ###############################################
    
    ### ---------------------------- 2. Manage annotations -------------------------- ###
    
    ## to improve: Instead of just adding the annotations we should use the MNC and compare them and validate with the ms2.
    ### Add the ms1 annotations
    # Identify the columns
    philic_col <- grep("philic", tolower(colnames(merged_data)))
    hmdb_col <- grep("hmdb", tolower(colnames(merged_data)))
    emdb_col <- grep("emdb", tolower(colnames(merged_data)))
    leachables_col <- grep("leachables", tolower(colnames(merged_data)))
    
    # Create the 'name' column using coalesce, with warnings for missing columns
    if (length(philic_col) == 1 & length(hmdb_col) == 1) {
      merged_data <- merged_data %>%
        mutate(name = coalesce(
          .[[philic_col]], 
          .[[hmdb_col]], 
          if (length(emdb_col) == 1) .[[emdb_col]] else NA, 
          if (length(leachables_col) == 1) .[[leachables_col]] else NA
        ))
      
    } else {
      # Warnings for missing or ambiguous EMDB and Leachables columns
      if (length(emdb_col) == 0) {
        warning("The 'emdb' column was not found.")
      } else if (length(emdb_col) > 1) {
        warning("More than one column contains 'emdb'. Please check the column names.")
      }
      
      if (length(leachables_col) == 0) {
        warning("The 'leachables' column was not found.")
      } else if (length(leachables_col) > 1) {
        warning("More than one column contains 'leachables'. Please check the column names.")
      }
      
      # # Separate warnings for missing or ambiguous "philic" and "hmdb" columns
      # if (length(philic_col) == 0) {
      #   warning("The 'philic' column was not found.")
      # } else if (length(philic_col) > 1) {
      #   warning("More than one column contains 'philic'. Please check the column names.")
      # }
      
      if (length(hmdb_col) == 0) {
        warning("The 'hmdb' column was not found.")
      } else if (length(hmdb_col) > 1) {
        warning("More than one column contains 'hmdb'. Please check the column names.")
      }
    }
    
    
    ### Add the ms2 annotations
    vault_col <- grep("vault", tolower(colnames(merged_data)))
    cloud_col <- grep("cloud", tolower(colnames(merged_data)))
    # Second block: Handle vault and cloud
    if (length(vault_col) == 1) {
      merged_data <- merged_data %>%
        mutate(name = coalesce(name, .[[vault_col]]))
    } else {
      if (length(vault_col) == 0) {
        warning("The 'vault' column was not found.")
      } else if (length(vault_col) > 1) {
        warning("More than one column contains 'vault'. Please check the column names.")
      }
    }
    
    if (length(cloud_col) == 1) {
      merged_data <- merged_data %>%
        mutate(name = coalesce(name, .[[cloud_col]]))
    } else {
      if (length(cloud_col) == 0) {
        warning("The 'cloud' column was not found.")
      } else if (length(cloud_col) > 1) {
        warning("More than one column contains 'cloud'. Please check the column names.")
      }
    }
    
    
    ### ---------------------- 3. Make Raw data dataframe --------------------------- ###
    
    # Read an Excel file
    df <- readxl::read_excel(file_listData[1])
    # Load file and parse the data
    if("Name" %in% colnames(df) & "m/z"%in% colnames(df) & "RT [min]"%in% colnames(df)){
      df <- df %>%
        select(-"Name") %>%
        rename("mz"= "m/z",
               "RT"= "RT [min]") %>% 
        select("mz","RT","Reference Ion", grep("Area",colnames(.)))  %>%
        select(-grep("(Max.)",colnames(.)))
    }
    
    #Merge the data with the annotations
    final_data <- merge(merged_data%>% select(name, mz, RT), df, by = c("mz", "RT")) 
    # fix column names
    colnames(final_data) <-sub(paste0(".*", Folder_paths$ProjectCode), Folder_paths$ProjectCode,  colnames(final_data))
    colnames(final_data)  <- sub(".raw.*", "",  colnames(final_data) )
    # Remove blank samples
    final_data <- final_data[,!grepl("blank|blk|blnk", tolower(colnames(final_data) ))] # remove any blank samples
    final_data <- final_data[order(final_data$mz),]
    
    # Save data
    writexl::write_xlsx(final_data, file.path(Folder_paths$PreBinnerFolder, "Raw_Data_Annotated.xlsx", sep = ""))#,showNA = TRUE)
    # Save Rdata
    save(final_data, file = file.path(Folder_paths$PreBinnerRDS, paste("MakeRawDataAnnotations", OutputName,".Rdata", sep = "")))
  } # End of calculations
  # Return data
  return(Raw_Data_Annotated = final_data) 
} # Script End