

RawDataPreparation <- function(Folder_paths,
                               scriptpath,
                               CDInputFolder = NULL, # This is the added folder name inside the '3. CD data' folder in order to take the correct CD data
                               Input_SettingsInfo = c(Conditions = "Conditions",PoolSamples = "Pool"),
                               OutputName = NULL,
                               ForceRun = FALSE
                               ){
  # # Define the output file path
  # output_file <- file.path(Folder_paths$PreBinnerFolder, "RawDataPreparation_output.log")
  # 
  # # Open a connection to the file
  # log_connection <- file(output_file, open = "wt")
  # 
  # # Redirect stdout to the file and console
  # sink(log_connection, split = TRUE)
  # 
  # # Redirect stderr to the file
  # log_connection_stderr <- file(output_file, open = "wt")
  # sink(log_connection_stderr, type = "message")
  # 
  

####------------ Create Raw Data with Combined Annotations from CD -----------####
# usethis::edit_file(file.path(scriptpath,"MakeRawDataAnnotations.R"))
source(file.path(scriptpath,"MakeRawDataAnnotations.R"))
Raw_Data_Annotated <- MakeRawDataAnnotations(Folder_paths=Folder_paths,
                                             CDInputFolder = CDInputFolder,
                                             OutputName = OutputName,
                                             ForceRun = ForceRun)
  
####------------------ Remove feature splits from CD ------------------####
# usethis::edit_file(file.path(scriptpath,"FixFeatureSplits.R"))
source(file.path(scriptpath,"FixFeatureSplits.R"))
Raw_Data_Annotated_SplitRemoved <- FixFeatureSplits(Folder_paths=Folder_paths,
                                                    OutputName = OutputName,
                                                    ForceRun = ForceRun,
                                                    complementary_score_threshold = 0.9,
                                                    Median_Intensity_distance_threshold = 0.2,
                                                    Input_Data = Raw_Data_Annotated)
  
  
####-------------- Check for duplicated features. Same mz and rt. And Make unique ---------####
if(sum(duplicated(paste0(Raw_Data_Annotated_SplitRemoved$mz,"_" ,Raw_Data_Annotated_SplitRemoved$RT)))>0){ 
  message("There are ", sum(duplicated(paste0(Raw_Data_Annotated_SplitRemoved$mz,"_" ,Raw_Data_Annotated_SplitRemoved$RT))) , " duplicated features in the data.")
  # Add unique IDs to features = mz _ rt _ dupication
  rownames(Raw_Data_Annotated_SplitRemoved) <- make.unique(paste0(Raw_Data_Annotated_SplitRemoved$mz,"_" ,Raw_Data_Annotated_SplitRemoved$RT))
}else{
  rownames(Raw_Data_Annotated_SplitRemoved) <- make.unique(paste0(Raw_Data_Annotated_SplitRemoved$mz,"_" ,Raw_Data_Annotated_SplitRemoved$RT))
}

# Make rownames into column in order to save
Raw_Data_Annotated_SplitRemoved <- tibble::rownames_to_column(Raw_Data_Annotated_SplitRemoved, var = "ID")

# Save data
writexl::write_xlsx(Raw_Data_Annotated_SplitRemoved, file.path(Folder_paths$PreBinnerFolder, "Raw_Data_Ready.xlsx", sep = ""))#,showNA = TRUE)
# Save Rdata
save(Raw_Data_Annotated_SplitRemoved, file = file.path(Folder_paths$PreBinnerRDS, paste("Raw_Data_Annotated_Ready", OutputName,".Rdata", sep = "")))


# # Stop redirection
# sink()               # Stop stdout redirection
# sink(type = "message")  # Stop stderr redirection
# 
# # Close the file connections
# close(log_connection)
# close(log_connection_stderr)
# 
# # Notify user
# cat("Console output has been saved to:", output_file, "\n")


return(Raw_Data_Annotated_SplitRemoved)
} # END
