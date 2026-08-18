

FixFeatureSplits <- function(Folder_paths,
                             OutputName = NULL,
                             ForceRun = FALSE,
                             
                             complementary_score_threshold = 0.95,
                             # Correlation_threshold = 0.8, not used
                             Median_Intensity_distance_threshold = 0.2,
                             Input_Data = NULL
){

  ## 1. ------------ Check outputs ----------- ##
  if (is.null(OutputName)) {OutputName = ""} else{OutputName = paste("_",OutputName, sep="")}
  if (file.exists(file.path(Folder_paths$PreBinnerRDS, paste( "FixFeatureSplits", OutputName,".Rdata", sep = ""))) == TRUE & ForceRun == FALSE) { # Check if RDS exists and forceRun is false then skip script and load result from RDS, otherwise run it
    message("Loading already existing ", paste( "FixFeatureSplits", OutputName,".Rdata", sep = ""))
    load(file.path(Folder_paths$PreBinnerRDS, paste( "FixFeatureSplits", OutputName,".Rdata", sep = "")), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running FixFeatureSplits script main body.")
    
  final_data <-  Input_Data
  #final_data <- readxl::read_excel( file.path(Folder_paths$PreBinnerFolder,"Raw_Data_Annotated.xlsx"))
  final_data <- final_data[order(final_data$mz),]
  
  ### ---------------------- 3. Fix feature duplications --------------------------- ###
  #fix issues 1. feature split in multiple rows (same mz) ,2. Select only 1 from all adducts for a feature, 3. metabolite found in both modes (same mw)
  
  
  ####  1. Subset data into lists with the same mz and rt (same mz)------------ ####
  
  # Initialize an empty list to store results
  duplicate_list_full <- list()
  
  # Loop through each unique rounded mz value
  for (mz_value in unique(round(final_data$mz, 3))) {   #    mz_value = unique(round(final_data$mz, 3))[153]
    # Filter rows with the same rounded mz value
    sublist <- final_data[round(final_data$mz, 3) == mz_value, ]
    
    # If there are at least 2 rows, check RT differences
    if (nrow(sublist) > 1) {
      # Get RT values
      rt_values <- sublist$RT
      
      # Calculate pairwise absolute differences
      rt_diffs <- outer(rt_values, rt_values, FUN = function(x, y) abs(x - y))
      
      # Get indices of pairs where the RT difference is less than 0.5 (excluding the diagonal)
      indices <- which(rt_diffs < 0.1, arr.ind = TRUE)
      
      # Unique pairs of rows that have RT differences less than 0.5
      unique_pairs <- unique(indices[indices[, 1] < indices[, 2], ]) %>% as.matrix()
      
      if(length(unique_pairs)==2){
        unique_pairs <- unique_pairs %>% t()
      }
      
      # Check if unique_pairs has rows
      if (length(unique_pairs) > 0) {
        # Create a list for each unique RT pair
        for (pair in 1:nrow(unique_pairs)) {
          row1 <- unique_pairs[pair,1]
          row2 <- unique_pairs[pair,2]
          rt_sublist <- sublist[c(row1, row2), ]
          
          # Use mean mz and mean rt as the list name
          list_name <- paste("mz =", round(mean(rt_sublist$mz), 3), "RT =", round(mean(rt_sublist$RT), 1))
          duplicate_list_full[[list_name]] <- rt_sublist
        }
      }
    }
  }
  
  ## apply some tests to find the actual duplicates features
  duplicate_list_short <- list()
  # loop though each list and determine if its duplicationor not
  for (duplicate in names(duplicate_list_full)){ #  duplicate =  names(duplicate_list_full)[26]
    list_selected <- duplicate_list_full[[duplicate]]
    
    ### Calculate how complementary they are
    #Count col values
    values_score <- sum(colSums(!is.na(list_selected  %>% select(-mz, -RT, -name, -`Reference Ion`)))==1) / length(colnames(list_selected  %>% select(-mz, -RT, -name, -`Reference Ion`)))
    #Count col nas
    nas_score <-  sum(colSums(is.na(list_selected  %>% select(-mz, -RT, -name, -`Reference Ion`)))==1) / length(colnames(list_selected  %>% select(-mz, -RT, -name, -`Reference Ion`)))
    
    #test complementarity
    if(values_score>complementary_score_threshold & nas_score > complementary_score_threshold){
      complementarity_TestPass=TRUE
    }else{
      complementarity_TestPass=FALSE
    }
    
    ### Compare the max intensity if they are roughly the same
    # Calculate the row maximum for the selected columns, excluding specified columns
    row_max1 <- list_selected[1, ] %>% select(-mz, -RT, -name, -`Reference Ion`)%>% unlist()%>% median(., na.rm=T)
    row_max2 <-list_selected[2, ] %>% select(-mz, -RT, -name, -`Reference Ion`)%>% unlist()%>% median(., na.rm=T)
    
    Median_Int_distance <- abs(row_max1-row_max2)/row_max1
    
    # test max Intensity distance
    if(Median_Int_distance < Median_Intensity_distance_threshold ){
      Median_Intensity_distance=TRUE
    }else{
      Median_Intensity_distance=FALSE
    }
    
    # Features have to be on the same mode
    
    
    if( length( grep('\\]\\-' , list_selected[1,"Reference Ion"])) ==    length( grep('\\]\\-' , list_selected[2,"Reference Ion"])) &
        length(   grep('\\]\\+' , list_selected[1,"Reference Ion"])) ==   length(   grep('\\]\\+' , list_selected[2,"Reference Ion"])) ){
      Same_mode = TRUE
    }else{
      Same_mode = FALSE
    }
    
    # if a feature pais passes the tests then they are added to the short list
    if (complementarity_TestPass==TRUE & Median_Intensity_distance==TRUE &  Same_mode == TRUE){
      duplicate_list_short[[duplicate]] <- list_selected
    }
  }# end of looping though duplicated list
  
  
  # # Plots
  for (duplicate in names(duplicate_list_short)){ #  duplicate =  names(duplicate_list_short)[1]
    list_selected <- duplicate_list_short[[duplicate]]
    
    #If mz and rt are exactly the same add a "_2" to the name in order to plot
    if(list_selected[1,2]== list_selected[2,2]){
      list_selected[2,2] <- list_selected[2,2] +0.0001
    }
    
     # Makie plots
     ## Plots
     #Make visual representation of sample complmentarization ?
     # Step 3: Create a binary dataframe to represent NA (1 for present, 0 for missing)
     df_na <- list_selected %>%
       pivot_longer(cols = starts_with(Folder_paths$ProjectCode), names_to = "Variable", values_to = "Value") %>%
       mutate(NA_status = ifelse(is.na(Value), 0, 1))  # 0 for NA, 1 for present

     # Step 4: Create the plot
     featureSplitHeatmap <- ggplot(df_na, aes(x = Variable, y = factor(RT), fill = factor(NA_status))) +
       geom_tile(color = "white") +  # Create heatmap tiles
       scale_fill_manual(values = c("0" = "red", "1" = "blue"),
                         labels = c("Missing", "Present"),
                         name = "Data Status") +  # Red for missing, blue for present
       labs(x = "Variables", y = "Row RT", title = paste("mz = ", round(mean(list_selected$mz), 4))) +
       theme_minimal() +
       theme(legend.position = "top",
           axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels by 45 degrees
     # save plot
     ggsave(file.path(Folder_paths$FeatureSplitsFolder, paste0(list_selected[1,1],"_", list_selected[1,2],"_", "Heatmap" ,".svg")),plot=featureSplitHeatmap, width=10, height=8) # save the volcano plot
      ## Plots
     # Make plot for visual check of values
     df_long <- list_selected %>%
       pivot_longer(cols = starts_with(Folder_paths$ProjectCode), names_to = "Variable", values_to = "Value")

     # Create the plot
     featureSplitDotplot <- ggplot(df_long, aes(x = Variable, y = Value, color = factor(RT))) +
       geom_line(aes(group = RT), size = 1) +  # Connect points with lines
       geom_point(size = 3) +                  # Add points
       #  scale_color_manual(values = c("red", "blue", "yellow")) +  # Set colors
       labs(x = "Variables", y = "Values", color = "Row RT", title = paste("mz = ", round( mean(df_long$mz),4))) +
       theme_minimal() +
       theme(legend.position = "top",
             axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels by 45 degrees
     # save plot
     ggsave(file.path(Folder_paths$FeatureSplitsFolder, paste0(list_selected[1,1],"_", list_selected[1,2],"_" , "Dotplot",".svg")),plot=featureSplitDotplot, width=10, height=8) # save the volcano plot
     
     }
  
  ############################################################
  #### Here add Theo's code to create a filter set for CD ####
  ############################################################
  
  fixed_duplicated_features <- list()
  removed_duplicated_features <- list()
  ## fix feature duplications
  for (duplicate in names(duplicate_list_short)){ #  duplicate =  names(duplicate_list_short)[1]
    list_selected <- duplicate_list_short[[duplicate]]
    
    # if features are found in positive mode
    if(length(grep('\\]\\+', list_selected[,"Reference Ion"]) )>0){
      # Look for M+-H and keep that annotation and replace NAs for the values of the other row
      if(sum(list_selected$`Reference Ion` == '[M+H]+1')==1){
        list_selected[which(list_selected$`Reference Ion` == '[M+H]+1'), ] <- ifelse(is.na(list_selected[which(list_selected$`Reference Ion` == '[M+H]+1'), ]), 
                                                                                     list_selected[which(list_selected$`Reference Ion` != '[M+H]+1'), ], 
                                                                                     list_selected[which(list_selected$`Reference Ion` == '[M+H]+1'), ])
        
        
        feature_corrected <-list_selected[which(list_selected$`Reference Ion` == '[M+H]+1'), ]
        feature_removed <-list_selected[which(list_selected$`Reference Ion` != '[M+H]+1'), ]
        
        
        fixed_duplicated_features[[paste("mz =", feature_corrected$mz, "RT =",feature_corrected$RT)]] <- feature_corrected
        removed_duplicated_features[[paste("mz =", feature_removed$mz, "RT =",feature_removed$RT)]] <- feature_removed
        
      }else{# if both are M+H or non
        #if both areM+H then keep the one with least NAs
        if(sum(is.na(list_selected[1,]))>  sum(is.na(list_selected[2,]) )){
          keep_row =2
        }else{
          keep_row =1
        }
        list_selected[keep_row,] <- ifelse(is.na( list_selected[keep_row,] ), 
                                           list_selected[-keep_row,] , 
                                           list_selected[keep_row,] )
        
        feature_corrected <-list_selected[keep_row, ]
        feature_removed <-list_selected[-keep_row, ]
        
        fixed_duplicated_features[[paste("mz =", feature_corrected$mz, "RT =",feature_corrected$RT)]] <- feature_corrected
        removed_duplicated_features[[paste("mz =", feature_removed$mz, "RT =",feature_removed$RT)]] <- feature_removed
        
      }
    }else{ #if mode = negative
      # Look for M+-H and keep that annotation and replace NAs for the values of the other row
      if(sum(list_selected$`Reference Ion` == '[M-H]-1')==1){
        list_selected[which(list_selected$`Reference Ion` == '[M-H]-1'), ] <- ifelse(is.na(list_selected[which(list_selected$`Reference Ion` == '[M-H]-1'), ]), 
                                                                                     list_selected[which(list_selected$`Reference Ion` != '[M-H]-1'), ], 
                                                                                     list_selected[which(list_selected$`Reference Ion` == '[M-H]-1'), ])
        
        
        feature_corrected <-list_selected[which(list_selected$`Reference Ion` == '[M-H]-1'), ]
        feature_removed <-list_selected[which(list_selected$`Reference Ion` != '[M-H]-1'), ]
        
        
        fixed_duplicated_features[[paste("mz =", feature_corrected$mz, "RT =",feature_corrected$RT)]] <- feature_corrected
        removed_duplicated_features[[paste("mz =", feature_removed$mz, "RT =",feature_removed$RT)]] <- feature_removed
        
      }else{# if both are M-H or non
        #if both areM+H then keep the one with least NAs
        if(sum(is.na(list_selected[1,]) )>  sum(is.na(list_selected[2,]) )){
          keep_row =2
        }else{
          keep_row =1
        }
        list_selected[keep_row,] <- ifelse(is.na( list_selected[keep_row,] ), 
                                           list_selected[-keep_row,] , 
                                           list_selected[keep_row,] )
        
        feature_corrected <-list_selected[keep_row, ]
        feature_removed <-list_selected[-keep_row, ]
        
        fixed_duplicated_features[[paste("mz =", feature_corrected$mz, "RT =",feature_corrected$RT)]] <- feature_corrected
        removed_duplicated_features[[paste("mz =", feature_removed$mz, "RT =",feature_removed$RT)]] <- feature_removed
      }
    }
  }
  
  
  # Substitute the corrected duplicated features in the data
  for(feature in names(fixed_duplicated_features)){  # feature = names(fixed_duplicated_features)[1]
    feat_data <- fixed_duplicated_features[[feature]]
    final_data[ which(final_data$mz==  feat_data$mz & final_data$RT==  feat_data$RT & final_data$`Reference Ion`==  feat_data$`Reference Ion`),] <- feat_data
  }
  
  # Remove the removed duplicated features in the data
  for(feature in names(removed_duplicated_features)){  # feature = names(removed_duplicated_features)[1]
    feat_data <- removed_duplicated_features[[feature]]
    final_data <- final_data[ -(which(final_data$mz==  feat_data$mz & final_data$RT==  feat_data$RT & final_data$`Reference Ion`==  feat_data$`Reference Ion`)),]
    final_data <- as.data.frame(final_data)
  }
  
  # # Save data
  # writexl::write_xlsx(final_data, file.path(Folder_paths$PreBinnerFolder, "Raw_Data_Annotated_SplitRemoved.xlsx", sep = ""))#,showNA = TRUE)
  # # Save Rdata
   save(final_data, file = file.path(Folder_paths$PreBinnerRDS, paste("FixFeatureSplits", OutputName,".Rdata", sep = "")))
  }  # End of calculations
  # Return data
  return(Raw_Data_Annotated_SplitRemoved = final_data)
} # Script End


# 
# if(sum(duplicated(paste0(data$`m/z`,"_" ,data$RT)))>0){ # check for duplicated feature names mz_RT
#   message("There are ", sum(duplicated(paste0(data$`m/z`,"_" ,data$RT))) , " duplicated features in the data.")
# }
# rownames(data) <- make.unique(paste0(data$`m/z`,"_" ,data$RT))
# 
# 
# # add valine-d8
# if(OutputName == "pos"){
#   if(data[which.min(abs(data$`m/z` - 126.13645)),3] - 5.67< 0.2){
#     valine <- rownames(data)[which.min(abs(data$`m/z` - 126.13645))]
#     rownames(data)[which.min(abs(data$`m/z` - 126.13645))] <- "Valine-d8"
#   }
# }


# Output_data <- final_data





# # first reload a df and get the mw column
# df <- readxl::read_excel(file)
# if("Name" %in% colnames(df) & "m/z"%in% colnames(df) & "RT [min]"%in% colnames(df)){
#   df <- df %>%
#     select("Name", "m/z", "RT [min]", "Calc. MW") %>%
#     rename(!!file_name := "Name",
#            "mz"= "m/z",
#            "RT"= "RT [min]",
#            "mw"= "Calc. MW") 
# }
#   
# # 1 add the mw column
# adduct_data <- merge(df%>% select(mz, RT, mw), final_data, by=c("mz", "RT"))
# # 2 separate the pos and neg dataset
# pos_adduct_data <- adduct_data[-grep("]-", adduct_data$`Reference Ion`),]
# neg_adduct_data <- adduct_data[grep("]-", adduct_data$`Reference Ion`),]  
# 
# Multiple_neg_adduct_data <- list()
# # loop though each list and determine if its duplicationor not
# for (mw_value in unique(round(neg_adduct_data$mw, 3))) { #  mw_value = unique(round(neg_adduct_data$mw, 3))[1]
#   sublist <- neg_adduct_data[round(neg_adduct_data$mw, 3) == mw_value, ]
#   if(nrow(sublist)>1 & sublist[1,"Reference Ion"] != sublist[2,"Reference Ion"]){
#     Multiple_neg_adduct_data[[paste0("mw = ", mw_value)]] <- sublist
#     
#       cor_matrix <- Hmisc::rcorr(  t(sublist %>% select(-mz, -RT, -mw,-name, -`Reference Ion`)), type = "spearman")
#       cor_matrix <- as.data.frame(cor_matrix$r)
#       print(cor_matrix)
#   }
# }
# 
# Multiple_pos_adduct_data <- list()
# # loop though each list and determine if its duplicationor not
# for (mw_value in unique(round(pos_adduct_data$mw, 3))) { #  mw_value = unique(round(pos_adduct_data$mw, 3))[1]
#   sublist <- pos_adduct_data[round(pos_adduct_data$mw, 3) == mw_value, ]
#   if(nrow(sublist)>1 & sublist[1,"Reference Ion"] != sublist[2,"Reference Ion"]){
#     
#     
#     Multiple_pos_adduct_data[[paste0("mw = ", mw_value)]] <- sublist
#     
#   }
# }
# 

# 3 find the same mw and decide what to keep


## compare pos and neg and keep a metabolite only one mode





# ### Take correlations
# # Calculate the row maximum for the selected columns, excluding specified columns
# # sample_data[] <- lapply(sample_data, function(x) as.numeric(as.character(x)))
# 
# list_selected %>% select(-mz, -RT, -name, -`Reference Ion`)
# 
# list_selected <- a
# #We need at least 2 pairs of values to take correlations
# if(sum(colSums(!is.na(list_selected  %>% select(-mz, -RT, -name, -`Reference Ion`)))==2) >2){
#   cor_matrix <- Hmisc::rcorr(  t(list_selected %>% select(-mz, -RT, -name, -`Reference Ion`)), type = "spearman")
#   cor_matrix <- as.data.frame(cor_matrix$r)
# 
#   # test correlation
#   if(cor_matrix[1,2] > Correlation_threshold ){
#     Correlation_TestPass=TRUE
#   }else{
#     Correlation_TestPass=FALSE
#   }
# }else{
# 
# }


