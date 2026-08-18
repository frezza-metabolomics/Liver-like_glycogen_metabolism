
#' @param ProjectFolder String with path to a project folder. (ie. "Z:\1. Mass Spec Projects\RUGARLI\Rrejusha PARAYIL\2023_Rrejusha_RP01")
#' @param FolderName String which is added to the results and RDS output folder names. (ie. "test" results in Results_test and RDS_test)
#' @param OutputName String which is added to the RSD and output files of the function.  (ie. "test" results in LoadDataRes_test.RSD)
#' @param ForceRun TRUE or FALSE. If TRUE the script runs and overwrites results and RDS. IF FALSE and RDS exists script doenot run and loads RDS.


PoolEstimation <- function(Result_Folder,
                           RDS_folder,
                           scriptpath,
                           
                           OutputName = NULL,
                           ForceRun = TRUE,
                           COLORS = "Conditions",
                           Input_data,
                           Input_SettingsFile = NULL,
                           Input_SettingsInfo = NULL,
                           Sample_Order = NULL,
                           
                           Remove_unstable_features = TRUE,
                           Threshold_cv = 30,
                           
                           Save_as_Results = "xlsx", 
                           Save_as_Plot = "svg"
                           ){
  # Load dependencies
  # usethis::edit_file(file.path(scriptpath,"Functions/PCA.R"))   
  source(file.path(scriptpath,"PCA.R"))
  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)) {OutputName = ""} else{OutputName = paste("_",OutputName, sep="")}
  if (file.exists(file.path(RDS_folder, paste( "PoolEstimationRes", OutputName,".Rdata", sep = ""))) == TRUE & ForceRun == FALSE) { # Check if RDS exists and forceRun is false then skip script and load result from RDS, otherwise run it
    message("Loading already existing ", paste( "PoolEstimationRes", OutputName,".Rdata", sep = ""))
    load(file.path(RDS_folder, paste( "PoolEstimationRes", OutputName,".Rdata", sep = "")), .GlobalEnv)
    
  ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running PoolEstimationRes script main body.")
  
  ## 5. -------- Do Pool Coefficient of variation Estimation ------ ##
  pool_data <- Input_data[rownames(Input_SettingsFile)[Input_SettingsFile$Conditions == "Pool"], ]
  # Calculate the Coefficient of Variation
  result_df <- apply(pool_data, 2,  function(x) { (sd(x, na.rm =T)/  mean(x, na.rm =T))*100 }  ) %>% t()%>% as.data.frame()
  rownames(result_df)[1] <- "CV"
  # reorder 
  result_df <- result_df[, names(pool_data)]
  
  # calculate the percentage of  NAs
  NAvector <- apply(pool_data, 2,  function(x) { (sum(is.na(x))/length(x))*100 }  )#%>% t()%>% as.data.frame()
  # Make the final dataframe
  result_df_final <- result_df %>% t() %>% as.data.frame()%>% mutate(HighVar = CV > Threshold_cv)# %>% as.data.frame()
  result_df_final$MissingValuePercentage <- NAvector
  rownames(result_df_final)<- colnames(Input_data)
  result_df_final <- rownames_to_column(result_df_final,"Metabolite" )
  
  # Save results
  if(is.null(Save_as_Results)==FALSE){
    if (Save_as_Results == "xlsx") {
      xlsDMA <- file.path(Result_Folder, paste("Pool_CV_table", OutputName, ".xlsx", sep = ""))   # Save the DMA results table
      writexl::write_xlsx(result_df_final, xlsDMA, col_names = TRUE) # save the DMA result DF
    } else if (Save_as_Results == "csv") {
      csvDMA <- file.path(Result_Folder, paste("Pool_CV_table", OutputName, ".csv", sep = ""))
      write.csv(result_df_final, csvDMA, row.names = FALSE, na = "NA") # save the DMA result DF with NAs as NA
    } else if (Save_as_Results == "txt") {
      txtDMA <- file.path(Result_Folder, paste("Pool_CV_table", OutputName, ".txt", sep = ""))
      write.table(result_df_final, txtDMA, col.names = TRUE, row.names = FALSE, na = "NA") # save the DMA result DF with NAs as NA
    }
    
  }

  # Make unstable metabolites/High var filtered dataset
  filtered_Input_data <- Input_data
  if(Remove_unstable_features ==TRUE){
    
    singleValue_metabs <- result_df_final[is.na(result_df_final$HighVar),"Metabolite"]
    if(length(singleValue_metabs)>0){
      message(length(singleValue_metabs), " features have a single value, thus are removed.")
      write.table(singleValue_metabs,row.names =  FALSE,col.names = FALSE, file = paste(Result_Folder,"/SingleValueFilteredFeatures", OutputName,".csv",sep =  ""))
      filtered_Input_data <- filtered_Input_data %>% select(!singleValue_metabs)
      }
    
    result_df_final <- result_df_final[is.na(result_df_final$HighVar)==FALSE,] # remove the singlevalue features
    unstable_metabs <- result_df_final[result_df_final$HighVar== TRUE,"Metabolite"] # get metabs with high var

    
    if(length(unstable_metabs)>0){
      message(paste(length(unstable_metabs),"metabolites have high variance and are removed."))
      write.table(unstable_metabs,row.names =  FALSE,col.names = FALSE, file = paste(Result_Folder,"/HighVarianceFilteredFeatures", OutputName,".csv",sep =  ""))
     
      filtered_Input_data <- filtered_Input_data %>% select(!unstable_metabs)
      
      }else{
      message("No metabolites had high variance")
    }
  }
  

  # library(tidyverse)
  # generatedMatrix <- as.matrix(t(Input_data))
  # # Replace NA with 0
  # generatedMatrix[is.na(generatedMatrix)] <- 1
  # 
  # # Convert all non-NA values to 1
  # generatedMatrix[generatedMatrix != 1] <- 0
  # 
  # generatedMatrix %>% as.vector %>% 
  #   tibble(value = ., row = rep(1:nrow(generatedMatrix), times = ncol(generatedMatrix)),
  #          col = rep(1: ncol(generatedMatrix), each = nrow(generatedMatrix))) %>%
  #   ggplot(aes(x = row, y = col, fill = value)) +
  #   geom_tile(size = 2) +
  #   scale_fill_gradient(low = 'black',high = 'white')+
  #   theme_minimal() +
  #   theme(legend.position = 'none')
  
  
  
  
  ## 6. ---------------------- Make plots PCA, Hist, Violin -------------------------- ##
  # Start QC plot list
  pool_plot_list <- list()
  
  # All sample PCA
  pca_data <- merge(Input_SettingsFile %>% select(Conditions), filtered_Input_data, by = 0) %>% column_to_rownames("Row.names")
  if(length( as.vector(  which(is.na(colSums(pca_data[,-1])))))>0){ # Ignore column with NAs from PCA
    warning(paste(length( as.vector(  which(is.na(colSums(pca_data[,-1]))))), "out of",dim(pca_data)[2]-1 , "metabolites were ignored for PCA as they contained missing values."))
    pca_data <- pca_data %>% select(-c(as.vector(  which(is.na(colSums(pca_data[,-1])))  ) +1) )
    }
  pca_data <- pca_data %>% 
    mutate(Sample_type = case_when(Conditions == Input_SettingsInfo[["PoolSamples"]] ~ "Pool", TRUE ~ "Sample"))
  if (ncol(pca_data %>% dplyr::select(where(is.numeric))) == 0){
    print("empty pca, no metabolites with values across all samples")
    plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="PC2 0%", xlab="PC1 0%", main = "empty plot as a placeholder")
  }else if (ncol(pca_data %>% dplyr::select(where(is.numeric))) == 1){
  
    colnames(pca_data)[2] <- "Value"
    pca_plot <-  pca_data %>% ggplot(aes(x="", y=Value), color = Sample_type) +
      geom_violin() +
      geom_jitter(aes(color = Sample_type))
      
  
  } else {
    pca_plot <- suppressWarnings(PCA(Input_data=pca_data %>% dplyr::select(where(is.numeric)), Plot_SettingsInfo= c(color=COLORS),
                                     Plot_SettingsFile = pca_data, OutputPlotName = "QC PCA", Save_as_Plot =  NULL))
  }
  suppressWarnings(  pool_plot_list[["QC_PCA"]] <- pca_plot )# add plot to plotlist
  
  if (is.null(Save_as_Plot) == FALSE) {
    suppressWarnings(ggsave(filename = paste0(Result_Folder,"/QC_PCA", OutputName,".",Save_as_Plot, sep = ""), plot = pca_plot, width = 10,  height = 8))
    }
  
  # Pool CV Histogram
  HistCV <- suppressWarnings(invisible(ggplot(result_df_final, aes(CV)) +
                                        geom_histogram(aes(y=after_stat(density)), color="black", fill="white")+
                                        geom_vline(aes(xintercept=Threshold_cv),
                                                   color="darkred", linetype="dashed", size=1)+
                                        geom_density(alpha=.2, fill="#FF6666") +
                                        labs(title="Coefficient of Variation for metabolites of Pool samples",x="Coefficient of variation (CV%)", y = "Frequency", caption = paste0(length(result_df_final$CV), " total features"))+
                                        theme_classic()))
  if (length(grep("aline-d8", result_df_final$Metabolite ))>0){
    HistCV <- HistCV + 
      geom_vline(aes(xintercept=  result_df_final[ grep("aline-d8", result_df_final$Metabolite ),2]),color="blue", linetype="dashed", size=1) + 
      scale_x_continuous(breaks=c(30,50,100,result_df_final[ grep("aline-d8",  result_df_final$Metabolite ),2]), labels=c("30",'50','100',"Valine-d8 (IS)"))
    
  }
 
  pool_plot_list[["Pool_CV_Histogram"]] <- HistCV # add plot to plotlist
  
  if (is.null(Save_as_Plot) == FALSE){
    suppressMessages(ggsave(filename = paste0(Result_Folder, "/Pool_CV_Histogram",OutputName,".",Save_as_Plot, sep = ""), plot = invisible(HistCV), width = 8,  height = 8))
    }
  
  ### 7. -------------------------------------- Signal drift ------------------------- ###
  if(is.null(Sample_Order)==FALSE){
    Drift_data <- merge(Input_SettingsFile %>% select(Conditions), filtered_Input_data, by = 0, sort=FALSE)
    colnames(Drift_data)[1] <- "File.Name"
    
    rownames(Drift_data) <- Drift_data$File.Name
    
    
    order_selected <- as.vector(as.factor(Sample_Order[Sample_Order %in% Drift_data$File.Name])) # get the sample order vector
    indices <- match(order_selected, Drift_data$File.Name) # find the matching order
    Drift_data_ordered <- Drift_data[indices, ] # Reorder the dataframe 
    Drift_data_ordered$File.Name <- NULL
    Drift_data_ordered[, -1] <- sapply(Drift_data_ordered[, -1], as.numeric) # make all values numeric
    row_medians <- matrixStats::rowMedians(as.matrix(Drift_data_ordered[, -1]), na.rm = TRUE)    # Get sample median
    
    ## Signal drift all data bar plot
    plotdata <- data.frame("SampleMedian" = row_medians, SampleGroup = Drift_data_ordered$Conditions) # get the data to plot
    plotdata <- rownames_to_column(plotdata, var = "RowNames") # add row names as a separate column
    plotdata$RowNames <- factor(plotdata$RowNames, levels = plotdata$RowNames)
    
    signaldrift_plot_all <- ggplot(plotdata, aes(x = RowNames, y = SampleMedian, fill = SampleGroup)) +
      geom_col(position = "dodge") +
      labs(title = "Signal drift, Sample Groups", x = "Samples by acquirement order", y = "Sample Medians") +
      theme_minimal()  + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
      theme(plot.title = element_text(size=22))
    
    signaldrift_plot_all  <- signaldrift_plot_all + geom_smooth(data = filter(plotdata, SampleGroup == "Pool"), # add trend line to plot
                                  aes(x = as.numeric(RowNames), y = SampleMedian),
                                  se = FALSE,method = "lm",  color = "red")
    
    pool_plot_list[["Signal drift, Sample Groups"]] <- signaldrift_plot_all
    
    if (is.null(Save_as_Plot) == FALSE){ # save plot
      suppressMessages(suppressWarnings(ggsave(filename = paste0(Result_Folder, "/SignalDrift_SampleGroup",OutputName,".",Save_as_Plot, sep = ""), plot = signaldrift_plot_all, width = 12,  height = 10)))
    }
    
    ## Median Intensity by group
    totalInt_data <- plotdata[order(plotdata$SampleGroup), ]
    totalInt_data$RowNames <- factor(totalInt_data$RowNames, levels = totalInt_data$RowNames)
    mean_median_by_group <- totalInt_data %>% group_by(SampleGroup) %>% dplyr::summarize(mean_median = mean(SampleMedian))
    
    mean_median_by_group <-  mean_median_by_group[match(  as.vector(unique(totalInt_data$SampleGroup)), mean_median_by_group$SampleGroup), ]
 
    samplesingroups <- table(totalInt_data$SampleGroup)
    cumulative_sum <- cumsum(c(1,unname(samplesingroups)))
    
    x_values_start <-as.vector( cumulative_sum[-length(cumulative_sum)]-0.5)
    x_values_end <- as.vector( x_values_start+unname(samplesingroups))
    
    
    MedianIntensity_plot_all <- ggplot(totalInt_data, aes(x = RowNames, y = SampleMedian, fill = factor(SampleGroup)))+
      geom_col(position = "dodge") +
      geom_segment(data = mean_median_by_group, aes(x = x_values_start, xend = x_values_end, y = mean_median, yend = mean_median), color = "black", size = 1) +
      labs(title = "Sample medians by group", x = "Sample groups", y = "Sample Medians") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 0.9)) +
      theme(plot.title = element_text(size = 22))
    
    pool_plot_list[["Median Intensity, Group Means"]] <- MedianIntensity_plot_all
    
    if (is.null(Save_as_Plot) == FALSE){ # save plot
      suppressMessages(suppressWarnings(ggsave(filename = paste0(Result_Folder, "/MedianIntensity_GroupMeans",OutputName,".",Save_as_Plot, sep = ""), plot = MedianIntensity_plot_all, width = 12,  height = 10)))
    }
    
    
    ## Signal drift pool samples
    plotdata[plotdata$SampleGroup != Input_SettingsInfo[["PoolSamples"]],3] <- "Sample"
    signaldrift_plot_pools <- ggplot(plotdata, aes(x = RowNames, y = SampleMedian, fill = SampleGroup)) +
      geom_col(position = "dodge") +
      labs(title = "Signal drift, Pool Samples", x = "Samples by acquirement order", y = "Sample Medians") +
      theme_minimal()  + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
      theme(plot.title = element_text(size=22))
    
    signaldrift_plot_pools  <- signaldrift_plot_pools + geom_smooth(data = filter(plotdata, SampleGroup == "Pool"), # add trend line
                                  aes(x = as.numeric(RowNames), y = SampleMedian),
                                  se = FALSE,method = "lm", color = "red")
    pool_plot_list[["Signal drift, Pool Samples"]] <- signaldrift_plot_pools
    
    if (is.null(Save_as_Plot) == FALSE){ # save the plot
      suppressMessages(suppressWarnings(ggsave(filename = paste0(Result_Folder, "/SignalDrift_PoolSamples",OutputName,".",Save_as_Plot, sep = ""), plot = signaldrift_plot_pools, width = 12,  height = 10)))
    }

    ## Signal drift Valine d8 
    if(length(grep("d8", colnames(Drift_data_ordered)))>0){
      plotdataIS <- data.frame("SampleMedian" = Drift_data_ordered %>% select(  grep("d8", colnames(Drift_data_ordered))), SampleGroup = Drift_data_ordered$Conditions) # get the data to plot
      plotdataIS <- rownames_to_column(plotdataIS, var = "RowNames")
      plotdataIS$RowNames <- factor(plotdataIS$RowNames, levels = plotdataIS$RowNames)
      colnames(plotdataIS)[2] <- "ValineIS"
      plotdataIS$Trend <- 1
      
      # calculate Valided8 sample CV 
      ValSampCV <- (sd(plotdataIS$ValineIS, na.rm =T)/  mean(plotdataIS$ValineIS, na.rm =T))*100
      
      # Create a bar plot using ggsignaldrift_plot_valine
      signaldrift_plot_valine <- ggplot(plotdataIS, aes(x = RowNames, y = ValineIS, fill = SampleGroup)) +
        geom_col(position = "dodge") +
        labs(title = "Signal drift, Valine-d8 (IS)", x = "Samples by acquirement order", y = "Peak Area", subtitle =  paste0("PoolCV = ",  round( result_df_final[ grep("aline-d8",  result_df_final$Metabolite ),2],3), ", SampleCV = ",round(ValSampCV,4))) +
        theme_minimal()  + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
        theme(plot.title = element_text(size=22))
      
      signaldrift_plot_valine <- signaldrift_plot_valine + geom_smooth(data = plotdataIS, # add trend line to plot
                                                                       aes(x = as.numeric(RowNames), y = ValineIS, group = Trend),
                                                                       se = FALSE, method = "lm", color = "red")
      
      pool_plot_list[["Signal drift, Valine-d8 (IS)"]] <- signaldrift_plot_valine
      
      if (is.null(Save_as_Plot) == FALSE){ # save the plot
        suppressMessages(suppressWarnings(ggsave(filename = paste0(Result_Folder, "/SignalDrift_Valine-d8(IS)",OutputName,".",Save_as_Plot, sep = ""), plot = signaldrift_plot_valine, width = 12,  height = 10)))
      }
      
      
      #Valine by sample group
      plotdataIS <- plotdataIS[order(plotdataIS$SampleGroup), ]
      plotdataIS$RowNames <- factor(plotdataIS$RowNames, levels = plotdataIS$RowNames)
      
      
      Valine_by_group <- plotdataIS %>% group_by(SampleGroup) %>% dplyr::summarize(mean_median = mean(ValineIS))
      
      Valine_by_group <-  Valine_by_group[match(  as.vector(unique(plotdataIS$SampleGroup)), Valine_by_group$SampleGroup), ]
      
      samplesingroups <- table(plotdataIS$SampleGroup)
      cumulative_sum <- cumsum(c(1,unname(samplesingroups)))
      
      x_values_start <-as.vector( cumulative_sum[-length(cumulative_sum)]-0.5)
      x_values_end <- as.vector( x_values_start+unname(samplesingroups))
      
      
      # Create a bar plot using ggsignaldrift_plot_valine
      signaldrift_plot_valine <- ggplot(plotdataIS, aes(x = RowNames, y = ValineIS, fill = SampleGroup)) +
        geom_col(position = "dodge") +
        geom_segment(data = Valine_by_group, aes(x = x_values_start, xend = x_values_end, y = mean_median, yend = mean_median), color = "black", size = 1)+
        labs(title = "Valine-d8 (IS) Sample group", x = "Samples by Group", y = "Peak Area", subtitle =  paste0("PoolCV = ",  round( result_df_final[ grep("aline-d8",  result_df_final$Metabolite ),2],3), ", SampleCV = ",round(ValSampCV,4))) +
        theme_minimal()  + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
        theme(plot.title = element_text(size=22))
      
      pool_plot_list[["SampleGroup, Valine-d8 (IS)"]] <- signaldrift_plot_valine
      
      if (is.null(Save_as_Plot) == FALSE){ # save the plot
        suppressMessages(suppressWarnings(ggsave(filename = paste0(Result_Folder, "/SampleGroup_Valine-d8(IS)",OutputName,".",Save_as_Plot, sep = ""), plot = signaldrift_plot_valine, width = 12,  height = 10)))
      }

    }
  }
  
  
  # Prepare output list
  DF_list <- list("Filtered_Input_data" = filtered_Input_data, "CV_result" = result_df_final)
  Pool_Estimation_res <- list("DF" = DF_list,"Plot" = pool_plot_list)
  # Save Rdata
  save(Pool_Estimation_res, file = file.path(RDS_folder, paste( "PoolEstimationRes", OutputName,".Rdata", sep = "")))
  
  }
  invisible(return(Pool_Estimation_res))
}

