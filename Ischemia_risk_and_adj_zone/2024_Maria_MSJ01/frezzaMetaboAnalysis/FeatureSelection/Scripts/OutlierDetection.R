OutlierDetection <- function(Result_Folder ,
                             RDS_folder ,
                             scriptpath,
                             FolderName = NULL,
                             OutputName = NULL,
                             ForceRun = FALSE,
                             
                             data_norm,
                             Input_SettingsFile,
                             Input_SettingsInfo,
                             
                             OutlierLoop = 2,
                             npcs=NULL,
                             HotellinsConfidence = 0.99,
                             
                             Save_as_Results = "csv",
                             Save_as_Plot = "svg"
                             ){
  
  
  
  # Sources
  # usethis::edit_file(file.path(scriptpath,"Functions/PCA.R"))   
  source(file.path(scriptpath,"PCA.R"))
  # usethis::edit_file(file.path(scriptpath,"Functions/ZeroVarianceCheck.R"))   
  source(file.path(scriptpath,"ZeroVarianceCheck.R"))
  
  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)) {OutputName = ""} else{OutputName = paste("_",OutputName, sep="")}
  if (file.exists(file.path(RDS_folder, paste("OutlierDetectionRes", OutputName,".Rdata", sep = ""))) == TRUE & ForceRun == FALSE) { # Check if RDS exists and forceRun is false then skip script and load result from RDS, otherwise run it
    message("Loading already existing ", paste("OutlierDetectionRes", OutputName,".Rdata", sep = ""))
    load(file.path(RDS_folder, paste("OutlierDetectionRes", OutputName,".Rdata", sep = "")), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running Outlier Detection script main body.")
    
    # # # ignore pool samples
    pool_data <- data_norm %>% filter(Input_SettingsFile[,"Conditions"] == Input_SettingsInfo[["PoolSamples"]])
    pool_SettingsFile <- Input_SettingsFile %>% filter(Input_SettingsFile[,"Conditions"] == Input_SettingsInfo[["PoolSamples"]])
    data_norm <- data_norm %>% filter(!Input_SettingsFile[,"Conditions"] == Input_SettingsInfo[["PoolSamples"]])
    Input_SettingsFile <- Input_SettingsFile %>% filter(!Input_SettingsFile[,"Conditions"] == Input_SettingsInfo[["PoolSamples"]])
     
    initial_data = data_norm
    Conditions <- Input_SettingsFile[,Input_SettingsInfo[["Conditions"]]]
    #message("Identification of outlier samples is performed using Hotellin's T2 test to define sample outliers in a mathematical way (Confidence = 0.99 ~ p.val < 0.01) REF: Hotelling, H. (1931), Annals of Mathematical Statistics. 2 (3), 360–378, doi:https://doi.org/10.1214/aoms/1177732979.")
    Outlier_filtering_loop = OutlierLoop
    
    sample_outliers <- list()
    scree_plot_list <- list()
    outlier_plot_list <- list()
    metabolite_zero_var_total_list <- list()
    #zero_var_metab_warning = FALSE
    for (loop in 1:Outlier_filtering_loop){ # loop = 2 here we do 10 rounds of hotelling filtering
      
      ####----------------- Zero variance check for features --------------####
      ZeroVarCheckRes <- ZeroVarCheck(Result_Folder,
                                      Input_data= data_norm, 
                                      OutputName = paste0("Hotelling loop ",loop, OutputName),
                                      RemoveFeatures = TRUE)
      
      
      data_norm <-  ZeroVarCheckRes$Input_data_filtered
      
      ### ### PCA  ### ###
      PCA.res <- prcomp(data_norm, center =  TRUE, scale. =  TRUE)
      outlier_PCA_data <- data_norm
      outlier_PCA_data <- merge(outlier_PCA_data, Input_SettingsFile %>% select(Conditions),by=0, sort = FALSE)
      rownames(outlier_PCA_data)<- outlier_PCA_data$Row.names
      outlier_PCA_data$Row.names <- NULL
      
      ###=== Add a loop for more component plots? ###
      dev.new()
      pca_outlier <- suppressWarnings(PCA(Input_data=data_norm, Plot_SettingsInfo= c(color="Conditions"),
                                             Plot_SettingsFile= outlier_PCA_data, OutputPlotName = paste("PCA outlier test filtering round ",loop),
                                             Save_as_Plot =  NULL))
      
      suppressWarnings(plot(pca_outlier))
      suppressWarnings(outlier_plot_list[[paste("PCA_round",loop,sep="")]] <- recordPlot())
      dev.off()
      ### ### Scree plot ### ### what is the assumption for the knee (% of variance expained)
      # get Scree plot values for inflection point calculation
      inflect_df<- as.data.frame(c(1:length(PCA.res$sdev)))
      colnames(inflect_df)<- "x"
      inflect_df$y<- summary(PCA.res)$importance[2,]
      inflect_df$Cumulative<- summary(PCA.res)$importance[3,]
      #make cumulative variation labels for plot
      screeplot_cumul<- format(round(inflect_df$Cumulative[1:20]*100, 1), nsmall = 1)
      # Calculate the knee and select optimal number of components
      knee=inflection::uik(inflect_df$x,inflect_df$y)
      if (is.null(npcs)){
        npcs = knee -1 #Note: we subtract 1 components from the knee cause the root of the knee is the PC that does not add something. npcs=30
      }
      
      # Make a scree plot with the selected component cut-off for HotellingT2 test
      screeplot <- factoextra::fviz_screeplot(PCA.res, main = paste("PCA Explained variance plot filtering round ",loop, sep = ""),
                                              addlabels = TRUE,
                                              ncp = 20,
                                              geom = c("bar", "line"),
                                              barfill = "grey",
                                              barcolor = "grey",
                                              linecolor = "black",linetype = 1) + theme_classic()+ geom_vline(xintercept = npcs+0.5, linetype = 2, color = "red") +
        annotate("text", x = c(1:20),y = -0.8,label = screeplot_cumul,col = "black", size = 3)
      
      dev.new()
      plot(screeplot)
      outlier_plot_list[[paste("ScreePlot_round",loop,sep="")]] <- recordPlot() # save plot
      dev.off()
      
      ### ### HotellingT2 test for outliers ### ###
      data_hot <- as.matrix(PCA.res$x[,1:npcs])
      message("***Checking for outliers***")
      hotelling_qcc <- qcc::mqcc(data_hot, type = "T2.single",labels = rownames(data_hot),confidence.level = HotellinsConfidence, title = paste("Outlier filtering via HotellingT2 test filtering round ",loop,", with ",HotellinsConfidence, "% Confidence",  sep = ""), plot = FALSE)
      HotellingT2plot_data <- as.data.frame(hotelling_qcc$statistics)
      HotellingT2plot_data <- rownames_to_column(HotellingT2plot_data, "Samples")
      colnames(HotellingT2plot_data) <- c("Samples", "Group summary statisctics")
      
      outlier <- HotellingT2plot_data %>% filter(HotellingT2plot_data$`Group summary statisctics`>hotelling_qcc$limits[2])
      limits <- as.data.frame(hotelling_qcc$limits)
      legend <- colnames(HotellingT2plot_data[2])
      LegendTitle = "Limits"
      
      HotellingT2plot <- ggplot(HotellingT2plot_data, aes(x = Samples, y = `Group summary statisctics`, group = 1, fill = ))
      HotellingT2plot <- HotellingT2plot +
        geom_point(aes(x = Samples,y = `Group summary statisctics`), color = 'blue', size = 2) +
        geom_point(data = outlier, aes(x = Samples,y = `Group summary statisctics`), color = 'red',size = 3) +
        geom_line(linetype = 2)
      
      #draw the horizontal lines corresponding to the LCL,UCL
      HotellingT2plot <- HotellingT2plot + geom_hline(aes(yintercept = limits[,1]), color = "black", data = limits,  show.legend = F) +
        geom_hline(aes(yintercept = limits[,2], linetype = "UCL"), color = "red", data = limits, show.legend = T) +
        #only the LCl and UCL to be shown in y axis
        scale_y_continuous(breaks = sort(c(ggplot_build(HotellingT2plot)$layout$panel_ranges[[1]]$y.major_source, c(limits[,1],limits[,2]))))
      
      HotellingT2plot <- HotellingT2plot + theme_classic()
      HotellingT2plot <- HotellingT2plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
      HotellingT2plot <- HotellingT2plot + ggtitle(paste("Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"))
      HotellingT2plot <- HotellingT2plot + scale_linetype_discrete(name = LegendTitle,)
      HotellingT2plot <- HotellingT2plot + theme(plot.title = element_text(size = 13))+#, face = "bold")) +
        theme(axis.text = element_text(size = 12))
      
      dev.new()
      plot(HotellingT2plot)
      outlier_plot_list[[paste("HotellingsPlot_round",loop,sep="")]] <- recordPlot()
      dev.off()
      
      
      
      ### Save the outlier detection plots in the outlier detection folder
      suppressWarnings(ggsave(filename = paste(Result_Folder, "/PCA_OD_round_" ,loop , OutputName,".", Save_as_Plot, sep = ""),
                              plot = pca_outlier, width = 10,height = 8))
      ggsave(filename = paste(Result_Folder, "//Scree_plot_OD_round_" ,loop , OutputName,".",Save_as_Plot, sep = ""),
             plot = screeplot, width = 10,height = 8)
      ggsave(filename = paste(Result_Folder, "/Hotelling_OD_round_" ,loop , OutputName,".", Save_as_Plot, sep = ""),
             plot = HotellingT2plot, width = 10,height = 8)
      
      # Here for the outliers we use confidence of 0.999 and p.val < 0.01.
      if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 0){ # loop for outliers until no outlier is detected
        data_norm <- data_norm  # filter the selected outliers from the data
        break
      } else if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 1){
        data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]
        Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]
        
        # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
        hotelling_qcc[["violations"]][["beyond.limits"]][1] <-   rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][1]]
        sample_outliers[loop] <- list(hotelling_qcc[["violations"]][["beyond.limits"]])
        
      } else {
        data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]
        Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]
        
        # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
        sm_out <- c() # list of outliers samples
        for (i in 1:length(hotelling_qcc[["violations"]][["beyond.limits"]])){
          sm_out <-  append(sm_out, rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][i]]  )
        }
        sample_outliers[loop] <- list(sm_out )
      }
    }
    
    # Save outlier result
    pdf(file = paste(Result_Folder, "/Outlier_testing", OutputName,".pdf", sep = ""), onefile = TRUE ) # or Outlier detection related plots
    suppressWarnings(for (plot in outlier_plot_list) {
      replayPlot(plot)
    })
    dev.off()
    
    # Print outlier result
    if (length(sample_outliers)>0){
    outlierSamples <- vector(mode = "list", length  = length(sample_outliers))
    namesoutlierSamples <- paste0("Filtering Round: ", 1:length(sample_outliers))
    outlierSamples <- setNames(outlierSamples, namesoutlierSamples)}
    if(length(sample_outliers)>0){
      for (i in 1:length(sample_outliers)  ){
        message('Sample outliers were identified:')
        message("Filtering round ",i ,", Outlier Samples:", paste( head(sample_outliers[[i]]) ," "))
        outlierSamples[i] <- as.data.frame(paste(paste0( head(sample_outliers[[i]]))), col.names = NULL, row.names	= NULL)
      }
      #export the outliers in a text format
      
      sink(file.path(Result_Folder, "Outliers.txt"))
      print(outlierSamples)
      sink()
    }else{  message("No sample outliers were found.")}
    
    
    # Make putput table
    total_outliers <- hash::hash() # make a dictionary
    if(length(sample_outliers) > 0){ # Create columns with outliers to merge to output dataframe
      for (i in 1:length(sample_outliers)  ){
        total_outliers[[paste("Outlier_filtering_round_",i, sep = "")]] <- sample_outliers[i]
      }
    }
    
    data_norm_filtered_full <- as.data.frame(initial_data)
    
    if(length(total_outliers) > 0){  # add outlier information to the full output dataframe
      data_norm_filtered_full$Outliers <- "no"
      for (i in 1:length(total_outliers)){
        for (k in 1:length( hash::values(total_outliers)[i] ) ){
          data_norm_filtered_full[as.character(hash::values(total_outliers)[[i]]) , "Outliers"] <- hash::keys(total_outliers)[i]
        }
      }
    }else{
      data_norm_filtered_full$Outliers <- "no"
    }
    
    data_norm_filtered_full <- data_norm_filtered_full %>% relocate(Outliers) #Put Outlier columns in the front
    data_norm_filtered_full <- merge(Input_SettingsFile, data_norm_filtered_full,  by = 0, sort = FALSE) # add the design in the output df (merge by rownames/sample names)
    rownames(data_norm_filtered_full) <- data_norm_filtered_full$Row.names
    data_norm_filtered_full$Row.names <- c()
    
    dtp <- data_norm_filtered_full %>%
      select(Conditions, Outliers) %>%
      mutate(Outliers = case_when(Outliers == "no" ~ 'no',
                                  Outliers == "Outlier_filtering_round_1" ~ ' Outlier_filtering_round = 1',
                                  Outliers == "Outlier_filtering_round_2" ~ ' Outlier_filtering_round = 2',
                                  Outliers == "Outlier_filtering_round_3" ~ ' Outlier_filtering_round = 3',
                                  Outliers == "Outlier_filtering_round_4" ~ ' Outlier_filtering_round = 4',
                                  TRUE ~ 'Outlier_filtering_round = or > 5'))
    dtp$Outliers <- relevel( as.factor(dtp$Outliers), ref="no")
    
    dev.new()
    suppressWarnings(pca_QC <-PCA(Input_data=as.data.frame(initial_data), Plot_SettingsInfo= c(color="Conditions", shape = "Outliers"),
                                     Plot_SettingsFile= dtp,OutputPlotName = "Quality Control PCA Condition clustering and outlier check",
                                     Save_as_Plot =  NULL))
    dev.off()
    
    suppressWarnings(ggsave(filename = paste0(Result_Folder, "/HotellingPCA", OutputName,".",Save_as_Plot), plot = pca_QC, width = 10,  height = 8))
    
    
    
  # Save results
  #########################################################
  ### ### ###  Make list with output dataframes ### ### ###
 # data_norm_filtered_full <- data_norm_filtered_full[-grep("pool",tolower(data_norm_filtered_full$Conditions), ignore.case = T),]
    
    
  preprocessing_output_list <- list("Processed_data" = data_norm_filtered_full)
  ##Write to file
  preprocessing_output_list_out <- lapply(preprocessing_output_list, function(x) rownames_to_column(x, "Sample_ID")) #  # use this line to make a sample_ID column in each dataframe
  # save result
  writexl::write_xlsx(preprocessing_output_list_out, paste0(Result_Folder, "/ProcessedData", OutputName,".xlsx", sep = ""))#,showNA = TRUE)
  # save result csv for Binner
  write.csv(preprocessing_output_list_out,paste0(Result_Folder, "/ProcessedData", OutputName,".csv", sep = ""), row.names = FALSE)
  
 
  # Make Result list
  OutlierDetectionRes = list("DFs" = preprocessing_output_list, "Plots" = list("HotellingPCA" =  pca_QC))
  # Save output RDS
  save(OutlierDetectionRes , file = file.path(RDS_folder,  paste("OutlierDetectionRes", OutputName,".Rdata", sep = "")))
  }
  return(OutlierDetectionRes)
}