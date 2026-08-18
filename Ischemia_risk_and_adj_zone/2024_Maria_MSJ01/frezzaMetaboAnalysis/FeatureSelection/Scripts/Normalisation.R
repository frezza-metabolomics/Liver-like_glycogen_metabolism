
Normalisation <- function(Result_Folder,
                          RDS_folder ,
                          scriptpath,
                          FolderName = NULL,
                          OutputName = NULL,
                          ForceRun = TRUE,
                          CoRe=FALSE,
                          
                          NA_removed_matrix,
                          Input_SettingsFile,
                          Input_SettingsInfo,
                          
                          Method = "PQN", 
                          PQN_selected_persentage=100,
                          
                          Save_as_Plot = "svg"
                          ){
  
  
  # Load dependencies
  # usethis::edit_file(file.path(scriptpath,"Functions/PCA.R"))   
  source(file.path(scriptpath,"PCA.R"))
  

  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)) {OutputName = ""} else{OutputName = paste("_",OutputName, sep="")}
  if (file.exists(file.path(RDS_folder, paste("NormalisationRes", OutputName,".Rdata", sep = ""))) == TRUE & ForceRun == FALSE) { # Check if RDS exists and forceRun is false then skip script and load result from RDS, otherwise run it
    message("Loading already existing ", paste("NormalisationRes", OutputName,".Rdata", sep = ""))
    load(file.path(RDS_folder, paste("NormalisationRes", OutputName,".Rdata", sep = "")), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running Normalisation script main body.")
    Norm_data <- NA_removed_matrix 
    Norm_SettingsFile <- Input_SettingsFile

    
    if (Method == "TIC"){ #Normalise for total ion counts
      message("Total Ion Count (TIC) normalization is used to reduce the variation from non-biological sources, while maintaining the biological variation. REF: Wulff et. al., (2018), Advances in Bioscience and Biotechnology, 9, 339-351, doi:https://doi.org/10.4236/abb.2018.98022")
      RowSums <- rowSums(Norm_data)
      Median_RowSums <- median(RowSums) #This will built the median
      Data_TIC_Pre <- apply(Norm_data, 2, function(i) i/RowSums) #This is dividing the ion intensity by the total ion count
      Data_TIC <- Data_TIC_Pre*Median_RowSums #Multiplies with the median metabolite intensity
      Data_TIC<- as.data.frame(Data_TIC)
      
    }else if (Method == "Log2"){ 
      Data_TIC<- as.data.frame(log2(Norm_data+1))
    }else if (Method == "PQN"){ 
      temp <- Norm_data
      print(unbal_reg(t(temp)))
      data_raw_no_zero <- temp
      data_raw_no_zero[data_raw_no_zero==0] <- 1
      #if (data_raw_no_zero[data_raw_no_zero==0] ) warning("matrix includes 0s, need to be imputed or properly treated")
      resamp_mswsd(metabo.data = t(data_raw_no_zero),Results_folder_Preprocessing_folder=Result_Folder,OutputName= OutputName)
      
      # Define ANSI escape code for blue text
      blue_text <- "\033[34m"
      reset_text <- "\033[0m"
      
      # Display the prompt in blue using cat
      cat(blue_text, "Please enter a percentage for PQN normalization.\n",
          "Check the MSWSD plot here: ", paste(Result_Folder, "/resamp_mswsd", OutputName, ".pdf", sep=""), "\n",
          "The optimal percentage is the maximum percentage where MSWSD approaches a nearly constant minimal value and its variation does not decrease further.\n",
          "The percenrage must be numeric between 10 and 100", reset_text, sep = "")
      
      # Use readline to get the user input
      selected_percentage <- readline(prompt = "Enter percentage here: ")
      #  selected_percentage = 15
      selected_percentage <- as.numeric(selected_percentage)
      
      suppressMessages(  data_norm_PQN <- norm_unbal(t(data_raw_no_zero),selected_percentage,"PQN",Results_folder_Preprocessing_folder=Result_Folder)
      )
      
      Data_TIC <- as.data.frame(t(data_norm_PQN))
    }
    
  
    if (CoRe ==  TRUE){
      # # If we have more than 1 media group
      # if(length(unique(Input_SettingsFile[ grepl(Input_SettingsInfo["CoRe_media"],Input_SettingsFile$Conditions),"Conditions"]))>1 ){
      # for each media condition
      CoReNormList <- list()
      for(mediaCond in  unique(Input_SettingsFile[ grepl(Input_SettingsInfo["CoRe_media"],Input_SettingsFile$Conditions),"Conditions"])){           # mediaCond = unique(Input_SettingsFile[ grepl("BasalMedia",Input_SettingsFile$Conditions),"Conditions"])[1]
        # Subset the Data
        CoReSubGroupData <- Data_TIC[grep(strsplit(mediaCond, "_")[[1]][2], Input_SettingsFile$Conditions ),]
        CoReSubGroupDesign <-  Input_SettingsFile[grep(strsplit(mediaCond, "_")[[1]][2], Input_SettingsFile$Conditions ),]
        
        # Get media samples
        CoRe_medias <-  CoReSubGroupData[grep(mediaCond, CoReSubGroupDesign$Conditions),] %>% as.data.frame()
        # IF we have only 1 media sample
        if(dim(CoRe_medias)[1]==1){
          warning("Only 1 CoRe_media sample was found. Thus, the consistency of the CoRe_media samples cannot be checked. It is assumed that the CoRe_media samples are already summed.")
          CoRe_media_df <- CoRe_medias %>% t() %>% as.data.frame()
          colnames(CoRe_medias) <- "CoRe_mediaMeans"
        }else{
          # If we have more than 1 media sample
          
          ## Media_control PCA
          media_pca_data <- merge(CoReSubGroupDesign %>% select(Conditions), CoReSubGroupData, by=0) %>%
            column_to_rownames("Row.names") %>%
            mutate(Sample_type = case_when(Conditions == mediaCond~ "CoRe_media",
                                           TRUE ~ "Sample"))
          
          ####----------------- Zero variance check for features --------------####
          # usethis::edit_file(file.path(scriptpath,"ZeroVarianceCheck.R"))   
          source(file.path(scriptpath,"ZeroVarianceCheck.R"))
          ZeroVarCheckRes <- ZeroVarCheck(Result_Folder = Result_Folder,
                                          Input_data= media_pca_data, 
                                          OutputName = OutputName,
                                          RemoveFeatures = TRUE)
          dev.new()
          pca_QC_media <-invisible(PCA(Input_data=ZeroVarCheckRes$Input_data_filtered %>%select(-Conditions, -Sample_type), Plot_SettingsInfo= c(color="Sample_type"),
                                       Plot_SettingsFile= media_pca_data, OutputPlotName = "QC Media_samples",
                                       Save_as_Plot =  NULL))
          dev.off()
          qc_plot_list[["CoRe PCA MediaSamples"]] <- pca_QC_media
          
          
          ## Check metabolite variance
          # Thresholds
          Threshold_cv = 30
          data_cv <- CoRe_medias
          
          ## Coefficient of Variation
          result_df <- apply(data_cv, 2,   function(x) { (sd(x, na.rm =T)/  mean(x, na.rm =T))*100 } ) %>% t()%>% as.data.frame()
          result_df[1, is.na(result_df[1,])]<- 0
          rownames(result_df)[1] <- "CV"
          
          result_df <- result_df %>% t()%>%as.data.frame() %>% rowwise() %>%
            mutate(HighVar = CV > Threshold_cv) %>% as.data.frame()
          rownames(result_df)<- colnames(data_cv)
          
          # calculate the NAs
          NAvector <- apply(data_cv, 2,  function(x) { (sum(is.na(x))/length(x))*100 }  )#%>% t()%>% as.data.frame()
          result_df$MissingValuePercentage <- NAvector
          
          cv_result_df <- result_df
          
          HighVar_metabs <- sum(result_df$HighVar == TRUE)
          if(HighVar_metabs>0){
            message(paste0(HighVar_metabs, " of variables have high variability in the CoRe_media samples. Consider checking the pooled samples to decide whether to remove these metabolites or not."))
          }
          ###########################
          #Make histogram of CVs
          HistCV <- invisible(ggplot(cv_result_df, aes(CV)) +
                                geom_histogram(aes(y=after_stat(density)), color="black", fill="white")+
                                geom_vline(aes(xintercept=Threshold_cv),
                                           color="darkred", linetype="dashed", size=1)+
                                geom_density(alpha=.2, fill="#FF6666") +
                                labs(title=paste0("CV for metabolites of ",mediaCond," samples"),x="Coefficient of variation (CV)", y = "Frequency")+
                                theme_classic())
          
          qc_plot_list[["CoRe_Media_CV_Hist"]] <- HistCV
          
          # Do sample outlier testing
          if(dim(CoRe_medias)[1]>=3){
            
            Outlier_data <- CoRe_medias
            Outlier_data <- Outlier_data %>% mutate_all(.funs = ~ FALSE)
            
            while(HighVar_metabs>0){
              
              #remove the furthest value from the mean
              if(HighVar_metabs>1){
                max_var_pos <-  data_cv[,result_df$HighVar == TRUE]  %>%
                  as.data.frame() %>%
                  mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
                  summarise_all(.funs = ~ which.max(abs(.)))
              }else{
                max_var_pos <-  data_cv[,result_df$HighVar == TRUE]  %>%
                  as.data.frame() %>%
                  mutate_all(.funs = ~ . - mean(., na.rm = TRUE)) %>%
                  summarise_all(.funs = ~ which.max(abs(.)))
                colnames(max_var_pos)<- colnames( data_cv)[result_df$HighVar == TRUE]
                
              }
              
              # Remove rows based on positions
              for(i in 1:length(max_var_pos)){
                data_cv[max_var_pos[[i]],names(max_var_pos)[i]] <- NA
                Outlier_data[max_var_pos[[i]],names(max_var_pos)[i]] <- TRUE
              }
              
              # ReCalculate coefficient of variation for each column in the filtered data
              result_df <- apply(data_cv, 2,   function(x) { sd(x, na.rm =T)/  mean(x, na.rm =T) } ) %>% t()%>% as.data.frame()
              result_df[1, is.na(result_df[1,])]<- 0
              rownames(result_df)[1] <- "CV"
              
              result_df <- result_df %>% t()%>%as.data.frame() %>% rowwise() %>%
                mutate(HighVar = CV > Threshold_cv) %>% as.data.frame()
              rownames(result_df)<- colnames(data_cv)
              
              HighVar_metabs <- sum(result_df$HighVar == TRUE)
            }
            
            
            data_cont <- Outlier_data %>% t() %>% as.data.frame()
            
            # List to store results
            fisher_test_results <- list()
            large_contingency_table <- matrix(0, nrow = 2, ncol = ncol(data_cont))
            
            for (i in 1:length(colnames(data_cont))) {
              sample = colnames(data_cont)[i]
              current_sample <- data_cont[, sample]
              
              contingency_table <- matrix(0, nrow = 2, ncol = 2)
              contingency_table[1, 1] <- sum(current_sample)
              contingency_table[2, 1] <- sum(!current_sample)
              contingency_table[1, 2] <- sum(rowSums(data_cont) - current_sample)
              contingency_table[2, 2] <- dim(data_cont %>% select(!all_of(sample)))[1]*dim(data_cont %>% select(!all_of(sample)))[2] -sum( rowSums(data_cont) - current_sample)
              
              # Fisher's exact test
              fisher_test_result <- fisher.test(contingency_table)
              fisher_test_results[[sample]] <- fisher_test_result
              
              # Calculate the sum of "TRUE" and "FALSE" for the current sample
              large_contingency_table[1, i] <- sum(current_sample)  # Sum of "TRUE"
              large_contingency_table[2, i] <- sum(!current_sample) # Sum of "FALSE"
              
            }
            
            # Convert the matrix into a data_contframe for better readability
            contingency_data_contframe <- as.data.frame(large_contingency_table)
            colnames(contingency_data_contframe) <- colnames(data_cont)
            rownames(contingency_data_contframe) <- c("HighVar", "Low_var")
            
            contingency_data_contframe <- contingency_data_contframe %>% mutate(Total = rowSums(contingency_data_contframe))
            contingency_data_contframe <- rbind(contingency_data_contframe, Total= colSums(contingency_data_contframe))
            
            
            different_samples <- c()
            for (sample in colnames(data_cont)) {
              p_value <- fisher_test_results[[sample]]$p.value
              if (p_value < 0.05) {  # Adjust the significance level as needed
                different_samples <- c(different_samples, sample)
              }
            }
            
            if(is.null(different_samples)==FALSE){
              warning("The CoRe_media samples ", paste(different_samples, collapse = ", "), " were found to be different from the rest. They will not be included in the sum of the CoRe_media samples.")
            }
            # Filter the CoRe_media samples
            CoRe_medias <- CoRe_medias %>% filter(!rownames(CoRe_medias) %in% different_samples)
          } # end of media outlier testing
          CoRe_media_df <- as.data.frame(data.frame("CoRe_mediaMeans"=  colMeans( CoRe_medias, na.rm = TRUE)))
          
          
          #cv_result_df <- rownames_to_column(cv_result_df, "Metabolite")
          write.table(cv_result_df,row.names =  FALSE, file = paste(Result_Folder,"/CV_table_CoRe",".csv",sep =  ""))
          if(dim(CoRe_medias)[1]>3){
            write.table(contingency_data_contframe,row.names =  TRUE, file = paste(Result_Folder,"/Contigency_table_CoRe",".csv",sep =  ""))
          }
        } # end of if we have more than 1 media
        
        # Remove CoRe_media samples from the data
        CoReSubGroupData <- CoReSubGroupData[-c(which(CoReSubGroupDesign$Conditions == mediaCond)),]
        
        Data_TIC_CoRe <- as.data.frame(t( apply(t(CoReSubGroupData),2, function(i) i-CoRe_media_df$CoRe_mediaMeans)))  #Subtract from each sample the CoRe_media mean
        
        CoReNormList[[mediaCond]] <- Data_TIC_CoRe
      }# end for each media cond
      # Bind all dataframes
      Data_TIC_CoRe <- bind_rows(CoReNormList)
      
      
      # end of more than 2 media conditions
      # }else{# If we have only 1 condition
      # }
      
      
      
      Input_SettingsFile <- Input_SettingsFile[rownames(Input_SettingsFile) %in% rownames(Data_TIC_CoRe),]
      Data_TIC_CoRe <- Data_TIC_CoRe[rownames(Input_SettingsFile), ]
      
      
      
      message("CoRe data are normalised using CoRe_norm_factor")
      Data_TIC <- apply(Data_TIC_CoRe, 2, function(i) i*Input_SettingsFile$NormVector) %>% as.data.frame()
      rownames(Data_TIC) <- rownames(Data_TIC_CoRe)
      
      if (var(Input_SettingsFile$NormVector) ==  0){
        warning("The growth rate or growth factor for normalising the CoRe result, is the same for all samples")
      }
      # # Remove CoRe_media samples from the data
      # Input_SettingsFile <- Input_SettingsFile[Input_SettingsFile$Conditions!=mediaCond,]
      
    }# end core
    
    # Start QC plot list
    qc_plot_list <- list()
    
    ## PLots
    ### Make normalization plot
    Norm_SettingsFile$Var2 <- rownames(Norm_SettingsFile)
    rawNormplotdata <- Norm_data
    rawNormplotdata[rawNormplotdata==0]<- 1
    # For Raw data
    YlinE <- log2(median(Biobase::rowMedians(as.matrix(rawNormplotdata))))
    meansOF <- log2(rowMeans(as.matrix(rawNormplotdata)))
    mdat <- reshape2::melt(t(log(as.matrix(rawNormplotdata),2)))  ## convert to long format
    mdat <- merge(mdat,Norm_SettingsFile %>% select(Conditions, Var2), by="Var2" )
    p1 <- ggplot2::ggplot(mdat,ggplot2::aes(x=factor(Var2),y=value, fill = Conditions))+
      ggplot2::geom_violin(size = 0.1)+
      ggplot2::geom_boxplot(ggplot2::aes(x = factor(Var2), y = value),outlier.colour = "darkred",
                            outlier.size = 0.9,width = 0.2, size = 0.3)+
      ggplot2::guides(fill=ggplot2::guide_legend(title="samples")) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = YlinE), size = 0.2) +
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.9), angle=45, hjust = 1, vjust = 1)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
      ggplot2::xlab("Samples") +
      ggplot2::ylab("log2 of Raw Intensity ")
    
    # For normalized data
    YlinE <- log2(median(Biobase::rowMedians(as.matrix(Data_TIC))))
    meansOF <- log2(rowMeans(as.matrix(Data_TIC)))
    mdat <- reshape2::melt(t(log(as.matrix(Data_TIC),2)))  ## convert to long format
    mdat <- merge(mdat,Norm_SettingsFile %>% select(Conditions, Var2), by="Var2" )
    p2 <- ggplot2::ggplot(mdat,ggplot2::aes(x=factor(Var2),y=value, fill = Conditions))+
      ggplot2::geom_violin(size = 0.1)+
      ggplot2::geom_boxplot(ggplot2::aes(x = factor(Var2), y = value),outlier.colour = "darkred",
                            outlier.size = 0.9,width = 0.2, size = 0.3)+
      ggplot2::guides(fill=ggplot2::guide_legend(title="samples")) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = YlinE), size = 0.2) +
      ggplot2::theme_classic()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.9), angle=45, hjust = 1, vjust = 1)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold")) +
      ggplot2::xlab("Samples") +
      ggplot2::ylab("log2 of Normalized Intensity ")
    #dev.off()
    
    library("patchwork")
    Normalization_plot <- p1 + p2
    Normalization_plot <- patchwork::wrap_plots(p1, p2, ncol = 1)
    Normalization_plot
    
    suppressMessages(suppressWarnings(ggsave(filename = paste0(Result_Folder, "/Normalization_plot", OutputName,".",Save_as_Plot), plot = Normalization_plot, width = 14,  height = 8)))
    
    ### Make PCA plot
    #Raw data
    # rawNormplotdata
    # Norm_SettingsFile
    
    pca_data <- merge(Norm_SettingsFile %>% select(Conditions), rawNormplotdata, by = 0, sort=FALSE) %>% column_to_rownames("Row.names")
    
    #which(sapply(pca_data, sd, na.rm=T)==0)
    
    
    pca_plot1 <- suppressWarnings(PCA(Input_data=pca_data %>%
                                        select(-Conditions), Plot_SettingsInfo= c(color="Conditions"),
                                      Plot_SettingsFile = pca_data, OutputPlotName = "PCA raw imputed data", Save_as_Plot =  NULL))
    
    # Norm data
    pca_data <- merge(Norm_SettingsFile %>% select(Conditions), Data_TIC, by = 0, sort=FALSE) %>% column_to_rownames("Row.names")
    pca_plot2 <- suppressWarnings(PCA(Input_data=pca_data %>%
                                        select(-Conditions), Plot_SettingsInfo= c(color="Conditions"),
                                      Plot_SettingsFile = pca_data, OutputPlotName = "PCA normalised imputed data", Save_as_Plot =  NULL))
    
    suppressWarnings(PCA_norm_plot <- pca_plot1 + pca_plot2)
    PCA_norm_plot <- patchwork::wrap_plots(pca_plot1, pca_plot2, ncol = 2)
    #PCA_norm_plot
    
    suppressMessages(suppressWarnings(ggsave(filename = paste0(Result_Folder, "/normalisationPCA", OutputName,".",Save_as_Plot), plot = PCA_norm_plot, width = 20,  height = 8)))

    
    NormalisationRes <- list("DF" = list("NormalizedData"= Data_TIC),"Plot" =  list("Normalization"= Normalization_plot))
    # Save output RDS
    save(NormalisationRes, file = file.path(RDS_folder,  paste("NormalisationRes", OutputName,".Rdata", sep = "")))

  }
  return(NormalisationRes)
}




###### ---------------------- The code below is the mswd_resamp_bubli.R script for PQN -------------------------------- #######
suppressMessages(library(multtest))

# Warranty
# THE AUTHORS MAKE NO WARRANTIES, EXPRESSED OR IMPLIED, REGARDING THE FITNESS OF OUR SOFTWARE FOR ANY PARTICULAR PURPOSE. THE AUTHORS CLAIM NO
# LIABILITY FOR DATA LOSS OR OTHER PROBLEMS CAUSED DIRECTLY OR INDIRECTLY BY THE SOFTWARE. THE USER IS ASSUMING THE ENTIRE RISK AS TO THE SOFTWARES
# QUALITY AND ACCURACY. 
# Parameters
# Note features should be in rows, samples in columns
mswsd <- function(x.data, n)
{
  use_percent = n
  #cat("#### percentage of used features ######","\n")
  #cat("percentage of used features", use_percent,"\n")
  # Set input data
  data <- x.data
  # determination of invariant features
  rem<-floor(nrow(data)-(nrow(data)*use_percent/100))
  # calculate variance
  vari<-apply(data,1,var) #calculate variance for each feature across samples
  # remove features with highest variance
  sel<-order(vari)[1:(nrow(data)-rem)]
  data2<-data[sel,] #data matrix that contains only features with low variance
  rm(vari)
  # Now we want to compute whether enough features have been removed
  # For this we calculate for each spectrum feature wise scaling factors to a baseline spectrum
  # In case that only non-regulated features remain the standard deviation of these scaling factors should
  # be minimal. Note this is done for each spectrum separately since the scaling factors could be
  # drastically different for different spectra. To obtain a value for the complete data set we 
  # calculate the mean over all spectra
  linear.baseline <- apply(abs(data2),1,median) #compute baseline
  # compute for each feature in each spectrum a scaling factor to the baseline
  fac_feat<-matrix(nrow=nrow(data2), ncol=ncol(data2))
  for  (i in 1:ncol(data2)) #samples
  {
    for (j in 1:nrow(data2)) # features
    {
      fac_feat[j,i]<-data2[j,i]/linear.baseline[j]
    }
  }
  # standard deviation of factors witin sample shhould be ideally small
  sd_fac_feat<-numeric(ncol(data2))
  for (i in 1:ncol(data2))
  {
    sd_fac_feat[i]<-sd(fac_feat[,i])
  }
  # mean standard deviation of factors across samples
  mean_sd_fac_feat=mean(sd_fac_feat)
  #cat("##### standard deviation across feature wise correction factors ####", "\n")
  #cat("mean standard deviation across correction factors", mean_sd_fac_feat,"\n")
  #cat("##### scaling factors ######","\n")
  
  return(mean_sd_fac_feat)
}
# Here comes the main function for resampling of mswsd values
# start with all features than in steps of 5% remove features until 10% of the original features remain
# start by removing the most variable ones first, in total we will have 19 different levels
# to test the stability of the mswsd values we employ a resampling approach with 100 iterations
resamp_mswsd<-function(metabo.data,Results_folder_Preprocessing_folder, OutputName)
{
  s <- seq(5, 100, 5) #prepare steps for feature reduction
  vari <- matrix(ncol = 20, nrow = 100) #to store obtained mswsd values
  no_samp<-ncol(metabo.data) #number of samples
  part_samp<-floor(no_samp*0.66) # use two thirds of samples
  
  for(i in 1:length(s)) #variables=5 i=20
  {
    variables <- s[i]
    for(k in 1:100)	# 100 times resampling # k=100 #sel=72
    {
      sel <- sample(1:no_samp, part_samp) # randomly select two thirds out of all samples
      vari[k,i] <- mswsd(metabo.data[,sel], variables)
    }
    # normality check on total spectral areas of remaining features
    rem<-floor(nrow(metabo.data)-(nrow(metabo.data)*variables/100))
    # calculate variance
    vari2<-apply(metabo.data,1,var) #calculate variance for each feature across samples
    # remove features with highest variance
    sel2<-order(vari2)[1:(nrow(metabo.data)-rem)]
    metabo.data2<-metabo.data[sel2,] #data matrix that contains only features with low variance
    total <- apply(metabo.data2, 2, sum)
    pp<-shapiro.test(total)$p.value
    cat("pecentage of features used=",variables,"p-value of total areas=",pp,"\n")
    rm(vari2,rem,sel2,metabo.data2)
  }
  # prepare plot
  
  pdf(paste(Results_folder_Preprocessing_folder,"/resamp_mswsd", OutputName,".pdf",sep=""), onefile = FALSE)
  box.mswsd <- boxplot(vari, use.cols = TRUE, names = as.character(s), xlab = "percentage of features", ylab = "mswsd", main = "Resampling of mswsd",  ylim = c(0, 100))
  dev.off()
  boxplot(vari, use.cols = TRUE, names = as.character(s), xlab = "percentage of features", ylab = "mswsd", main = "Resampling of mswsd",  ylim = c(0, 100))
}
# From the resampling of mswsd values plot you may have selected the ideal percentage of features used for data normalization
# With this percentage we apply now the actual data normalization
norm_unbal<-function(metabo.data, use_percent,flag,Results_folder_Preprocessing_folder)
{       
  metabo.data<-abs(metabo.data)
  # Identify normalization reference features
  cat("#### percentage of used features ######","\n")
  cat("percentage of used features", use_percent,"\n")
  # determination of invariant features
  rem<-floor(nrow(metabo.data)-(nrow(metabo.data)*use_percent/100))
  # calculate variance
  varian<-apply(metabo.data,1,var) #calculate variance for each feature across samples
  # remove features with highest variance
  sel<-order(varian)[1:(nrow(metabo.data)-rem)]
  metabo.data2<-metabo.data[sel,] #data matrix that contains only features with low variance
  
  # perform actual normalization
  # linear baseline normalization based on mean values
  if(flag=="LBME")
  {
    linear.baseline <- apply(metabo.data2,1,median) #compute baseline
    baseline.mean <- mean(linear.baseline)
    sample.means <- apply(metabo.data2,2,mean)
    linear.scaling <- baseline.mean/sample.means
    cat("mean lin sca=",mean(linear.scaling),"sd lin sca=",sd(linear.scaling),"rsd lin sca=",sd(linear.scaling)/mean(linear.scaling),"\n")
    cat("min lin sca=",min(linear.scaling),"max lin sca=",max(linear.scaling),"\n")
    norm.metabo.data <- t(t(metabo.data)*linear.scaling) 
  }
  # linear baseline normalization based on median values
  if (flag=="LBMD")
  {
    linear.baseline <- apply(metabo.data2,1,median) #compute baseline
    baseline.median<-median(linear.baseline)
    sample.medians<-apply(metabo.data2,2,median)
    linear.scaling<-baseline.median/sample.medians
    cat("mean lin sca=",mean(linear.scaling),"sd lin sca=",sd(linear.scaling),"rsd lin sca=",sd(linear.scaling)/mean(linear.scaling),"\n")
    cat("min lin sca=",min(linear.scaling),"max lin sca=",max(linear.scaling),"\n")
    norm.metabo.data <- t(t(metabo.data)*linear.scaling)
  }
  #VSN
  if(flag=="VSN")
  {
    library(vsn)
    vsn.model<-vsn2(metabo.data2)
    norm.metabo.data<-predict(vsn.model,(metabo.data))
  }
  #PQN
  if(flag=="PQN")
  {
    reference<-apply(metabo.data2,1,median)
    quotient<-metabo.data2/reference
    quotient.median<-apply(quotient,2,median)
    norm.metabo.data<-t(t(metabo.data)/quotient.median)
  }
  return(norm.metabo.data)
  rm(varian,sel,rem)
}
unbal_reg<-function(metabo.data)
{
  metabo.data<-abs(metabo.data)
  {
    cat("if p-value smaller than 0.05 total spectral areas are not normally distributed-> unbalanced regulation","\n")
    total <- apply(metabo.data, 2, sum)
    shapiro.test(total)
  }
}
