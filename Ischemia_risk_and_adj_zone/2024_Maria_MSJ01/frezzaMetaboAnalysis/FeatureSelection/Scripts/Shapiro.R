
#############################################################################################
### ### ### Shapiro function: Internal Function to perform Shapiro test and plots ### ### ###
#############################################################################################

#' @param DMA_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param DMA_SettingsFile DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param DMA_SettingsInfo \emph{Optional: } Named vector including the information about the Conditions column c(Conditions="ColumnName_Plot_SettingsFile"). Can additionally pass information on numerator or denominator c(numerator = "ColumnName_Plot_SettingsFile", denumerator = "ColumnName_Plot_SettingsFile") for specifying which comparison(s) will be done (one-vs-one, all-vs-one, all-vs-all). Using =NULL selects all the condition and performs multiple comparison all-vs-all. Log2FC are obtained by dividing the numerator by the denominator, thus positive Log2FC values mean higher expression in the numerator and are presented in the right side on the Volcano plot (For CoRe the Log2Distance). \strong{Default = c(Conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param STAT_pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test" or "cor.test", for one-vs-all or all-vs-all comparison choose aov (=annova), kruskal.test or lmFit (=limma) \strong{Default = "t-test"}
#' @param OutputName String which is added to the output files of the DMA.
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used \strong{Default = FALSE}
#' @param QQplots \emph {Optional: } TRUE or FALSE for whether QQ plots should be plotted  \strong{Default = TRUE}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{Default = "csv"}
#' @param Plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param FolderName {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @keywords Shapiro test,Normality testing, Density plot, QQplot
#' @export
#'

ShapiroTest <- function(Result_Folder,
                        RDS_folder,
                        scriptpath,
                        
                        DMA_data,
                   DMA_SettingsFile,
                   DMA_SettingsInfo = c(Conditions="Conditions", numerator = NULL, denumerator = NULL),
                   ForceRun = FALSE,
                   CoRe=FALSE,
                   QQplots=FALSE,
                   
                   Save_as_Results="csv",
                   Save_as_Plot="svg",
                   
                   FolderName = NULL,
                   OutputName = NULL
){
  

  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)){ OutputName <- "ShapiroTestRes.Rdata"
  }else{ OutputName <- paste("ShapiroTestRes_" ,OutputName, ".Rdata" ,sep = "" )}
  if (file.exists(file.path(RDS_folder,  OutputName))==TRUE  & ForceRun == FALSE){
    message("Loading already existing ", paste(OutputName))
    load(file.path(RDS_folder,  OutputName), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running Load data script main body.")
    
    
    ## Input checks
    ################################################################################################################################################################################################
    ## ------------ Check Input files ----------- ##
    #1. DMA_data and Conditions
    if(class(DMA_data) != "data.frame"){
      stop("DMA_data should be a data.frame. It's currently a ", paste(class(DMA_data), ".",sep = ""))
    }
    if(any(duplicated(row.names(DMA_data)))==TRUE){
      stop("Duplicated row.names of DMA_data, whilst row.names must be unique")
    } else{
      Test_num <- apply(DMA_data, 2, function(x) is.numeric(x))
      if((any(Test_num) ==  FALSE) ==  TRUE){
        stop("DMA_data needs to be of class numeric")
      } else{
        Test_match <- merge(DMA_SettingsFile, DMA_data, by.x = "row.names",by.y = "row.names", all =  FALSE, sort = FALSE) # Do the unique IDs of the "DMA_data" match the row names of the "DMA_SettingsFile"?
        if(nrow(Test_match) ==  0){
          stop("row.names DMA_data need to match row.names DMA_SettingsFile.")
        } else{
          DMA_data <- DMA_data
        }
      }
    }
    
    ## ------------ Check Input SettingsInfo ----------- ##
    #3. DMA_SettingsInfo
    if(DMA_SettingsInfo[["Conditions"]] %in% DMA_SettingsInfo==TRUE){
      if(DMA_SettingsInfo[["Conditions"]] %in% colnames(DMA_SettingsFile)== FALSE){
        stop("The ",DMA_SettingsInfo[["Conditions"]], " column selected as Conditions in DMA_SettingsInfo was not found in DMA_SettingsFile. Please check your input.")
      }else{# if true rename to Conditions
        DMA_SettingsFile<- DMA_SettingsFile%>%
          dplyr::rename("Conditions"= paste(DMA_SettingsInfo[["Conditions"]]) )
      }
    }else{
      stop("You have to provide a DMA_SettingsInfo for Conditions.")
    }
    
    ##########################
    if("denominator" %in% names(DMA_SettingsInfo)==TRUE){
      if(DMA_SettingsInfo[["denominator"]] %in% DMA_SettingsFile$Conditions==FALSE){
        stop("The ",DMA_SettingsInfo[["denominator"]], " column selected as denominator in DMA_SettingsInfo was not found in DMA_SettingsFile. Please check your input.")
      }else{
        denominator <- DMA_SettingsInfo[["denominator"]]
      }
    }
    if("numerator" %in% names(DMA_SettingsInfo)==TRUE){
      if(DMA_SettingsInfo[["numerator"]] %in% DMA_SettingsFile$Conditions  == FALSE){
        stop("The ",DMA_SettingsInfo[["numerator"]], " column selected as numerator in DMA_SettingsInfo was not found in DMA_SettingsFile. Please check your input.")
      }else{
        numerator <- DMA_SettingsInfo[["numerator"]]
      }
    }
    if("denominator" %in% names(DMA_SettingsInfo)==FALSE  & "numerator" %in% names(DMA_SettingsInfo) ==TRUE){
      stop("Check input. The selected denominator option is empty while ",paste(DMA_SettingsInfo[["numerator"]])," has been selected as a numerator. Please add a denminator for 1-vs-1 comparison or remove the numerator for all-vs-all comparison." )
    }
    
    ## ------------ Check Denominator/numerator ----------- ##
    #4.  Denominator and numerator: Define if we compare one_vs_one, one_vs_all or all_vs_all.
    if("denominator" %in% names(DMA_SettingsInfo)==FALSE  & "numerator" %in% names(DMA_SettingsInfo) ==FALSE){
      Conditions = DMA_SettingsFile$Conditions
      denominator <-unique(DMA_SettingsFile$Conditions)
      numerator <-unique(DMA_SettingsFile$Conditions)
    }else if("denominator" %in% names(DMA_SettingsInfo)==TRUE  & "numerator" %in% names(DMA_SettingsInfo)==FALSE){
      #all-vs-one: Generate the pairwise combinations
      Conditions = DMA_SettingsFile$Conditions
      denominator <- DMA_SettingsInfo[["denominator"]]
      numerator <-unique(DMA_SettingsFile$Conditions)
    }
    
    ## ------------ Check General parameters ----------- ##
    #6. General parameters
    Save_as_Plot_options <- c("svg","pdf","png")
    if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
      stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
    }
    Save_as_Results_options <- c("txt","csv", "xlsx" )
    if(Save_as_Results %in% Save_as_Results_options == FALSE){
      stop("Check input. The selected Save_as_Results option is not valid. Please select one of the folowwing: ",paste(Save_as_Results_options,collapse = ", "),"." )
    }
    
    #7. Are sample numbers enough?
    Num <- DMA_data %>%
      filter(DMA_SettingsFile$Conditions %in% numerator) %>%
      select_if(is.numeric)#only keep numeric columns with metabolite values
    Denom <- DMA_data %>%
      filter(DMA_SettingsFile$Conditions %in% denominator) %>%
      select_if(is.numeric)
    
    if(nrow(Num)==1){
      stop("There is only one sample available for ", numerator, ", so no statistical test can be performed.")
    } else if(nrow(Denom)==1){
      stop("There is only one sample available for ", denominator, ", so no statistical test can be performed.")
    }else if(nrow(Num)==0){
      stop("There is no sample available for ", numerator, ".")
    }else if(nrow(Denom)==0){
      stop("There is no sample available for ", denominator, ".")
    }
    
    ###############################################################################################################################################################################################################
    ## ------------ Check data normality and statistical test chosen and generate Output DF----------- ##
    # Before Hypothesis testing, we have to decide whether to use a parametric or a non parametric test. We can test the data normality using the Shapiro test.
    ##-------- First: Load the data and perform the shapiro.test on each metabolite across the samples of one condition. this needs to be repeated for each condition:
    #Prepare the input:
    DMA_shaptest <- replace(DMA_data, DMA_data==0, NA) %>% #Shapiro test ignores NAs!
      filter(DMA_SettingsFile$Conditions %in% numerator | DMA_SettingsFile$Conditions %in% denominator)%>%
      select_if(is.numeric)
    temp<- as.vector(sapply(DMA_shaptest, function(x) var(x)) == 0)#  we have to remove features with zero variance if there are any.
    DMA_shaptest <- DMA_shaptest[,!temp]
    DMA_shaptest_Cond <-merge(data.frame(Conditions = DMA_SettingsFile[, "Conditions", drop = FALSE]), DMA_shaptest, by=0, all.y=TRUE, sort = FALSE)
    
    UniqueConditions <- DMA_SettingsFile%>%
      subset(DMA_SettingsFile$Conditions %in% numerator | DMA_SettingsFile$Conditions %in% denominator, select = c(Conditions))
    UniqueConditions <- unique(UniqueConditions$Conditions)
    
    #Generate the results
    shapiro_results <- list()
    for (i in UniqueConditions) { # i = UniqueConditions[2]
      # Subset the data for the current condition
      subset_data <- DMA_shaptest_Cond%>%
        column_to_rownames("Row.names")%>%
        subset(Conditions == i, select = -c(1))
      
      sd_values <- sapply(subset_data, sd)
      sd_values <- as.vector(sd_values==0) 

      subset_data <- subset_data[,!sd_values]
      
      # Apply Shapiro-Wilk test to each feature in the subset
      shapiro_results[[i]] <- as.data.frame(sapply(subset_data, function(x) shapiro.test(x)))
    }
    
    #Make the output DF
    DF_shapiro_results <- as.data.frame(matrix(NA, nrow = length(UniqueConditions), ncol = ncol(DMA_shaptest)))
    rownames(DF_shapiro_results) <- UniqueConditions
    colnames(DF_shapiro_results) <- colnames(DMA_shaptest)
    
    for(k in 1:length(UniqueConditions)){# for each group  k=1 HH, k=2 WT
      for(l in 1:ncol(DMA_shaptest)){ # for each feature l=2
        data <- shapiro_results[[k]]
        if(colnames(DF_shapiro_results)[l] %in% colnames(data)){
          DF_shapiro_results[k, l]<- data[2,which(colnames(data) %in% colnames(DF_shapiro_results)[l])]
        }else{
          DF_shapiro_results[k, l]<- 1
        }

      }
    }
    colnames(DF_shapiro_results) <- paste("Shapiro p.val(", colnames(DF_shapiro_results),")", sep = "")
    
    ##------ Second: Give feedback to the user if the chosen test fits the data distribution. The data are normal if the p-value of the shapiro.test > 0.05.
    Density_plots <- list()
    if(QQplots==TRUE){
      QQ_plots <- list()
    }
    for(x in 1:nrow(DF_shapiro_results)){
      transpose <- as.data.frame(t(DF_shapiro_results[x,]))
      Norm <- format((round(sum(transpose[[1]] > 0.05)/nrow(transpose),4))*100, nsmall = 2) # Percentage of normally distributed metabolites across samples
      NotNorm <- format((round(sum(transpose[[1]] < 0.05)/nrow(transpose),4))*100, nsmall = 2) # Percentage of not-normally distributed metabolites across samples
      
      message("For the condition ", colnames(transpose) ," ", Norm, " % of the metabolites follow a normal distribution and ", NotNorm, " % of the metabolites are not-normally distributed according to the shapiro test. `shapiro.test` ignores missing values in the calculation.")
      
      
      # Assign the calculated values to the corresponding rows in result_df
      DF_shapiro_results$`Metabolites with normal distribution [%]`[x] <- Norm
      DF_shapiro_results$`Metabolites with not-normal distribution [%]`[x] <- NotNorm
      
      #reorder the DF:
      DF_shapiro_results<-DF_shapiro_results[,c(ncol(DF_shapiro_results)-1, ncol(DF_shapiro_results), 1:(ncol(DF_shapiro_results)-2))]
      
      DF_shapiro_results_out<- t(DF_shapiro_results)%>% as.data.frame()%>% rownames_to_column("Shapiro_p.val")
      DF_shapiro_results_out$Shapiro_p.val <-  str_replace_all(DF_shapiro_results_out$Shapiro_p.val, "Shapiro p.val", " ")
      DF_shapiro_results_out$Shapiro_p.val <-gsub("[[:punct:]]", " ", DF_shapiro_results_out$Shapiro_p.val)
      
      # Save the DF Shapiro
      if (Save_as_Results == "xlsx"){
        writexl::write_xlsx(DF_shapiro_results_out,paste(Result_Folder,"/DF_shapiro_results_",".",Save_as_Results,sep =  "")) # save the DMA result DF
      }else if (Save_as_Results == "csv"){
        write.csv(DF_shapiro_results_out,paste(Result_Folder,"/DF_shapiro_results",".",Save_as_Results,sep =  ""),row.names =FALSE) # save the DMA result DF
      }else if (Save_as_Results == "txt"){
        write.table(DF_shapiro_results_out,paste(Result_Folder,"/DF_shapiro_results",".",Save_as_Results,sep =  ""), col.names = TRUE, row.names = FALSE) # save the DMA result DF
      }
      
      ## Make Group wise data distribution plot and QQ plots
      subset_data <- DMA_shaptest_Cond%>%
        column_to_rownames("Row.names")%>%
        subset(Conditions ==  colnames(transpose), select = -c(1))
      all_data <- unlist(subset_data)
      
      plot <- ggplot(data.frame(x = all_data), aes(x = x)) +
        geom_histogram(aes(y=..density..), binwidth=.5, colour="black", fill="white")  +
        geom_density(alpha = 0.2, fill = "grey45")
      
      density_values <- ggplot_build(plot)$data[[2]]
      
      plot <- ggplot(data.frame(x = all_data), aes(x = x)) +
        geom_histogram(aes(y=..density..), binwidth=.5, colour="black", fill="white") +
        geom_density(alpha=.2, fill="grey45") +
        scale_x_continuous(limits = c(0, density_values$x[max(which(density_values$scaled >= 0.1))]))
      
      density_values2 <- ggplot_build(plot)$data[[2]]
      
      suppressWarnings( sampleDist <- ggplot(data.frame(x = all_data), aes(x = x)) +
                          geom_histogram(aes(x = all_data, y = ..density..), bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
                          geom_density(alpha=.2, fill="red") +
                          scale_x_continuous(limits = c(0, density_values$x[max(which(density_values$scaled >= 0.1))])) +
                          theme_minimal()+
                          # geom_vline(xintercept =median(all_data) , linetype = "dashed", color = "red")+
                          labs(title=paste("Data distribution ",  colnames(transpose)), subtitle = paste(NotNorm, "% of metabolites not normally distributed based on Shapiro test"),x="Normalized Intensity values", y = "Density")#+
                        # geom_text(aes(x = density_values2$x[which.max(density_values2$y)], y = 0, label = "Median"),  vjust = 0, hjust = -0.5, color = "red", size = 3.5)  # Add label for
      )
      
      plot(sampleDist)
      Density_plots[[paste(colnames(transpose))]] <- recordPlot()
      
      
      if(CoRe==TRUE){
        ggsave(filename = paste0(Result_Folder, "/Density_plot_", paste(colnames(transpose)),".",Save_as_Plot), plot = sampleDist, width = 10,  height = 8)
      }else{
        ggsave(filename = paste0(Result_Folder, "/Density_plot_", paste(colnames(transpose)),".",Save_as_Plot), plot = sampleDist, width = 10,  height = 8)
        
      }
      # QQ plots
      if(QQplots==TRUE){
        # Make folders !has to be moved on top!
        conds <- unique(c(numerator, denominator))
        for(x in conds){
          Result_Folder_Condition <- file.path(Result_Folder, paste(x)) # Make DMA results folder
          if (!dir.exists(Result_Folder_Condition)) {dir.create(Result_Folder_Condition)}
        }
        #QQ plots for each groups for each metabolite for normality visual check
        qq_plot_list <- list()
        for (col_name in colnames(subset_data)){
          qq_plot <- ggplot(data.frame(x = subset_data[[col_name]]), aes(sample = x)) +
            geom_qq() +
            geom_qq_line(color = "red") +
            labs(title = paste("QQPlot for", col_name),x = "Theoretical", y="Sample")+ theme_minimal()
          
          plot.new()
          plot(qq_plot)
          qq_plot_list[[col_name]] <-  recordPlot()
          
          col_name2 <- (gsub("/","_",col_name))#remove "/" cause this can not be safed in a PDF name
          col_name2 <- gsub("-", "", col_name2)
          col_name2 <- gsub("/", "", col_name2)
          col_name2 <- gsub(" ", "", col_name2)
          col_name2 <- gsub("\\*", "", col_name2)
          col_name2 <- gsub("\\+", "", col_name2)
          col_name2 <- gsub(",", "", col_name2)
          col_name2 <- gsub("\\(", "", col_name2)
          col_name2 <- gsub("\\)", "", col_name2)
          
          ggsave(paste0(Result_Folder, "/", paste(colnames(transpose)),"/",paste(col_name2),".",Save_as_Plot), plot = qq_plot, device = Save_as_Plot, width = 10,  height = 8)
          
          dev.off()
        }
        QQ_plots[[paste(colnames(transpose))]] <- qq_plot_list
      }
    }
    
    # Save results
    #  We save files inside the script writexl::write_xlsx(filtered_matrix, paste(Results_folder_Preprocessing_folder, "/Feature_Filtered_Matrix.xlsx", sep = ""))#,showNA = TRUE)
    # Make Result list
    if(QQplots==TRUE){
      ShapiroTestRes <- list("DF" = list("Shapiro_result"=DF_shapiro_results),"Plot"=list("Distributions"=Density_plots, "QQ_plots" = QQ_plots))
    }else{
      ShapiroTestRes <- list("DF" = list("Shapiro_result"=DF_shapiro_results),"Plot"=list("Distributions"=Density_plots))
    }
    # Save output RDS
    save(ShapiroTestRes, file = file.path(RDS_folder,  OutputName))
  }
  return(suppressWarnings(invisible(ShapiroTestRes)))
}

