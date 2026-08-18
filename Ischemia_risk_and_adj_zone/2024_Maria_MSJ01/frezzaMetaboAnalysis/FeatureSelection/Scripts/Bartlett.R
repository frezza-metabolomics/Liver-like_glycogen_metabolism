###########################################################################################
### ### ### Bartlett function: Internal Function to perform Bartlett test and plots ### ###
###########################################################################################

#' @param DMA_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param DMA_SettingsFile DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param DMA_SettingsInfo \emph{Optional: } Named vector including the information about the Conditions column c(Conditions="ColumnName_Plot_SettingsFile"). Can additionally pass information on numerator or denominator c(numerator = "ColumnName_Plot_SettingsFile", denumerator = "ColumnName_Plot_SettingsFile") for specifying which comparison(s) will be done (one-vs-one, all-vs-one, all-vs-all). Using =NULL selects all the condition and performs multiple comparison all-vs-all. Log2FC are obtained by dividing the numerator by the denominator, thus positive Log2FC values mean higher expression in the numerator and are presented in the right side on the Volcano plot (For CoRe the Log2Distance). \strong{Default = c(Conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param OutputName String which is added to the output files of the DMA.
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{Default = "csv"}
#' @param Plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param FolderName {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @keywords Bartlett test,Normality testing, Density plot, QQplot
#' @export
#'

BartlettTest <-function(Result_Folder = Folder_paths$DMAFolder,
                        RDS_folder = Folder_paths$DMARDS,
                        scriptpath= scriptpath,
                        DMA_data,
                        DMA_SettingsFile,
                        DMA_SettingsInfo = c(Conditions="Conditions", numerator = NULL, denumerator = NULL),
                        ForceRun = FALSE,
                        Plot=TRUE,
                        
                        Save_as_Results="csv",
                        Save_as_Plot="svg",
                        
                        FolderName = NULL,
                        OutputName = NULL
){
  


  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)){ OutputName <- "BarlettTestRes.Rdata"
  }else{ OutputName <- paste("BarlettTest_" ,OutputName, ".Rdata" ,sep = "" )}
  if (file.exists(file.path(RDS_folder,  OutputName))==TRUE & ForceRun == FALSE){
    message("Loading already existing ", paste(OutputName))
    load(file.path(RDS_folder,  OutputName), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running Load data script main body.")
    
    
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
    if(is.logical(Plot) == FALSE){
      stop("Check input. The plot value should be either =TRUE if a Volcano plot presenting the DMA results is to be exported or =FALSE if not.")
    }
    Save_as_Plot_options <- c("svg","pdf","png")
    if(Save_as_Plot %in% Save_as_Plot_options == FALSE){
      stop("Check input. The selected Save_as_Plot option is not valid. Please select one of the folowwing: ",paste(Save_as_Plot_options,collapse = ", "),"." )
    }
    Save_as_Results_options <- c("txt","csv", "xlsx" )
    if(Save_as_Results %in% Save_as_Results_options == FALSE){
      stop("Check input. The selected Save_as_Results option is not valid. Please select one of the folowwing: ",paste(Save_as_Results_options,collapse = ", "),"." )
    }
    
    
    
    DMA_data <- DMA_data %>% filter(DMA_SettingsFile$Conditions %in% numerator |   DMA_SettingsFile$Conditions %in% denominator )
    DMA_SettingsFile <- DMA_SettingsFile %>% filter(DMA_SettingsFile$Conditions %in% numerator |   DMA_SettingsFile$Conditions %in% denominator )  
    
    Conditions <- DMA_SettingsFile$Conditions
    
    # Use Bartletts test
    bartlett_res =  apply(DMA_data,2,function(x) bartlett.test(x~Conditions))
    
    #Make the output DF
    DF_bartlett_results <- as.data.frame(matrix(NA, nrow = ncol(DMA_data)), ncol = 1)
    rownames(DF_bartlett_results) <- colnames(DMA_data)
    colnames(DF_bartlett_results) <- "Bartlett p.val"
    
    
    for(l in 1:length(bartlett_res)){
      DF_bartlett_results[l, 1] <-bartlett_res[[l]]$p.value
    }
    DF_bartlett_results <- DF_bartlett_results %>% mutate(`Var homogeneity`= case_when(`Bartlett p.val`< 0.05~ FALSE,
                                                                                       `Bartlett p.val`>=0.05 ~ TRUE))
    # if p<0.05 then unequal variances
    paste("For",round(sum(DF_bartlett_results$`Var homogeneity`, na.rm = T)/  nrow(DF_bartlett_results), digits = 4) * 100, "% of metabolites the group variances are equal.")
    
    DF_bartlett_results <- DF_bartlett_results %>% rownames_to_column("Metabolite") %>% relocate("Metabolite")
    DF_Bartlett_results_out <- DF_bartlett_results
    
    # Save the DF Bartlett
    if (Save_as_Results == "xlsx"){
      writexl::write_xlsx(DF_Bartlett_results_out,paste(Result_Folder,"/DF_Bartlett_results_table",".",Save_as_Results,sep =  "")) # save the DMA result DF
    }else if (Save_as_Results == "csv"){
      write.csv(DF_Bartlett_results_out,paste(Result_Folder,"/DF_Bartlett_results_table",".",Save_as_Results,sep =  ""),row.names =FALSE) # save the DMA result DF
    }else if (Save_as_Results == "txt"){
      write.table(DF_Bartlett_results_out,paste(Result_Folder,"/DF_Bartlett_results_table",".",Save_as_Results,sep =  ""), col.names = TRUE, row.names = FALSE) # save the DMA result DF
    }
    
    # Make density plots
    invisible( Bartlettplot <-ggplot(data.frame(x = DF_Bartlett_results_out), aes(x =DF_Bartlett_results_out$`Bartlett p.val`)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 20)  +
      geom_density(alpha = 0.2, fill = "grey45")+ ggtitle("Bartlett's test p.value distribution") +
      xlab("p.value")+ geom_vline(aes(xintercept = 0.05, color="darkred"))+ labs(subtitle = paste("For ",round(sum(DF_bartlett_results$`Var homogeneity`, na.rm = T)/  nrow(DF_bartlett_results), digits = 4) * 100, "% of metabolites the group variances are equal."))
    )
    
    # Do we save the pvalue density plot?
    ggsave(filename = paste0(Result_Folder, "/Bartlett_Density_plot.",Save_as_Plot), plot = Bartlettplot, width = 10,  height = 8)
    
    if(Plot == TRUE){
      plot(Bartlettplot)
    }
    
    message("For ",round(sum(DF_bartlett_results$`Var homogeneity`, na.rm = T)/  nrow(DF_bartlett_results), digits = 4) * 100, "% of metabolites the group variances are equal.")
    
    
    # Save results
    # writexl::write_xlsx(filtered_matrix, paste(Results_folder_Preprocessing_folder, "/Feature_Filtered_Matrix.xlsx", sep = ""))#,showNA = TRUE)
    # Make Result list
    BartlettTestRes =  list("DF"=DF_Bartlett_results_out , "Plot"= Bartlettplot)
    # Save output RDS
    save(BartlettTestRes, file = file.path(RDS_folder,  OutputName))
  }
  return(suppressWarnings(invisible(BartlettTestRes)))
}


