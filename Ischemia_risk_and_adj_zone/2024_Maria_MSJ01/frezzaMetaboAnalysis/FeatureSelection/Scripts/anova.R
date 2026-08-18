
################################################################
### ### ### AOV: Internal Function to perform anova  ### ### ###
################################################################

#' @param DMA_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param Conditions Factor with sample group information.
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param Log2FC_table Table with Metabolites are rows and a Log2FC column \strong(Default = Log2FC_table)
#'
#' @keywords Kruskal test,Hypothesis testing, p.value
#' @export

AnovaTest <-function(DMA_data,
                     DMA_SettingsFile,
                     DMA_SettingsInfo = c(Conditions="Conditions", numerator = NULL, denumerator = NULL),
                     Log2FC_table,
                     PCutoff= 0.05,
                     FCcutoff= 1,
                     CoRe=FALSE,                     
                     STAT_padj = "fdr",
                     Labels = TRUE,
                     ForceRun = FALSE,
                     
                     Save_as_Results = "xlsx",
                     Save_as_Plot = "svg",
                     
                     ProjectFolder,
                     FolderName = NULL,
                     OutputName = NULL
){
  
  ## 1. ------------ Setup and installs ----------- ##
  
  RequiredPackages <- c("tidyverse" # general scripting
  )
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressWarnings(suppressMessages(library(tidyverse)))
  # usethis::edit_file(file.path(scriptpath,"Volcano.R"))
  source((file.path(scriptpath,"Volcano.R")))
  
  
  ## 1,5. -------- Input Checks ------------- ##
  ## ------------ Check Input SettingsInfo ----------- ##
  #3. DMA_SettingsInfo
  if(DMA_SettingsInfo[["Conditions"]] %in% DMA_SettingsInfo==TRUE){
    if(DMA_SettingsInfo[["Conditions"]] %in% colnames(DMA_SettingsFile)== FALSE){
      stop("The ",DMA_SettingsInfo[["Conditions"]], " column selected as Conditions in DMA_SettingsInfo was not found in DMA_SettingsFile. Please check your input.")
    }else{# if true rename to Conditions
      DMA_SettingsFile<- DMA_SettingsFile %>% dplyr::rename("Conditions"= paste(DMA_SettingsInfo[["Conditions"]]) )
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
    # all-vs-all: Generate all pairwise combinations
    Conditions = DMA_SettingsFile$Conditions
    denominator <-unique(DMA_SettingsFile$Conditions)
    numerator <-unique(DMA_SettingsFile$Conditions)
    comparisons <- combn(unique(Conditions), 2) %>% as.matrix()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = TRUE
  }else if("denominator" %in% names(DMA_SettingsInfo)==TRUE  & "numerator" %in% names(DMA_SettingsInfo)==FALSE){
    #all-vs-one: Generate the pairwise combinations
    Conditions = DMA_SettingsFile$Conditions
    denominator <- DMA_SettingsInfo[["denominator"]]
    numerator <-unique(DMA_SettingsFile$Conditions)
    # Remove denom from num
    numerator <- numerator[!numerator %in% denominator]
    comparisons  <- t(expand.grid(numerator, denominator)) %>% as.data.frame()
    #Settings:
    MultipleComparison = TRUE
    all_vs_all = FALSE
  }else if("denominator" %in% names(DMA_SettingsInfo)==TRUE  & "numerator" %in% names(DMA_SettingsInfo)==TRUE){
    # one-vs-one: Generate the comparisons
    comparisons <- matrix(c(numerator, denominator))
    Conditions = DMA_SettingsFile$Conditions
    #Settings:
    MultipleComparison = FALSE
    all_vs_all = FALSE
  }
  
  
  ## 2. ------------ Create Results output folder ------------ ##
  if(is.null(FolderName)){name <- paste("Results",sep = "" ) }else{ name <- paste("Results",FolderName,sep = "_" )}
  ProjectCode <- strsplit(ProjectFolder, "_")[[1]][[ length(strsplit(ProjectFolder, "_")[[1]])]] 
  Analysed_data_folder <- file.path(ProjectFolder,paste("5.", ProjectCode, "Analysed Data", sep = " ")) 
  RDS_folder <- file.path(ProjectFolder,paste("5.", ProjectCode, "Analysed","Data"), str_replace(name, "Results", "RDS")) 
  if (!dir.exists(RDS_folder)) {dir.create(RDS_folder)} 
  
  Results_folder <- file.path(Analysed_data_folder, name) # assign results folder
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_DMA_folder <- file.path(Results_folder,"DMA") # Make DMA results folder
  if (!dir.exists(Results_folder_DMA_folder)) {dir.create(Results_folder_DMA_folder)}
  
  if(MultipleComparison==TRUE){
    if (length(denominator)==1){
      Results_folder_Conditions <- file.path(Results_folder_DMA_folder,  paste0("All","_vs_",toString(denominator)))
    }else{
      Results_folder_Conditions <- file.path(Results_folder_DMA_folder,  paste0("All","_vs_","All"))
    }
  }else{
    Results_folder_Conditions <- file.path(Results_folder_DMA_folder,  paste0(toString(numerator),"_vs_",toString(denominator)))
  }
  if (!dir.exists(Results_folder_Conditions)) {dir.create(Results_folder_Conditions)}
  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)){ OutputName <- "AnovaTestRes.Rdata"
  }else{ OutputName <- paste("AnovaTest_" ,OutputName, ".Rdata" ,sep = "" )}
   if (file.exists(file.path(RDS_folder,  OutputName))==TRUE & ForceRun == FALSE){
     message("Loading already existing ", paste(OutputName))
     load(file.path(RDS_folder,  OutputName), .GlobalEnv)
    
  #   ## 4. ------------ Do calculations ------------- ##
   }else{
    message("Running Load data script main body.")
    
    ## 1. Anova p.val
    aov.res= apply(DMA_data,2,function(x) aov(x~Conditions))
    
    # Extract p-values
    p_values <- sapply(aov.res, function(x) summary(x)[[1]][["Pr(>F)"]][1])
    message(paste("using the annova test there are",sum(p_values < 0.05),"out of",length(names(aov.res)), "with < 0.05"))
    ## 2. Tukey test p.adj
    posthoc.res = lapply(aov.res, TukeyHSD, conf.level=0.95)
    Tukey_res <- do.call('rbind', lapply(posthoc.res, function(x) x[1][[1]][,'p adj']))
    Tukey_res <- as.data.frame(Tukey_res)
    
    comps <-   paste(comparisons[1, ], comparisons[2, ], sep="-")# normal
    opp_comps <-  paste(comparisons[2, ], comparisons[1, ], sep="-")
    
    if(sum(opp_comps %in%  colnames(Tukey_res))>0){# if opposite comparisons is true
      for (comp in 1: length(opp_comps)){
        colnames(Tukey_res)[colnames(Tukey_res) %in% opp_comps[comp]] <-  comps[comp]
      }
    }
    
    Tukey_res[Tukey_res==0]<- 0.000000000000001
    
    ## 3. t.val
    Tukey_res_diff <- do.call('rbind', lapply(posthoc.res, function(x) x[1][[1]][,'diff'])) %>% as.data.frame()
    
    if (sum(opp_comps %in%  colnames(Tukey_res_diff))>0){# if oposite comparisons is true
      for (comp in 1: length(opp_comps)){
        colnames(Tukey_res_diff)[colnames(Tukey_res_diff) %in% opp_comps[comp]] <-  comps[comp]
      }
    }
    
    #Make output DFs:
    Pval_table <- Tukey_res
    Pval_table <- rownames_to_column(Pval_table,"Metabolite")
    
    Tval_table <- rownames_to_column(Tukey_res_diff,"Metabolite")
    
    common_col_names <- setdiff(names(Tukey_res_diff), "row.names")#Here we need to adapt for one_vs_all or all_vs_all
    
    results_list <- list()
    for(col_name in common_col_names){
      # Create a new data frame by merging the two data frames
      merged_df <- merge(Pval_table[,c("Metabolite",col_name)], Tval_table[,c("Metabolite",col_name)], by="Metabolite", all=TRUE, sort = FALSE)%>%
        dplyr::rename("p.adj"=2,
                      "t.val"=3)
      
      #We need to add _vs_ into the comparison col_name
      pattern <- paste(Conditions, collapse = "|")
      conditions_present <- unique(unlist(regmatches(col_name, gregexpr(pattern, col_name))))
      # modified_col_name <- paste(conditions_present[1], "vs", conditions_present[2], sep = "_")
      modified_col_name <- paste(conditions_present[1], conditions_present[2], sep = "-")
      
      # Add the new data frame to the list with the column name as the list element name
      results_list[[modified_col_name]] <- merged_df
    }
    
    # Merge the data frames in list1 and list2 based on the "Metabolite" column
    list_names <-  names(results_list)
    
    #
    colnames(Log2FC_table) <-  str_replace( colnames(Log2FC_table), "Log2FC_", "")
    
    merged_list <- list()
    for(name in list_names){ # name = list_names[121]
      # Check if the data frames exist in both lists
      if(name %in% names(results_list) && name %in% names(Log2FC_table)){
        merged_df <- merge(results_list[[name]], Log2FC_table %>% select(Metabolite,name), by = "Metabolite", all = TRUE, sort = FALSE)
        merged_df <- merged_df[,c(1,4,2,3)]#reorder the columns
        colnames(merged_df)[2] <- "Log2FC"
        merged_list[[name]] <- merged_df
      }
    }
    
    # # Make sure the right comparisons are returned:
    # if(all_vs_all==TRUE){
    #   STAT_C1vC2 <- merged_list
    # }else if(all_vs_all==FALSE){
    #   #remove the comparisons that are not needed:
    #   modified_df_list <- list()
    #   for(df_name in names(merged_list)){
    #     if(endsWith(df_name, DMA_SettingsInfo[["denominator"]])){
    #       modified_df_list[[df_name]] <- merged_list[[df_name]]
    #     }
    #   }
    #   STAT_C1vC2 <- modified_df_list
    # }
    
    DMA_Output <-  merged_list
    
    for(DF in names(DMA_Output)){
      DMA_Output_Save <- DMA_Output[[DF]]
      DF_save <- gsub("[^A-Za-z0-9._-]", "_", DF)## Remove special characters and replace spaces with underscores
      
      if(is.null(Save_as_Results)==FALSE){
      if (Save_as_Results == "xlsx"){
        xlsDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",DF_save, ".xlsx"))   # Save the DMA results table
        writexl::write_xlsx(DMA_Output_Save,xlsDMA, col_names = TRUE) # save the DMA result DF
      }else if (Save_as_Results == "csv"){
        csvDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",DF_save, ".csv"))
        write.csv(DMA_Output_Save,csvDMA) # save the DMA result DF
      }else if (Save_as_Results == "txt"){
        txtDMA <- file.path(Results_folder_Conditions,paste0("DMA_Output_",DF_save, ".txt"))
        write.table(DMA_Output_Save,txtDMA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
      }
    }
    }
    
    #return(DMA_Output)
    
    if(CoRe==TRUE){
      x <- "Log2(Distance)"
      VolPlot_SettingsInfo= c(color="CoRe")
      VOlPlot_SettingsFile = DMA_Output
    }else{
      x <- "Log2FC"
      VolPlot_SettingsInfo= NULL
      VOlPlot_SettingsFile = NULL
    }
    
    #if (Labels==TRUE){Labels <- NULL}else{Labels <-""}
    
    ################################################################################################################################################################################################
    ###############  Plots ###############
    volplotList = list()
    
    if(is.null(Save_as_Plot)==FALSE){
      for(DF in names(DMA_Output)){ # DF = names(DMA_Output)[1]
        Volplotdata<- DMA_Output[[DF]]
        #Volplotdata<- DMA_Output$`Sicca_Lac-Sicca_Lac_Ab`
        Labels <- Volplotdata[which(Volplotdata$p.adj < PCutoff & abs(Volplotdata$Log2FC) > FCcutoff),]$Metabolite
        
        dev.new()
        VolcanoPlot <- invisible(Volcano(Plot_Settings="Standard",
                                         Input_data =Volplotdata,
                                         Plot_SettingsInfo=VolPlot_SettingsInfo,
                                         Plot_SettingsFile=VOlPlot_SettingsFile[[DF]],
                                         y= "p.adj",
                                         x= x,
                                         AdditionalInput_data= NULL,
                                         OutputPlotName= DF,
                                         Comparison_name= c(DMA_data="Cond1", AdditionalDMA_data= "Cond2"),
                                         xlab= "~Log[2]~FC",
                                         ylab= "~-Log[10]~p.adj",
                                         PCutoff= PCutoff,
                                         FCcutoff= FCcutoff,
                                         color_palette= NULL,
                                         shape_palette=NULL,
                                         SelectLab= Labels,
                                         Connectors=  FALSE,
                                         Subtitle=  paste("Anova test,","p.adj<",PCutoff,"log2 FC >",FCcutoff),
                                         Theme= NULL,
                                         Save_as_Plot= NULL))
        
        volplotList[[DF]]<- VolcanoPlot
        
        DF_save <- gsub("[^A-Za-z0-9._-]", "_", DF)## Remove special characters and replace spaces with underscores
        volcanoDMA <- file.path(Results_folder_Conditions,paste0( "Volcano_Plot_", DF_save,".",Save_as_Plot))
        ggsave(volcanoDMA,plot=VolcanoPlot, width=10, height=8) # save the volcano plot
        dev.off()
      }
    }

    
    # Save results
    # writexl::write_xlsx(filtered_matrix, paste(Results_folder_Preprocessing_folder, "/Feature_Filtered_Matrix.xlsx", sep = ""))#,showNA = TRUE)
    # Make Result list
    AnovaTestRes <- list("DF" = list("DMA_result"=DMA_Output),"Plot"=list("Volcano"=volplotList))
    # Save output RDS
    save(AnovaTestRes, file = file.path(RDS_folder,  OutputName))
    

  }
return(AnovaTestRes)
}
