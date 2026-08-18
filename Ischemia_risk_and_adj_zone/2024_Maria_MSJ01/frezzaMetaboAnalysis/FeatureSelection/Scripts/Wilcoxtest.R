

#STAT_C1vC2 <-MetaProViz:::DMA_Stat_single(C1=C1, C2=C2, Log2FC_table=Log2FC_table, STAT_pval=STAT_pval, STAT_padj=STAT_padj)



#################################data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==#########################################################
### ### ### DMA helper function: Internal Function to perform single comparison ### ### ###
##########################################################################################

#' @param C1 This is the C1 (Condition 1) DF generated within the DMA function.
#' @param C2 This is the C2 (Condition 2) DF generated within the DMA function.
#' @param Log2FC_table this is the Log2FC DF generated within the DMA function.
#' @param Metabolites_Miss these are the metabolites with missing values generated within the DMA function.
#' @param STAT_pval Passed to DMA
#' @param STAT_padj Passed to DMA
#'
#' @keywords DMA helper function
#' @noRd
#'

wTest <- function(Result_Folder = Folder_paths$DMAFolder,
                  RDS_folder = Folder_paths$DMARDS,
                  scriptpath= scriptpath,
                  DMA_data,
                  DMA_SettingsFile,
                  DMA_SettingsInfo = c(conditions="Conditions", numerator = NULL, denumerator = NULL),
                  Log2FC_table, 
                  
                  CoRe=FALSE,
                  STAT_padj = "fdr",
                  Labels = TRUE,
                  
                  Save_as_Results = "xlsx",
                  Save_as_Plot = "svg",
                  
                  FolderName = NULL,
                  OutputName = NULL
){
  
  # Source 
  # usethis::edit_file(file.path(scriptpath,"Volcano.R"))
  source((file.path(scriptpath,"Volcano.R")))
  
  ## 1,5. ------------ Check Inputs------------ ##
  #3. DMA_SettingsInfo
  if(DMA_SettingsInfo[["Conditions"]] %in% DMA_SettingsInfo==TRUE){
    if(DMA_SettingsInfo[["Conditions"]] %in% colnames(DMA_SettingsFile)== FALSE){
      stop("The ",DMA_SettingsInfo[["Conditions"]], " column selected as Conditions in DMA_SettingsInfo was not found in DMA_SettingsFile. Please check your input.")
    }else{# if true rename to Conditions
      DMA_SettingsFile<- DMA_SettingsFile%>%
        dplyr::rename("Conditions"= paste(DMA_SettingsInfo[["Conditions"]]) )
    }
  }else{
    stop("You have to provide a DMA_SettingsInfo for conditions.")
  }
  
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
  
  comparisons <- data.frame(c(numerator, denominator))
  
  
  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)){ OutputName <- "wTestRes.Rdata"
  }else{ OutputName <- paste("wTest_" ,OutputName, ".Rdata" ,sep = "" )}
  if (file.exists(file.path(RDS_folder,  OutputName))==TRUE){
    message("Loading already existing ", paste(OutputName))
    load(file.path(RDS_folder,  OutputName), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running Load data script main body.")

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
    
    
    
    ## ------------ Check Missingness ------------- ##
    #7.
    # If missing value imputation has not been performed the input data will most likely contain NA or 0 values for some metabolites, which will lead to Log2FC = NA.
    # Here we will check how many metabolites this affects in Num and Denom, and weather all replicates of a metabolite are affected.
    Num_Miss <- replace(Num, Num==0, NA)
    Num_Miss <- Num_Miss[, (colSums(is.na(Num_Miss)) > 0), drop = FALSE]
    
    Denom_Miss <- replace(Denom, Denom==0, NA)
    Denom_Miss <- Denom_Miss[, (colSums(is.na(Denom_Miss)) > 0), drop = FALSE]
    
    if((ncol(Num_Miss)>0 & ncol(Denom_Miss)==0)){
      message("In `numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
      Metabolites_Miss <- colnames(Num_Miss)
    } else if(ncol(Num_Miss)==0 & ncol(Denom_Miss)>0){
      message("In `denominator` ",paste0(toString(denominator)), ", NA/0 values exist in ", ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "), ". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
      Metabolites_Miss <- colnames(Denom_Miss)
    } else if(ncol(Num_Miss)>0 & ncol(Denom_Miss)>0){
      message("In `numerator` ",paste0(toString(numerator)), ", NA/0 values exist in ", ncol(Num_Miss), " Metabolite(s): ", paste0(colnames(Num_Miss), collapse = ", "), " and in `denominator`",paste0(toString(denominator)), " ",ncol(Denom_Miss), " Metabolite(s): ", paste0(colnames(Denom_Miss), collapse = ", "),". Those metabolite(s) will return p.val= NA, p.adj.= NA, t.val= NA. The Log2FC = Inf, if all replicates are 0/NA.")
      Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
      Metabolites_Miss <- unique(Metabolites_Miss)
    } else{
      message("There are no NA/0 values")
      Metabolites_Miss <- c(colnames(Num_Miss), colnames(Denom_Miss))
      Metabolites_Miss <- unique(Metabolites_Miss)
    }
    
    
    C1 <- DMA_data %>% # Numerator
      filter(DMA_SettingsFile$Conditions %in% comparisons[1,1]) %>%
      select_if(is.numeric)#only keep numeric columns with metabolite values
    C2 <- DMA_data %>% # Deniminator
      filter(DMA_SettingsFile$Conditions %in%  comparisons[2,1]) %>%
      select_if(is.numeric)
    
    
    ## ------------ Perform Hypothesis testing ----------- ##
    # For C1 and C2 we use 0, since otherwise we can not perform the statistical testing.
    C1[is.na(C1)] <- 0
    C2[is.na(C2)] <- 0
    
    # these are for fabri
    # DMA_SettingsFile <- DMA_SettingsFile[order(DMA_SettingsFile$Biological_Replicates, DMA_SettingsFile$Conditions),]
    # 
    # cancerOrder <- DMA_SettingsFile %>% filter(Conditions!="Healthy")%>% rownames()
    # cancerOrder <-  cancerOrder[-7]
    # healthyOrder <- DMA_SettingsFile %>% filter(Conditions=="Healthy")%>% rownames()
    # 
    # rownames(C1) %in% cancerOrder
    # 
    # 
    # C1 <- C1[match(cancerOrder, rownames(C1)), ]
    # C2 <- C2[match(healthyOrder, rownames(C2)), ]
    # C1 <- C1[complete.cases(C1),]
    # C2 <- C2[complete.cases(C2),]
    # 
    # C1set <-DMA_SettingsFile %>% filter(Conditions!="Healthy")%>% rownames()
    # C2set <-DMA_SettingsFile %>% filter(Conditions=="Healthy")

    
    #### 1. p.value and test statistics (=t.val)
    T_C1vC2 <-mapply(wilcox.test, x= as.data.frame(C2), y = as.data.frame(C1), exact=FALSE, SIMPLIFY = F)#, paired = TRUE)
    
    VecPVAL_C1vC2 <- c()
    VecTVAL_C1vC2 <- c()
    for(i in 1:length(T_C1vC2)){
      p_value <- unlist(T_C1vC2[[i]][3])
      t_value <- unlist(T_C1vC2[[i]])[1]   # Extract the t-value
      VecPVAL_C1vC2[i] <- p_value
      VecTVAL_C1vC2[i] <- t_value
    }
    Metabolite <- colnames(C2)
    PVal_C1vC2 <- data.frame(Metabolite, p.val = VecPVAL_C1vC2, t.val = VecTVAL_C1vC2)
    
    
    #### 2. p.adjusted
    #Split data for p.value adjustment to exclude NA
    PVal_NA <- PVal_C1vC2[is.na(PVal_C1vC2$p.val), c(1:3)]
    PVal_C1vC2 <-PVal_C1vC2[!is.na(PVal_C1vC2$p.val), c(1:3)]
    
    #perform adjustment
    VecPADJ_C1vC2 <- p.adjust((PVal_C1vC2[,2]),method = STAT_padj, n = length((PVal_C1vC2[,2]))) #p-adjusted
    Metabolite <- PVal_C1vC2[,1]
    PADJ_C1vC2 <- data.frame(Metabolite, p.adj = VecPADJ_C1vC2)
    STAT_C1vC2 <- merge(PVal_C1vC2,PADJ_C1vC2, by="Metabolite", sort = FALSE)
    
    #Add Metabolites that have p.val=NA back into the DF for completeness.
    if(nrow(PVal_NA)>0){
      PVal_NA$p.adj <- NA
      STAT_C1vC2 <- rbind(STAT_C1vC2, PVal_NA)
    }
    
    #Add Log2FC
    STAT_C1vC2 <- merge(Log2FC_table,STAT_C1vC2[,c(1:2,4,3)], by="Metabolite", sort = FALSE)
    
    #order for t.value
    STAT_C1vC2 <- STAT_C1vC2[order(STAT_C1vC2$t.val,decreasing=TRUE),] # order the df based on the t-value
    
    colnames(STAT_C1vC2)[2] <- "Log2FC"
    Output <- STAT_C1vC2
    
    
    # DMA_Output <- merge(savedMetaboliteNames, DMA_Output, by="Metabolite", sort = FALSE)
    # DMA_Output$Metabolite <- NULL
    # colnames(DMA_Output)[1] <- "Metabolite"
    
    
    if (Save_as_Results == "xlsx"){
      xlsDMA <- file.path(Result_Folder,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator), OutputName,".xlsx"))   # Save the DMA results table
      writexl::write_xlsx(Output,xlsDMA, col_names = TRUE) # save the DMA result DF
    }else if (Save_as_Results == "csv"){
      csvDMA <- file.path(Result_Folder,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator),OutputName,".csv"))
      write.csv(Output,csvDMA) # save the DMA result DF
    }else if (Save_as_Results == "txt"){
      txtDMA <- file.path(Result_Folder,paste0("DMA_Output_",toString(numerator),"_vs_",toString(denominator), OutputName,".txt"))
      write.table(Output,txtDMA, col.names = TRUE, row.names = FALSE) # save the DMA result DF
    }
    
    DMA_Output <- Output
    
    ## Make the volcano plot
    if(CoRe==TRUE){
      x <- "Log2(Distance)"
      VolPlot_SettingsInfo= c(color="CoRe")
      VOlPlot_SettingsFile = DMA_Output
    }else{
      x <- "Log2FC"
      VolPlot_SettingsInfo= NULL
      VOlPlot_SettingsFile = NULL
    }
    
    if (Labels==TRUE){Labels <- NULL}else{Labels <-""}
    
    ###############  Plots ###############
    volplotList = list()
    # Make a simple Volcano plot
    dev.new()
    VolcanoPlot <- invisible(Volcano(Plot_Settings="Standard",
                                        Input_data=DMA_Output,
                                        Plot_SettingsInfo=VolPlot_SettingsInfo,
                                        Plot_SettingsFile=VOlPlot_SettingsFile,
                                        y= "p.adj",
                                        x= x,
                                        AdditionalInput_data= NULL,
                                        OutputPlotName= paste0(toString(numerator)," versus ",toString(denominator)),
                                        Comparison_name= c(Input_data="Cond1", AdditionalInput_data= "Cond2"),
                                        xlab= NULL,#"~Log[2]~FC"
                                        ylab= NULL,#"~-Log[10]~p.adj"
                                      #  pCutoff= 0.05,
                                        FCcutoff= 1,
                                        color_palette= NULL,
                                        shape_palette=NULL,
                                        SelectLab= Labels,
                                        Connectors=  FALSE,
                                        Subtitle=  bquote(italic("Differential Metabolite Analysis, w-test")),
                                        Theme= NULL,
                                        Save_as_Plot= NULL))
    volplotList[[paste0(toString(numerator)," versus ",toString(denominator))]]<- VolcanoPlot
    dev.off()
    #plot(VolcanoPlot)
    
    volcanoDMA <- file.path(Result_Folder,paste0( "Volcano_Plot_",toString(numerator),"_versus_",toString(denominator),OutputName,".",Save_as_Plot))
    ggsave(volcanoDMA,plot=VolcanoPlot, width=10, height=8) # save the volcano plot
    
    
    # Save results
    # writexl::write_xlsx(filtered_matrix, paste(Results_folder_Preprocessing_folder, "/Feature_Filtered_Matrix.xlsx", sep = ""))#,showNA = TRUE)
    # Make Result list
    suppressWarnings(wTestRes <- list("DF" = list("DMA_result"=DMA_Output),"Plot"=list("Volcano"=volplotList)))
    # Save output RDS
    save(wTestRes, file = file.path(RDS_folder,  OutputName))
  }
  return(wTestRes)
}

