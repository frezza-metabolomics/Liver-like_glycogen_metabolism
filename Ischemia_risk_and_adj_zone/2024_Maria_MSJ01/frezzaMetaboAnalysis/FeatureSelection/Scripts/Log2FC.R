#' This script allows you to perform differential metabolite analysis to obtain a Log2FC, pval, padj and tval comparing two or multiple conditions.
#'
#' @param DMA_data DF with unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param DMA_SettingsFile DF which contains metadata information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames.
#' @param DMA_SettingsInfo \emph{Optional: } Named vector including the information about the conditions column c(conditions="ColumnName_Plot_SettingsFile"). Can additionally pass information on numerator or denominator c(numerator = "ColumnName_Plot_SettingsFile", denumerator = "ColumnName_Plot_SettingsFile") for specifying which comparison(s) will be done (one-vs-one, all-vs-one, all-vs-all). Using =NULL selects all the condition and performs multiple comparison all-vs-all. Log2FC are obtained by dividing the numerator by the denominator, thus positive Log2FC values mean higher expression in the numerator and are presented in the right side on the Volcano plot (For CoRe the Log2Distance). \strong{Default = c(conditions="Conditions", numerator = NULL, denumerator = NULL)}
#' @param STAT_pval \emph{Optional: } String which contains an abbreviation of the selected test to calculate p.value. For one-vs-one comparisons choose t.test, wilcox.test, "chisq.test" or "cor.test", for one-vs-all or all-vs-all comparison choose aov (=annova), kruskal.test or lmFit (=limma) \strong{Default = "t-test"}
#' @param STAT_padj \emph{Optional: } String which contains an abbreviation of the selected p.adjusted test for p.value correction for multiple Hypothesis testing. Search: ?p.adjust for more methods:"BH", "fdr", "bonferroni", "holm", etc.\strong{Default = "fdr"}
#' @param OutputName String which is added to the output files of the DMA.
#' @param DMA_MetaFile_Metab \emph{Optional: } DF which contains the metadata information , i.e. pathway information, retention time,..., for each metabolite. \strong{Default = NULL}
#' @param CoRe \emph{Optional: } TRUE or FALSE for whether a Consumption/Release  input is used \strong{Default = FALSE}
#' @param Save_as_Plot \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf. \strong{Default = svg}
#' @param Save_as_Results \emph{Optional: } File types for the analysis results are: "csv", "xlsx", "txt" \strong{Default = "csv"}
#' @param plot \emph{Optional: } TRUE or FALSE, if TRUE Volcano plot is saved as an overview of the results. \strong{Default = TRUE}
#' @param FolderName {Optional:} String which is added to the resulting folder name \strong(Default = NULL)
#'
#' @keywords Differential Metabolite Analysis, Multiple Hypothesis testing, Normality testing
#' @export


########################################################
### ### ### Differential Metabolite Analysis ### ### ###
########################################################

Log2FC <-function(Result_Folder,
                  RDS_folder,
                  scriptpath,
                  DMA_data,
                  DMA_SettingsFile,
                  DMA_SettingsInfo = c(conditions="Conditions", numerator = NULL, denumerator = NULL),
                  ForceRun = FALSE,
                  CoRe=FALSE,
                  Plot = TRUE,
                  Save_as_Results = "csv",
                  Save_as_Plot = "svg",
                 
                  FolderName = NULL,
                  OutputName = NULL
){
  

  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)){ OutputName <- "Log2FCRes.Rdata"
  }else{ OutputName <- paste("Log2FCRes_" ,OutputName, ".Rdata" ,sep = "" )}
  if (file.exists(file.path(RDS_folder,  OutputName))==TRUE & ForceRun == FALSE){
    message("Loading already existing ", paste(OutputName))
    load(file.path(RDS_folder,  OutputName), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running Load data script main body.")
    
    
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
      stop("You have to provide a DMA_SettingsInfo for conditions.")
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
      conditions = DMA_SettingsFile$Conditions
      denominator <-unique(DMA_SettingsFile$Conditions)
      numerator <-unique(DMA_SettingsFile$Conditions)
      comparisons <- combn(unique(conditions), 2) %>% as.matrix()
      #Settings:
      MultipleComparison = TRUE
      all_vs_all = TRUE
    }else if("denominator" %in% names(DMA_SettingsInfo)==TRUE  & "numerator" %in% names(DMA_SettingsInfo)==FALSE){
      #all-vs-one: Generate the pairwise combinations
      conditions = DMA_SettingsFile$Conditions
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
      comparisons <- matrix(c(denominator, numerator))
      #Settings:
      MultipleComparison = FALSE
      all_vs_all = FALSE
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
    
    
    Log2FC_table <- data.frame(Metabolite = colnames(DMA_data))
    for (column in 1:dim(comparisons)[2]){
      
      message(paste(column, " out of ", dim(comparisons)[2], sep=""))
      
      C1 <- DMA_data %>% # Numerator
        filter(DMA_SettingsFile$Conditions %in% comparisons[1,column]) %>%
        select_if(is.numeric)#only keep numeric columns with metabolite values
      C2 <- DMA_data %>% # Deniminator
        filter(DMA_SettingsFile$Conditions %in%  comparisons[2,column]) %>%
        select_if(is.numeric)
      
      ## ------------  Calculate Log2FC ----------- ##
      # For C1_Mean and C2_Mean use 0 to obtain values, leading to Log2FC=NA if mean = 0 (If one value is NA, the mean will be NA even though all other values are available.)
      C1_Zero <- C1
      C1_Zero[is.na(C1_Zero)] <- 0
      Mean_C1 <- C1_Zero %>%
        summarise_all("mean")
      
      C2_Zero <- C2
      C2_Zero[is.na(C2_Zero)] <- 0
      Mean_C2 <- C2_Zero %>%
        summarise_all("mean")
      
      if(CoRe==TRUE){#Calculate absolute distance between the means. log2 transform and add sign (-/+):
        #CoRe values can be negative and positive, which can does not allow us to calculate a Log2FC.
        Mean_C1_t <- as.data.frame(t(Mean_C1))%>%
          rownames_to_column("Metabolite")
        Mean_C2_t <- as.data.frame(t(Mean_C2))%>%
          rownames_to_column("Metabolite")
        Mean_Merge <-merge(Mean_C1_t, Mean_C2_t, by="Metabolite", all=TRUE, sort = FALSE)%>%
          rename("C1"=2,
                 "C2"=3)
        
        #Deal with NA/0s
        # Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss#Column to enable the check if mean values of 0 are due to missing values (NA/0) and not by coincidence
        
        if(any((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))==TRUE){
          Mean_Merge <- Mean_Merge%>%
            mutate(C1 = case_when(C2 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                  C1 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                  C2 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                  C1 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                  TRUE ~ paste(C1)))%>%
            mutate(C2 = case_when(C1 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                  C2 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                  C1 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                  C2 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                  TRUE ~ paste(C2)))%>%
            mutate(C1 = as.numeric(C1))%>%
            mutate(C2 = as.numeric(C2))
          
          X <- Mean_Merge%>%
            subset((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))
          message("We added +1 to the mean value of metabolite(s) ", paste0(X$Metabolite, collapse = ", "), ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
        }
        
        #Add the distance column:
        Mean_Merge$`Log2(Distance)` <-log2(abs(Mean_Merge$C1 - Mean_Merge$C2))
        
        Mean_Merge <- Mean_Merge%>%#Now we can adapt the values to take into account the distance
          mutate(`Log2(Distance)` = case_when(C1 > C2 ~ paste(`Log2(Distance)`*+1),#If C1>C2 the distance stays positive to reflect that C1 > C2
                                              C1 < C2 ~ paste(`Log2(Distance)`*-1),#If C1<C2 the distance gets a negative sign to reflect that C1 < C2
                                              TRUE ~ 'NA'))%>%
          mutate(`Log2(Distance)` = as.numeric(`Log2(Distance)`))
        
        
        #Log2FC_table <-Mean_Merge[,c(1,5)]
        # new down here
        #Add additional information:
        temp1 <- Mean_C1
        temp2 <- Mean_C2
        #Add Info of CoRe:
        CoRe_info <- rbind(temp1, temp2,rep(0,length(temp1)))
        for (i in 1:length(temp1)){
          if (temp1[i]>0 & temp2[i]>0){
            CoRe_info[3,i] <- "Released"
          }else if (temp1[i]<0 & temp2[i]<0){
            CoRe_info[3,i] <- "Consumed"
          }else if(temp1[i]>0 & temp2[i]<0){
            CoRe_info[3,i] <- paste("Released in" ,comparisons[1,column] , "and Consumed",comparisons[2,column] , sep=" ")
          } else if(temp1[i]<0 & temp2[i]>0){
            CoRe_info[3,i] <- paste("Consumed in" ,comparisons[1,column] , " and Released",comparisons[2,column] , sep=" ")
          }else{
            CoRe_info[3,i] <- "No Change"
          }
        }
        
        CoRe_info <- t(CoRe_info) %>% as.data.frame()
        CoRe_info <- rownames_to_column(CoRe_info, "Metabolite")
        names(CoRe_info)[2] <- paste("Mean",  comparisons[1,column], sep="_")
        names(CoRe_info)[3] <- paste("Mean",  comparisons[2,column], sep="_")
        names(CoRe_info)[4] <- "CoRe_specific"
        
        CoRe_info <-CoRe_info%>%
          mutate(CoRe = case_when(CoRe_specific == "Released" ~ 'Released',
                                  CoRe_specific == "Consumed" ~ 'Consumed',
                                  TRUE ~ 'Released/Consumed'))
        
        Log2FC_C1vC2 <-merge(Mean_Merge[,c(1,4)], CoRe_info[,c(1,4:5,2:3)], by="Metabolite", all.x=TRUE)
        Log2FC_table <- Log2FC_C1vC2
        
      }else if(CoRe==FALSE){
        #Mean values could be 0, which can not be used to calculate a Log2FC and hence the Log2FC(A versus B)=(log2(A+x)-log2(B+x)) for A and/or B being 0, with x being set to 1
        Mean_C1_t <- as.data.frame(t(Mean_C1))%>%
          rownames_to_column("Metabolite")
        Mean_C2_t <- as.data.frame(t(Mean_C2))%>%
          rownames_to_column("Metabolite")
        Mean_Merge <-merge(Mean_C1_t, Mean_C2_t, by="Metabolite", all=TRUE, sort = FALSE)%>%
          rename("C1"=2,
                 "C2"=3)
        Mean_Merge$`NA/0` <- Mean_Merge$Metabolite %in% Metabolites_Miss#Column to enable the check if mean values of 0 are due to missing values (NA/0) and not by coincidence
        
        Mean_Merge <- Mean_Merge%>%
          mutate(C1_Adapted = case_when(C2 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                        C1 == 0 & `NA/0`== TRUE ~ paste(C1),#Here we have a "true" 0 value due to 0/NAs in the input data
                                        C2 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                        C1 == 0 & `NA/0`== FALSE ~ paste(C1+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                        TRUE ~ paste(C1)))%>%
          mutate(C2_Adapted = case_when(C1 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                        C2 == 0 & `NA/0`== TRUE ~ paste(C2),#Here we have a "true" 0 value due to 0/NAs in the input data
                                        C1 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                        C2 == 0 & `NA/0`== FALSE ~ paste(C2+1),#Here we have a "false" 0 value that occured at random and not due to 0/NAs in the input data, hence we add the constant +1
                                        TRUE ~ paste(C2)))%>%
          mutate(C1_Adapted = as.numeric(C1_Adapted))%>%
          mutate(C2_Adapted = as.numeric(C2_Adapted))
        
        if(any((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))==TRUE){
          X <- Mean_Merge%>%
            subset((Mean_Merge$`NA/0`==FALSE & Mean_Merge$C1 ==0) | (Mean_Merge$`NA/0`==FALSE & Mean_Merge$C2==0))
          message("We added +1 to the mean value of metabolite(s) ", paste0(X$Metabolite, collapse = ", "), ", since the mean of the replicate values where 0. This was not due to missing values (NA/0).")
        }
        
        #Calculate the Log2FC
        Mean_Merge$FC_C1vC2 <- Mean_Merge$C1_Adapted/Mean_Merge$C2_Adapted #FoldChange
        Mean_Merge$Log2FC <- gtools::foldchange2logratio(Mean_Merge$FC_C1vC2, base=2)
        Log2FC_C1vC2 <-Mean_Merge[,c(1,8)]
        Log2FC_C1vC2$Log2FC <- as.numeric( Log2FC_C1vC2$Log2FC)
        colnames(Log2FC_C1vC2)[2] <-   paste("Log2FC", paste( comparisons[1,column], comparisons[2,column],sep="-"), sep="_")

        Log2FC_table <- merge(Log2FC_table, Log2FC_C1vC2, by= "Metabolite", sort = FALSE)

      }else{
        stop("Please choose CoRe= TRUE or CoRe=FALSE.")
      }
    }
    
    
    
    # Save results
    # Save nothing writexl::write_xlsx(filtered_matrix, paste(Results_folder_Preprocessing_folder, "/Feature_Filtered_Matrix.xlsx", sep = ""))#,showNA = TRUE)
    # Make Result list
    Log2FCRes = list("Log2FC_table"= Log2FC_table)
    # Save output RDS
    save(Log2FCRes, file = file.path(RDS_folder,  OutputName))
  }
  return(invisible(Log2FCRes))
}
