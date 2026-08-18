
# @param data Is the output of CD
# @param mode = ("pos", "neg", "both")


# Load BinneR function
BinneR  <- function(Result_Folder = Folder_paths$BinnerFolder,
                    RDS_folder = Folder_paths$RDSFolder,
                    OutputName = NULL,
                    ForceRun = FALSE,
                    data, 
                    CorrelationThreshold=0.5 ){
  
  
  ## 3. ------------ Check outputs ----------- ##
  if (is.null(OutputName)) {OutputName = ""} else{OutputName = paste("_",OutputName, sep="")}
  if (file.exists(file.path(RDS_folder, paste("BinnerRes", OutputName,".Rdata", sep = ""))) == TRUE & ForceRun == FALSE) { # Check if RDS exists and forceRun is false then skip script and load result from RDS, otherwise run it
    message("Loading already existing ", paste("BinnerRes", OutputName,".Rdata", sep = ""))
    load(file.path(RDS_folder, paste("BinnerRes", OutputName,".Rdata", sep = "")), .GlobalEnv)
    
    ## 4. ------------ Do calculations ------------- ##
  }else{
    message("Running BinneR script main body.")
    
    if(grepl("neg", OutputName)==T){
      mode = "neg"
    }else{
      mode="pos"
    }
    
  # Load Libraries
  library(tidyverse)
  suppressMessages(library(Hmisc)) # correlation
  suppressWarnings(library(cluster)) # clistering
  
  # Parse data
  data <- data %>% select("Compound", everything())
  
  # Sort the data by RT column
  data <- data %>% arrange( RT)
  data$Compound <- as.character(data$Compound)
  
  # Make output
  Output <- data %>% select(1:3)
  Output$Bin <- NA
  Output$Annotation <- NA
  # 5. Bin features by RT column
  rt_gap <- 0.04
  bins <- list()
  current_bin <- list()
  previous_rt <- data$RT[1]
  
  for (i in 1:nrow(data)) {
    current_rt <- data$RT[i]
    if ((current_rt - previous_rt) > rt_gap) {
      bins <- append(bins, list(current_bin))
      current_bin <- list()
    }
    current_bin <- append(current_bin, list(data[i, ]))
    previous_rt <- current_rt
  }
  bins <- append(bins, list(current_bin))
  
  for (k in 1:length(bins)){ # k = 3
    bin <- bins[[k]]
    for (l in 1:length(bin)){ # l = 1
      feature <- bin[[l]]
      
      if(feature$Compound %in% Output$Compound){
        Output[Output$Compound %in%   feature$Compound,"Bin"] <- k
      }
    }
  }
  
  
  # 6. Function to calculate Spearman correlation matrix for a bin
  calculate_spearman_correlation <- function(bin) {
    bin_df <- do.call(rbind, bin)
    if (nrow(bin_df) < 2) {
     # print("Too few samples")
      return(NULL)
    }else{
     # print("Works")
      sample_data <- bin_df[,-c(1:3)]
      sample_data[] <- lapply(sample_data, function(x) as.numeric(as.character(x)))
      cor_matrix <- rcorr(as.matrix(t(sample_data)), type = "spearman")
      cor_matrix <- as.data.frame(cor_matrix$r)
      cor_matrix[is.na(cor_matrix)] <- 0
      rownames(cor_matrix) <- bin_df$Compound
      colnames(cor_matrix) <- bin_df$Compound
      return(cor_matrix) 
    }
  }
  
  # 8. Calculate Spearman correlation coefficients and cluster features for each bin
  results <- lapply(bins, function(bin) { # bin = bins[[1]]
    cor_matrix <- calculate_spearman_correlation(bin)
  })
  
  
  # 10. Load annotation file
  annotation <- read.csv(file.path(scriptpath,"Functions/annotationfileFrezza.txt"), sep = "\t")
  
  if(mode == "pos"){
    annotation <- annotation %>% filter(Mode != "Negative")
    sign <- "+"
  }else if (mode == "neg"){
    annotation <- annotation %>% filter(Mode != "Positive")
    sign <- "-"
  }
  
  # count if a match was found
  match=0
  matchUpdated=0
  # to update annotation
  M=0
  # 11. Identify most abundant feature in each bin and test adducts
  for (result in 1:length(results)){ # result =4
    
    # for each bin
    bin =results[[result]]
    fullCorrTable =results[[result]]
    bindata <- data[data$Compound %in%  rownames(bin),]
    # while the bin has more than 1 feature
    if (is.null(bindata[[1]])!=T ){
      while(dim(bindata)[1]>0){
        
        # Get highest intensity feature 
        max <- which.max(rowMeans(bindata %>% select(-c(1,2,3)), na.rm = T) ) %>% names()# rowname of max in cluster
        parent <- bindata[ rownames(bindata)==max,"Compound"]
        parentMZ <- bindata[rownames(bindata)==max,"mz"]
        
        if(mode == "pos"){
          parentMW = parentMZ - 1.007276
        }else{
          parentMW = parentMZ + 1.007276
        }
        # for each starting M
        for (X in c(1,2,3,4)){# X=1  for (X in c(1,2,3,4,5,6,7,8)
          parentMWcheck <- X*as.numeric(parentMW)
          
          # test each adduct in annotation list
          for (adduct in 1:length(annotation$Annotation)){ #  adduct <- 1
            adduct <- annotation[adduct,]
            featurewithadduct <-  parentMWcheck+adduct$Mass
            
            # check if the adduct is found in data 
            if(sum(abs(bindata %>% filter(mz != parentMZ) %>% select(mz) %>% as.vector()%>% unlist() - featurewithadduct) < 0.005)>0){# mass error allowed
              
              #abs(bindata %>% filter(mz != parentMZ) %>% select(mz) %>% as.vector()%>% unlist() - featurewithadduct) < 0.005
              
              parentlessbindata <- bindata[bindata$mz!=parentMZ,]
              child <- parentlessbindata[abs(bindata %>% filter(mz != parentMZ) %>% select(mz) %>% as.vector()%>% unlist() - featurewithadduct) < 0.005,"Compound"]
              for(kid in child){
                # if found check correlation corr >0.6 then add the annotation
                if(bin[rownames(bin) %in% kid,colnames(bin) %in% parent]>CorrelationThreshold){
                  message("Match found")
                  # parent annot
                  Output[Output$Compound %in% parent,"Annotation"] <- paste0("M",M)
                  # child annot
                  if(X==1){
                    Output[Output$Compound %in% kid,"Annotation"] <- paste0("M",M,sign,adduct$Annotation )
                  }else{
                    Output[Output$Compound %in% kid,"Annotation"] <- paste0(X,"M",M,sign,adduct$Annotation )
                  }
                  
                  # remove child from bindata
                  bindata <- bindata %>% filter(Compound != kid)
                  # update match
                  matchUpdated = match+1
                }
              }
            } # end of checking mzs
          } # end of testing adducts
          
          # here test for double adducts
          # test each adduct in annotation list
          for (add1 in 1:length(annotation$Annotation)){ #  add1 <- 3
            for (add2 in 1:length(annotation$Annotation)){ # add2 <- 6
              if(add1 != add2){
                adduct1 <- annotation[add1,]
                adduct2 <- annotation[add2,]
                featurewithadduct <-  parentMWcheck+adduct1$Mass +adduct2$Mass
                
                if(adduct1$Mass>0){
                  sign1 = "+"
                }else(
                  sign1 = "-"
                )
                
                if(adduct2$Mass>0){
                  sign2 = "+"
                }else(
                  sign2 = "-"
                )
                
                # if the feature + adduct was found check correlation score
                if(sum(abs(bindata %>% filter(mz != parentMZ) %>% select(mz) %>% as.vector()%>% unlist() - featurewithadduct) < 0.005)>0){# mass error allowed
                  
                  parentlessbindata <- bindata[bindata$mz!=parentMZ,]
                  child <- parentlessbindata[abs(bindata %>% filter(mz != parentMZ) %>% select(mz) %>% as.vector()%>% unlist() - featurewithadduct) < 0.005,"Compound"]
                  for(kid in child){
                    # if corr >0.6 then add the annotation
                    if(bin[rownames(bin) %in% kid,colnames(bin) %in% parent]>CorrelationThreshold){
                      message("Match found")
                      # parent annot
                      Output[Output$Compound %in% parent,"Annotation"] <- paste0("M",M)
                      if(X==1){
                        # child annot
                        Output[Output$Compound %in% kid,"Annotation"] <- paste0("M",M,sign1,adduct1$Annotation,sign2,adduct2$Annotation )
                      }else{
                        # child annot
                        Output[Output$Compound %in% kid,"Annotation"] <- paste0(X,"M",M,sign1,adduct1$Annotation,sign2,adduct2$Annotation2 )
                      }
                      
                      # remove child from bindata
                      bindata <- bindata %>% filter(Compound != kid)
                      # update match
                      matchUpdated = match+1
                    }
                  }
                
                } # end of checking mzs
              } #  end of if(adduct1 != adduct2)
            }
          } # end of testing double adducts
          
          
        } # end of X x M
        
        # remove parent from data
        bindata <- bindata[!bindata$Compound %in% parent,]
        if (matchUpdated >match){
          M=M+1   
          match <- matchUpdated
        }
        
      } # end of while
    } # end of if >1 obs in bin
    print(paste0("Result: ", result, "/", length(results)))
  }
  
  # save result
  writexl::write_xlsx(Output, paste0(Result_Folder, "/Binner_result", OutputName,".xlsx", sep = ""))#,showNA = TRUE)
  
  # Save the clustering results
  save(Output, file = file.path(RDS_folder, paste("BinnerRes",OutputName,".Rdata", sep = "")))
  
  
  }
  return(Output)
}
