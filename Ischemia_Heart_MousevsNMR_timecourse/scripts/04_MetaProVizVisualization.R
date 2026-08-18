#' Metabolomics pre-processing pipeline
#' @author Prymidis Dimitrios, Schmidt Christina
#' Date: "2022-10-28"
#'
#' This script allows you to perform differential metabolite analysis
#'
#' @param Input_data1 Data matrix which contains samples in rows and metabolites, Log fold changes, pvalus and padjustee values in columns
#' @param Input_data2 Data matrix dame as data1 for another comparison
#' @param Output_Name String which contains the name of the output file of the Metabolic Clusters
#' @param Condition1 String which contains the name of the first condition
#' @param Condition2 String which contains the name of the second condition
#' @param pCutoff Number of the desired p value cutoff for assessing significance
#' @param FCcutoff Number of the desired log fold change cutoff for assessing significance
#' @param test String which selects pvalue or padj for significance
#'
#' @export
#'
#'
# Load libraries
library(tidyverse) # general scripting
library(ggfortify) # For visualization PCA
#library(reshape2) # for melting the data
#library(EnhancedVolcano) # for Volcano plots
#library(kableExtra) # used for table
#library(alluvial) # For alluvial plots
#library(ggpubr) # For general plots for adding statistics
#library(ggbeeswarm) # for superplots

#################################
### ### ### PCA Plots ### ### ###
#################################
#Notes
#  The x=0 and y=0 black lines in PCA will always be there regardless the theme change. I cannot yet make it to be there by default and not be there when you change theme.
# Well It can be done but it requires a lot additional of work. So for now the lines will be there
# Palette changing is still missing
# To do: select a better palette and add option to the user to change the palette to whatever they like
#data=Processed_data
MetaProVizPCA <- function(data= data,
                          Design= NULL,
                          Color = FALSE,
                          Shape = FALSE,
                          Fill = FALSE,
                          Show_Loadings = FALSE,
                          Scaling= TRUE,
                          loadingsSCALE = 5, # for the loading to properly ploted
                          # Palette= "Set2".
                          k=1,  # plot the PCk, PCk+1
                          Theme=theme_classic(), # theme_bw()
                          OutputPlotName= 'PCA_plot',
                          Save_as = svg){  #Save_as = "svg"
  
  data <- select_if(as.data.frame(data), is.numeric)
  # Make arguments into strings
  # Color<- deparse(substitute(Color)) # make the input variables into strings
  # Shape<- deparse(substitute(Shape))
  # 
  
  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create PCA plots folder in  result directory ###
  Results_folder_plots_PCA_folder = paste(Results_folder,"/PCA_Plots",  sep="")
  if (!dir.exists(Results_folder_plots_PCA_folder)) {dir.create(Results_folder_plots_PCA_folder)}
  
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  Save_as= deparse(substitute(Save_as))
  
  ### select plot based on arguments
  #1
  if (Fill != FALSE){
    #data=Processed_data
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    pc_data <- prcomp(as.matrix(select_if(mdata, is.numeric)),scale. = TRUE)
    
    #add the Color, Shape, Fill
    #Color
    data <- pc_data$x
    data<- merge(data, mdata %>% select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- column_to_rownames(data, "Row.names")
    #Shape
    data<- merge(as.data.frame(data), mdata %>%select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    data[,Shape] <- data[,Shape] # make the shape into a factor to be discrete
    #Fill <- "Conditions"
    data<- merge(as.data.frame(data), mdata %>%select(Fill), by=0) # merge the data with the design to get only the kept samples
    data <- column_to_rownames(data, "Row.names")
    #Define the type of Color
    # For color if we have character we go with discrete. If we have numeric we to discrete until 4groupd. if we have more we go for continuous
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    #find if you need more than 6 shapes
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    
    #this is a small function to get the asked PCAs
    #k=1
    b <- k
    c <- k+1
    #print(colnames(data))
    colnames(data)[1:2] <- c("s1","s2")
    #print(colnames(data))
    colnames(data)[b:c] <- c("PC1","PC2")
    #print(colnames(data))
	PC1val <- summary(pc_data)$importance[2,][[b]]*100
	PC2val <- summary(pc_data)$importance[2,][[c]]*100
    PCA <- ggplot(data,aes(PC1,PC2, color =  data[[Color]], shape = data[[Shape]])) + 
      geom_point(size=3)+
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)+
      ggtitle(paste(OutputPlotName)) +
      xlab(paste0("PC",b, "(",PC1val,"%)")) +
      ylab(paste0("PC",c,"(",PC2val,"%)")) 
    PCA <- PCA + ggConvexHull::geom_convexhull(aes(PC1,PC2, fill=data[[Fill]]),
                                               alpha = 0.3, inherit.aes = F)+ 
      labs(color = Color, shape = Shape, fill=Fill)   +
    ggrepel::geom_text_repel(aes(label = row.names(data)))
    #PCA
    loading_data <- pc_data
    
    # PCab <- paste0("PC",b)
    # PCac <- paste0("PC",c)
    # ggplot(as.data.frame(pc_data$rotation[,b:c]), aes(x= PCab, y= PCac))+geom_point()
    # scale.default
    #Scale <- 3
    # Create Biplot
    scores<- data.frame(pc_data$x)
    scores[] <- lapply(scores, function(x) x / sqrt(sum((x - mean(x))^2)))
    #loadings <- as.data.frame(pc_data$rotation)
    
    loadings <- select_if(as.data.frame(pc_data$rotation), is.numeric)
    #loadings$Species <- rownames(loadings)
    
    #print(colnames(loadings))
    colnames(loadings)[1:2] <- c("s1","s2")
    #print(colnames(loadings))
    colnames(loadings)[b:c] <- c("PC1","PC2")
    #print(colnames(loadings))
    #print(colnames(scores))
    colnames(scores)[1:2] <- c("s1","s2")
    #print(colnames(scores))
    colnames(scores)[b:c] <- c("PC1","PC2")
    #print(colnames(scores))
    # library(ggplot2)
    scale <- min(max(abs(scores$PC1))/max(abs(loadings$PC1)),
                 max(abs(scores$PC2))/max(abs(loadings$PC2))) * 10
    Biplot <- ggplot(data = data[,c("PC1","PC2")], aes(x = PC1, y = PC2)) +
      geom_point(aes(color = data[[Color]], group = data[[Shape]]), size = 2, shape = 19) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)+
      ggtitle(paste(OutputPlotName))
    Biplot <- Biplot + ggConvexHull::geom_convexhull(aes(PC1,PC2, fill=data[[Fill]]),
                                    alpha = 0.3, inherit.aes = F)+ 
      labs(color = Color, shape = Shape, fill=Fill)
    #loadingsSCALE = 5
    Biplot <- Biplot + geom_segment(
      data = loadings*scale*loadingsSCALE, aes(
        x = 0, y = 0,
        xend = PC1, yend = PC2
      ),
      arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 25),
      linewidth = 0.3, color = "darkred"
    ) +
      xlab(paste0("PC",b, "(",PC1val,"%)")) +
      ylab(paste0("PC",c,"(",PC2val,"%)"))
    
    #add the boxes
    #Biplot <- Biplot + ggrepel::geom_label_repel(
    Biplot <- Biplot + ggrepel::geom_text_repel(
      data = loadings,
      #box.padding = 0.1, point.padding = 0.3, max.overlaps = 150,
      aes(
        label = rownames(loadings),
        x = PC1 * scale*loadingsSCALE,
        y = PC2 * scale*loadingsSCALE
      ),
      size = 3, # Change the font size of the text here
      color = "black", # Change the color of the text here
      #arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 15), # arrow does not look good
      force = 4)
    #Biplot
    
  }else{
  if (Color != FALSE & Shape != FALSE & Show_Loadings == TRUE & Scaling== TRUE){
    print("yes color yes shape yes loadings yes scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>% select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- tibble::column_to_rownames(data, "Row.names")
    
    data<- merge(data, mdata %>%select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- tibble::column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete
    
    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Shape, -Color)),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   fill = Color,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F,
                   loadings=T, #draws Eigenvectors
                   loadings.label = F,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select(-Shape, -Color)),scale. = TRUE)
  #  ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName,"_coT_shT_loT_scT",  ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #2
  }else if(Color != FALSE & Shape != FALSE & Show_Loadings == TRUE & Scaling== FALSE){
    print("yes color yes shape yes loadings no scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>% select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- tibble::column_to_rownames(data, "Row.names")
    
    data<- merge(data, mdata %>% select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- tibble::column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete
    
    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Shape, -Color)),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   fill = Color,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F,
                   loadings=T, #draws Eigenvectors
                   loadings.label = F,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select(-Shape, -Color)),scale. = FALSE)
   # ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coT_shT_loT_scF", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    
    #3 Done
  }else if(Color != FALSE & Shape != FALSE & Show_Loadings == FALSE & Scaling== TRUE){
    print("yes color yes shape no loadings yes scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>%select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- tibble::column_to_rownames(data, "Row.names")
    
    data<- merge(data, mdata %>%select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- tibble::column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete
    
    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Shape, -Color)),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   fill = Color,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F) +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select(-Shape, -Color)),scale. = TRUE)
  #  ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coT_shT_loF_scT", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #4 Done
  } else if(Color != FALSE & Shape != FALSE & Show_Loadings == FALSE & Scaling== FALSE){
    print("yes color yes shape no loadings no scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>%select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- tibble::column_to_rownames(data, "Row.names")
    
    data<- merge(data, mdata %>%select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- tibble::column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete
    
    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Shape, -Color)),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   fill = Color,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F) +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select(-Shape, -Color)),scale. = FALSE)
   # ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coT_shT_loF_scF", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #5 Done
  } else if (Color == FALSE & Shape != FALSE & Show_Loadings == TRUE & Scaling== TRUE){
    print("no color yes shape yes loadings yes scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>%select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- tibble::column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete
    
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Shape)),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F,
                   loadings=T, #draws Eigenvectors
                   loadings.label = F,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select(-Shape)),scale. = TRUE)
  #  ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coF_shT_loT_scT",".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    
    
    #6 Done
  }else if(Color == FALSE & Shape != FALSE & Show_Loadings == TRUE & Scaling== FALSE){
    print("no color yes shape yes loadings no scale")
    #data <-Processed_data
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    
    data<- merge(data, mdata %>%select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- tibble::column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete
    
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Shape)),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F,
                   loadings=T, #draws Eigenvectors
                   loadings.label = F,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select(-Shape)),scale. = FALSE)
   # ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coF_shT_loT_scF", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    
    #7 Done
  }else if(Color == FALSE & Shape != FALSE & Show_Loadings == FALSE & Scaling== TRUE){
    print("no color yes shape no loadings yes scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>%select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- tibble::column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete
    
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Shape)),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F) +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select(-Shape)),scale. = TRUE)
   # ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coF_shT_loF_scT", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #8 Done
  } else if(Color == FALSE & Shape != FALSE & Show_Loadings == FALSE & Scaling== FALSE){
    print("no color yes shape no loadings no scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>%select(Shape), by=0) # merge the data with the design to get only the kept samples
    data <- tibble::column_to_rownames(data, "Row.names")
    data[,Shape] <- as.factor(data[,Shape]) # make the shape into a factor to be discrete
    
    if(length(unique(data[,Shape]))>6){
      stop("Error. You tried to plot more than 6 shapes. It would be preferable to use color instead of shape" )
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Shape)),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   shape = Shape,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F) +
      scale_shape_manual(values=c(22,21,24,23,25,7,8,11,12))+ #needed if more than 6 shapes are in place
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select(-Shape)),scale. = FALSE)
   # ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coF_shT_loF_scF", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    
    #9 Done
  } else if (Color != FALSE & Shape == FALSE & Show_Loadings == TRUE & Scaling== TRUE){
    print("yes color no shape yes loadings yes scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>%select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- tibble::column_to_rownames(data, "Row.names")
    
    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Color)),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F,
                   loadings=T, #draws Eigenvectors
                   loadings.label = F,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select( -Color)),scale. = TRUE)
  #  ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName,"_coT_shF_loT_scT", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #10 Done
  }else if(Color != FALSE & Shape == FALSE & Show_Loadings == TRUE & Scaling== FALSE){
    print("yes color no shape yes loadings no scale")
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>%select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- tibble::column_to_rownames(data, "Row.names")
    
    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Color)),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F,
                   loadings=T, #draws Eigenvectors
                   loadings.label = F,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select( -Color)),scale. = FALSE)
  #  ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName,"_coT_shF_loT_scF", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #11 Done
  }else if(Color != FALSE & Shape == FALSE & Show_Loadings == FALSE & Scaling== TRUE){
    print("yes color no shape no loadings yes scale")
    
    mdata<- merge(as.data.frame(data), as.data.frame(Design), by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    
    data<- merge(data, mdata %>%select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- tibble::column_to_rownames(data, "Row.names")
    
    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Color)),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F) +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select( -Color)),scale. = TRUE)
  #  ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName,"_coT_shF_loF_scT", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #12 Done
  } else if(Color != FALSE & Shape == FALSE & Show_Loadings == FALSE & Scaling== FALSE){
    print("yes color no shape no loadings no scale")
    
    
    mdata<- merge(data, Design, by=0) # merge the data with the design to get only the kept samples
    mdata <- tibble::column_to_rownames(mdata, "Row.names")
    data<- merge(data, mdata %>%select(Color), by=0) # Add the selected columns(the ones wanted to plot) and add them to the dat
    data <- tibble::column_to_rownames(data, "Row.names")
    
    # For color if we have character we go with disctere. If we have numeriic we to discrete untill 4groupd. if e have more we go for continouus
    if(is.numeric(data[,Color])==TRUE | is.integer(data[,Color])==TRUE){
      if(length(unique(data[,Color]))>4){ # change this to change the number after which color becomes from distinct to continuous
        data[,Color] <- as.numeric(data[,Color])
      }else(data[,Color] <- as.factor(data[,Color]))
    }
    
    PCA<- autoplot(prcomp(as.matrix(data%>% select(-Color)),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   colour = Color,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F) +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data%>% select( -Color)),scale. = FALSE)
  #  ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName,"_coT_shF_loF_scF", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #13 Done
  } else if (Color == FALSE & Shape == FALSE & Show_Loadings == TRUE & Scaling== TRUE){
    print("no color no shape yes loadings yes scale")
    
    mdata<- data
    PCA<- autoplot(prcomp(as.matrix(data),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F,
                   loadings=T, #draws Eigenvectors
                   loadings.label = F,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data),scale. = TRUE)
   # ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coF_shF_loT_scT", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #14 Done
  }else if(Color == FALSE & Shape == FALSE & Show_Loadings == TRUE & Scaling== FALSE){
    print("no color no shape yes loadings no scale")
    
    mdata<- data
    PCA<- autoplot(prcomp(as.matrix(data),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F,
                   loadings=T, #draws Eigenvectors
                   loadings.label = F,
                   loadings.label.vjust = 1.2,
                   loadings.label.size=2.5,
                   loadings.colour="grey10",
                   loadings.label.colour="grey10") +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data),scale. = FALSE)
   # ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName,"_coF_shF_loT_scF", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #15 Done
  }else if(Color == FALSE &  Shape == FALSE & Show_Loadings == FALSE & Scaling== TRUE){
    print("no color no shape no loadings yes scale")
    
    mdata<- data
    PCA<- autoplot(prcomp(as.matrix(data),scale. = TRUE),   # Run and plot PCA
                   data= data,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F) +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data),scale. = TRUE)
   # ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName,"_coF_shF_loF_scT", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    #16 Done
  } else if(Color == FALSE &  Shape == FALSE & Show_Loadings == FALSE & Scaling== FALSE){
    print("no color no shape no loadings no scale")
    
    mdata<- data
    PCA<- autoplot(prcomp(as.matrix(data),scale. = FALSE),   # Run and plot PCA
                   data= data,
                   size = 3,
                   alpha = 0.8,
                   label=FALSE,
                   label.size=2.5,
                   label.repel = F) +
      ggtitle(paste(OutputPlotName)) +
      Theme +
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)
    
    loading_data <- prcomp(as.matrix(data),scale. = FALSE)
    #ggsave(file=paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coF_shF_loF_scF", ".",Save_as, sep=""), plot=PCA, width=10, height=10)
    
    
    # For the selected few
  }  else {
    print("What are you doing? You've got it all wrong!")
  }
  
  # svg(filename = paste("Results_", Sys.Date(), "/PCA_plots/", OutputPlotName, "_coF_shF_loF_scF", ".svg", sep=""),
  #     width = 10,
  #     height = 8)
  # plot(PCA)
  # dev.off()
  }
  loading_data_table <-as.data.frame(loading_data$rotation)
  loading_data_table <- loading_data_table[,1:2]
  loading_data_table <- tibble::rownames_to_column(loading_data_table, "Metabolite")
  writexl::write_xlsx(loading_data_table, paste(Results_folder_plots_PCA_folder,"/", OutputPlotName, "_Loadings.xlsx", sep=""), col_names = TRUE)
  
  if (Fill != FALSE){
  ggsave(file=paste(Results_folder_plots_PCA_folder,"/", OutputPlotName,"PC",b, "PC",c,".",Save_as, sep=""), plot=PCA, width=10, height=10)
  ggsave(file=paste(Results_folder_plots_PCA_folder,"/Biplot_", OutputPlotName, "PC",b, "PC",c,".",Save_as, sep=""), plot=Biplot, width=10, height=10)
  }else{
    ggsave(file=paste(Results_folder_plots_PCA_folder,"/", OutputPlotName,".",Save_as, sep=""), plot=PCA, width=10, height=10)
  }
}

##########-------------------------
### Use function
#plotPCA(data= preprocessing_output[["data_processed"]], Scaling =TRUE, Show_Loadings = FALSE, Design= preprocessing_output[["Experimental_design"]], Color = Conditions, Shape = time, OutputPlotName= 'PCAftw')
#plotPCA(data= preprocessing_output[["data_processed"]], Scaling =TRUE, Show_Loadings = FALSE)
##########--------------------------



#####################################
### ### ### Volcano Plots ### ### ###
#####################################
#' Notes.
#' careful with x and y limits.
#' save as pdf or png user decides
#' save in results_volcano plots, and then plots
#' also after normal volcano add volcano multiple that it plots 2 plots Together i.e. c1vsko1 and c2vsko2 in the same plots
#' Change the shapes to something normal

MetaProVizVolcano <-function(data =data,
                            pathway= FALSE,
                            test = "p.adj",
                            pCutoff= 0.05,
                            FCcutoff=0.5,
                            OutputPlotName= 'Volcano',
                            Theme=theme_classic(),
                            xlab =NULL, ylab=NULL, data2=NULL,
                            Cond1name="Comparisson 1", Cond2name="Comparisson 2",
                            connectors = FALSE,
                            Save_as = svg){
  
  
  
  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_Volcano_folder = paste(Results_folder,"/Volcano_Plots",  sep="")
  if (!dir.exists(Results_folder_plots_Volcano_folder)) {dir.create(Results_folder_plots_Volcano_folder)}
  
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  Save_as= deparse(substitute(Save_as))
  
  Multiple = FALSE
  if(is.null(data2)!="TRUE"){
    Multiple = TRUE
  }
  
  if(Multiple == FALSE & pathway==FALSE){
    if(test=="p.adj"){
      # Change plot labs if the user has put the input
      if(is.null(xlab)){
        xlab=bquote(~Log[2]~ FC)
      }
      if(is.null(ylab)){
        ylab=bquote(~-Log[10]~p.adj)
      }
      Plot<- EnhancedVolcano::EnhancedVolcano(data,
                              lab = data$Metabolite,#Metabolite name
                              x = "Log2FC",#Log2FC
                              y = "p.adj",#p-value or q-value
                              xlab  =xlab,
                              ylab =ylab,#(~-Log[10]~adjusted~italic(P))
                              pCutoff = pCutoff,
                              FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                              pointSize = 3,
                              labSize = 3,
                              titleLabSize = 16,
                              # colCustom = c("black", "grey", "grey", "red"),
                              colAlpha = 1,
                              title= paste(OutputPlotName),
                              subtitle = bquote(italic("Differential metabolomics analysis")),
                              caption = paste0("total = ", nrow(data), " Metabolites"),
                              # xlim = c((ceiling(Reduce(min,data$Log2FC))-0.2),(ceiling(Reduce(max,data$Log2FC))+0.2)),
                              
                              xlim =  c(min(data$Log2FC[is.finite(data$Log2FC )])-0.2,max(data$Log2FC[is.finite(data$Log2FC )])+0.2  ),
                              
                              ylim = c(0,(ceiling(-log10(Reduce(min,data$p.adj))))),
                              cutoffLineType = "dashed",
                              cutoffLineCol = "black",
                              cutoffLineWidth = 0.5,
                              legendLabels=c('No changes',paste(FCcutoff,"< |Log2FC|"),paste("p.adj <",pCutoff) , paste('p.adj<',pCutoff,' &',FCcutoff,"< |Log2FC|")),
                              legendPosition = 'right',
                              legendLabSize = 8,
                              legendIconSize =4,
                              gridlines.major = FALSE,
                              gridlines.minor = FALSE,
                              drawConnectors = connectors)
      Plot <- Plot+Theme
     # ggsave(file=paste(Results_folder_plots_Volcano_folder, "/", OutputPlotName, ".", Save_as, sep=""), plot=Plot, width=12, height=9)
      
    }else if(test=="p.val"){
      # Change plot labs if the user has put the input
      if(is.null(xlab)){
        xlab=bquote(~Log[2]~ FC)
      }
      if(is.null(ylab)){
        ylab=bquote(~-Log[10]~p.val)
      }
      Plot<- EnhancedVolcano::EnhancedVolcano (data,
                              lab = data$Metabolite,#Metabolite name
                              x = "Log2FC",#Log2FC
                              y = "p.val",#p-value or q-value
                              xlab = xlab,
                              ylab = ylab,#(~-Log[10]~adjusted~italic(P))
                              pCutoff = pCutoff,
                              FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                              pointSize = 3,
                              labSize = 3,
                              titleLabSize = 16,
                              #    colCustom = c("black", "grey", "grey", "red"),
                              colAlpha = 1,
                              title= paste(OutputPlotName),
                              subtitle = bquote(italic("Differential metabolomics analysis")),
                              caption = paste0("total = ", nrow(data), " Metabolites"),
                              #xlim = c((ceiling(Reduce(min,data$Log2FC))-0.2),(ceiling(Reduce(max,data$Log2FC))+0.2)),
                              xlim =  c(min(data$Log2FC[is.finite(data$Log2FC )])-0.2,max(data$Log2FC[is.finite(data$Log2FC )])+0.2  ),
                              ylim = c(0,(ceiling(-log10(Reduce(min,data$p.val))))),
                              cutoffLineType = "dashed",
                              cutoffLineCol = "black",
                              cutoffLineWidth = 0.5,
                              legendLabels=c('No changes',paste(FCcutoff,"< |Log2FC|"),paste("p.val <",pCutoff) , paste('p.val<',pCutoff,' &',FCcutoff,"< |Log2FC|")),
                              legendPosition = 'right',
                              legendLabSize = 8,
                              legendIconSize =4,
                              gridlines.major = FALSE,
                              gridlines.minor = FALSE,
                              drawConnectors = connectors)
      Plot <- Plot+Theme
     # ggsave(file=paste("Results_", Sys.Date(), "/Volcano_plots/", OutputPlotName,  ".", Save_as , sep=""), plot=Plot, width=12, height=9)
    
    }else{
      print("Please select a significance test p.val or p.adj")
    }
  }else if(Multiple == FALSE & pathway==TRUE){
    if(test=="p.adj"){
      # Change plot labs if the user has put the input
      if(is.null(xlab)){
        xlab=bquote(~Log[2]~ FC)
      }
      if(is.null(ylab)){
        ylab=bquote(~-Log[10]~p.adj)
      }
      DMA_PathwaysPlot_IEC <- data
      
      #Make a list of metabolites that we want to see on the plot:
      Labels <- subset(DMA_PathwaysPlot_IEC,DMA_PathwaysPlot_IEC$Pathway != "unknown")
      Labels <-Labels[,1]
      
      #Prepare new colour scheme:
      keyvals <- ifelse(
        DMA_PathwaysPlot_IEC$Pathway == "amino acid and conjugate", "blue",
        ifelse(DMA_PathwaysPlot_IEC$Pathway == "glycolysis and PPP", "red",
               ifelse(DMA_PathwaysPlot_IEC$Pathway == "long chain acylcarnitine", "gold4",
                      ifelse(DMA_PathwaysPlot_IEC$Pathway == "nucleotides", "seagreen",
                             ifelse(DMA_PathwaysPlot_IEC$Pathway == "purine metabolism", "darkorchid1",
                                    ifelse(DMA_PathwaysPlot_IEC$Pathway == "pyrimidine metabolism", "darkorchid4",
                                           ifelse(DMA_PathwaysPlot_IEC$Pathway == "redox homeostasis", "orange",
                                                  ifelse(DMA_PathwaysPlot_IEC$Pathway == "short chain acylcarnitine", "gold1",
                                                         ifelse(DMA_PathwaysPlot_IEC$Pathway == "medium chain carnitine", "gold3",
                                                                ifelse(DMA_PathwaysPlot_IEC$Pathway == "TCA cycle", "firebrick4",
                                                                       "gray"))))))))))
                                                                       keyvals[is.na(keyvals)] <- 'gray'
                                                                         names(keyvals)[keyvals == 'gray'] <- 'Other'
                                                                         names(keyvals)[keyvals == 'blue'] <- "Amino acid and Conjugates"
                                                                         names(keyvals)[keyvals == 'red'] <- "Glycolysis and PPP"
                                                                         names(keyvals)[keyvals == 'gold4'] <- "Long chain Acylcarnitines"
                                                                         names(keyvals)[keyvals == 'seagreen'] <- "Nucleotides"
                                                                         names(keyvals)[keyvals == 'darkorchid1'] <- "Purine metabolism"
                                                                         names(keyvals)[keyvals == 'darkorchid4'] <- "Pyrimidine metabolism"
                                                                         names(keyvals)[keyvals == 'orange'] <- "Redox homeostasis"
                                                                         names(keyvals)[keyvals == 'gold1'] <- "Short chain Acylcarnitines"
                                                                         names(keyvals)[keyvals == 'firebrick4'] <- "TCA cycle"
                                                                         names(keyvals)[keyvals == 'gold3'] <- "Medium chain Carnitines"
                                                                         
                                                                         #Plot
                                                                         Plot<- EnhancedVolcano::EnhancedVolcano (DMA_PathwaysPlot_IEC,
                                                                                                 lab = DMA_PathwaysPlot_IEC$Metabolite,#Metabolite name
                                                                                                 selectLab =Labels,
                                                                                                 x = "Log2FC",#Log2FC
                                                                                                 y = "p.adj",#p-value or q-value
                                                                                                 xlab = xlab,
                                                                                                 ylab = ylab,#(~-Log[10]~adjusted~italic(P))
                                                                                                 pCutoff = pCutoff,
                                                                                                 FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                                                                 pointSize = 3,
                                                                                                 labSize = 3,
                                                                                                 titleLabSize = 16,
                                                                                                 colCustom = keyvals,
                                                                                                 colAlpha = 1,
                                                                                                 title= paste(OutputPlotName),
                                                                                                 subtitle = bquote(italic("Differential metabolomics analysis")),
                                                                                                 caption = paste0("total = ", nrow(DMA_PathwaysPlot_IEC), " Metabolites"),
                                                                                                 # xlim = c((ceiling(Reduce(min,data$Log2FC))-0.2),(ceiling(Reduce(max,data$Log2FC))+0.2)),
                                                                                                 xlim =  c(min(DMA_PathwaysPlot_IEC$Log2FC[is.finite(DMA_PathwaysPlot_IEC$Log2FC )])-0.2,max(DMA_PathwaysPlot_IEC$Log2FC[is.finite(DMA_PathwaysPlot_IEC$Log2FC )])+0.2  ),
                                                                                                 ylim = c(0,(ceiling(-log10(Reduce(min,data$p.adj))))),
                                                                                                 #xlim = c(-5,10),
                                                                                                 #ylim = c(0,65),
                                                                                                 #drawConnectors = TRUE,
                                                                                                 #widthConnectors = 0.5,
                                                                                                 #colConnectors = "black",
                                                                                                 #arrowheads=FALSE,
                                                                                                 cutoffLineType = "dashed",
                                                                                                 cutoffLineCol = "black",
                                                                                                 cutoffLineWidth = 0.5,
                                                                                                 legendLabels=c('No changes',paste(FCcutoff,"< |Log2FC|"),paste("p.adj <",pCutoff) , paste('p.adj<',pCutoff,' &',FCcutoff,"< |Log2FC|")),
                                                                                                 legendPosition = 'right',
                                                                                                 legendLabSize = 8,
                                                                                                 legendIconSize =4,
                                                                                                 gridlines.major = FALSE,
                                                                                                 gridlines.minor = FALSE,
                                                                                                 drawConnectors = connectors)
                                                                         Plot <- Plot+Theme  + labs(color='Pathway') 
                                                                       #  ggsave(file=paste("Results_", Sys.Date(), "/Volcano_plots/", OutputPlotName,  ".", Save_as , sep=""), plot=Plot, width=12, height=9)
                                                                        
    }else if(test=="p.val"){
      # Change plot labs if the user has put the input
      if(is.null(xlab)){
        xlab=bquote(~Log[2]~ FC)
      }
      if(is.null(ylab)){
        ylab=bquote(~-Log[10]~p.val)
      }
      DMA_PathwaysPlot_IEC <- data
      
      #Make a list of metabolites that we want to see on the plot:
      Labels <- subset(DMA_PathwaysPlot_IEC,DMA_PathwaysPlot_IEC$Pathway != "unknown")
      Labels <-Labels[,1]
      
      #Prepare new colour scheme:
      keyvals <- ifelse(
        DMA_PathwaysPlot_IEC$Pathway == "amino acid and conjugate", "blue",
        ifelse(DMA_PathwaysPlot_IEC$Pathway == "glycolysis and PPP", "red",
               ifelse(DMA_PathwaysPlot_IEC$Pathway == "long chain acylcarnitine", "gold4",
                      ifelse(DMA_PathwaysPlot_IEC$Pathway == "nucleotides", "seagreen",
                             ifelse(DMA_PathwaysPlot_IEC$Pathway == "purine metabolism", "darkorchid1",
                                    ifelse(DMA_PathwaysPlot_IEC$Pathway == "pyrimidine metabolism", "darkorchid4",
                                           ifelse(DMA_PathwaysPlot_IEC$Pathway == "redox homeostasis", "orange",
                                                  ifelse(DMA_PathwaysPlot_IEC$Pathway == "short chain acylcarnitine", "gold1",
                                                         ifelse(DMA_PathwaysPlot_IEC$Pathway == "medium chain carnitine", "gold3",
                                                                ifelse(DMA_PathwaysPlot_IEC$Pathway == "TCA cycle", "firebrick4",
                                                                       "gray"))))))))))
                                                                       keyvals[is.na(keyvals)] <- 'gray'
                                                                         names(keyvals)[keyvals == 'gray'] <- 'Other'
                                                                         names(keyvals)[keyvals == 'blue'] <- "Amino acid and Conjugates"
                                                                         names(keyvals)[keyvals == 'red'] <- "Glycolysis and PPP"
                                                                         names(keyvals)[keyvals == 'gold4'] <- "Long chain Acylcarnitines"
                                                                         names(keyvals)[keyvals == 'seagreen'] <- "Nucleotides"
                                                                         names(keyvals)[keyvals == 'darkorchid1'] <- "Purine metabolism"
                                                                         names(keyvals)[keyvals == 'darkorchid4'] <- "Pyrimidine metabolism"
                                                                         names(keyvals)[keyvals == 'orange'] <- "Redox homeostasis"
                                                                         names(keyvals)[keyvals == 'gold1'] <- "Short chain Acylcarnitines"
                                                                         names(keyvals)[keyvals == 'firebrick4'] <- "TCA cycle"
                                                                         names(keyvals)[keyvals == 'gold3'] <- "Medium chain Carnitines"
                                                                         
                                                                         #Plot
                                                                         Plot<- EnhancedVolcano::EnhancedVolcano (DMA_PathwaysPlot_IEC,
                                                                                                 lab = DMA_PathwaysPlot_IEC$Metabolite,#Metabolite name
                                                                                                 selectLab =Labels,
                                                                                                 x = "Log2FC",#Log2FC
                                                                                                 y = "p.val",#p-value or q-value
                                                                                                 xlab = xlab,
                                                                                                 ylab = ylab,#(~-Log[10]~adjusted~italic(P))
                                                                                                 pCutoff = pCutoff,
                                                                                                 FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                                                                 pointSize = 3,
                                                                                                 labSize = 3,
                                                                                                 titleLabSize = 16,
                                                                                                 colCustom = keyvals,
                                                                                                 colAlpha = 1,
                                                                                                 title= paste(OutputPlotName),
                                                                                                 subtitle = bquote(italic("Differential metabolomics analysis")),
                                                                                                 caption = paste0("total = ", nrow(DMA_PathwaysPlot_IEC), " Metabolites"),
                                                                                                 #xlim = c((ceiling(Reduce(min,data$Log2FC))-0.2),(ceiling(Reduce(max,data$Log2FC))+0.2)),
                                                                                                 xlim =  c(min(DMA_PathwaysPlot_IEC$Log2FC[is.finite(DMA_PathwaysPlot_IEC$Log2FC )])-0.2,max(DMA_PathwaysPlot_IEC$Log2FC[is.finite(DMA_PathwaysPlot_IEC$Log2FC )])+0.2  ),
                                                                                                 ylim = c(0,(ceiling(-log10(Reduce(min,data$p.val))))),
                                                                                                 #xlim = c(-5,10),
                                                                                                 #ylim = c(0,65),
                                                                                                 #drawConnectors = TRUE,
                                                                                                 #widthConnectors = 0.5,
                                                                                                 #colConnectors = "black",
                                                                                                 #arrowheads=FALSE,
                                                                                                 cutoffLineType = "dashed",
                                                                                                 cutoffLineCol = "black",
                                                                                                 cutoffLineWidth = 0.5,
                                                                                                 legendLabels=c('No changes',paste(FCcutoff,"< |Log2FC|"),paste("p.val <",pCutoff) , paste('p.val<',pCutoff,' &',FCcutoff,"< |Log2FC|")),
                                                                                                 legendPosition = 'right',
                                                                                                 legendLabSize = 8,
                                                                                                 legendIconSize =4,
                                                                                                 gridlines.major = FALSE,
                                                                                                 gridlines.minor = FALSE,
                                                                                                 drawConnectors = connectors)
                                                                         Plot <- Plot+Theme
                                                                       #  ggsave(file=paste("Results_", Sys.Date(), "/Volcano_plots/", OutputPlotName,  ".", Save_as, sep=""), plot=Plot, width=12, height=9)
                                                                         
    }else{
      print("Please select a significance test p.val or p.adj")
    }
  }else if(Multiple == TRUE & pathway==FALSE){
    print("Multiple yes pathway no")
    InputData_1 <- data
    Condition_1<- Cond1name
    InputData_2<- data2
    Condition_2<-  Cond2name
    
    
    #1. Include a column naming the set Proteomics or RNAseq:
    InputData_1[,"comparison"]  <- as.character("InputData1")
    InputData_2[,"comparison"]  <- as.character("InputData2")
    #2. Add the colour:
    InputData_1[,"colour"]  <- as.character("red")
    InputData_2[,"colour"]  <- as.character("blue")
    #3. Combine the files
    Combined_DMA <- rbind(InputData_1,InputData_2)
    #4.Prepare new colour scheme
    keyvals <- ifelse(
      Combined_DMA$colour == "red", "red",
      ifelse(Combined_DMA$colour == "blue", "blue",
             "black"))
             keyvals[is.na(keyvals)] <- 'black'
               names(keyvals)[keyvals == 'red'] <- paste(Condition_1)
               names(keyvals)[keyvals == 'blue'] <- paste(Condition_2)
               names(keyvals)[keyvals == 'black'] <- 'X'
               if(test=="p.adj"){
                 # Change plot labs if the user has put the input
                 if(is.null(xlab)){
                   xlab=bquote(~Log[2]~ FC)
                 }
                 if(is.null(ylab)){
                   ylab=bquote(~-Log[10]~p.adj)
                 }
                 Plot <- EnhancedVolcano::EnhancedVolcano (Combined_DMA,
                                          lab = Combined_DMA$Metabolite,#Metabolite name
                                          x = "Log2FC",#Log2FC
                                          y = "p.adj",#p-value or q-value
                                          xlab = xlab,
                                          ylab = ylab,#(~-Log[10]~adjusted~italic(P))
                                          pCutoff = pCutoff,
                                          FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                          pointSize = 3,
                                          labSize = 3,
                                          colCustom = keyvals,
                                          titleLabSize = 16,
                                          col=c("black", "grey", "grey", "purple"),#if you want to change colors
                                          colAlpha = 1,
                                          title=paste(OutputPlotName),
                                          subtitle = bquote(italic("Differential Metabolomics Analysis (DMA)")),
                                          caption = paste0("total = ", (nrow(Combined_DMA)/2), " Metabolites"),
                                          xlim =  c(min(Combined_DMA$Log2FC[is.finite(Combined_DMA$Log2FC )])-0.2,max(Combined_DMA$Log2FC[is.finite(Combined_DMA$Log2FC )])+0.2  ),
                                          ylim = c(0,(ceiling(-log10(Reduce(min,Combined_DMA$p.adj))))),
                                          #drawConnectors = TRUE,
                                          #widthConnectors = 0.5,
                                          #colConnectors = "black",
                                          cutoffLineType = "dashed",
                                          cutoffLineCol = "black",
                                          cutoffLineWidth = 0.5,
                                          legendLabels=c('No changes',paste(FCcutoff,"< |Log2FC|"),paste("p.adj <",pCutoff) , paste('p.adj<',pCutoff,' &',FCcutoff,"< |Log2FC|")),
                                          legendPosition = 'right',
                                          legendLabSize = 12,
                                          legendIconSize = 5.0,
                                          gridlines.major = FALSE,
                                          gridlines.minor = FALSE,
                                          drawConnectors = connectors)
                 Plot <- Plot+Theme
                # ggsave(file=paste("Results_", Sys.Date(), "/Volcano_plots/", OutputPlotName,  ".", Save_as, sep=""), plot=Plot, width=8, height=6)
                
               }else if(test=="p.val"){
                 if(is.null(xlab)){
                   xlab=bquote(~Log[2]~ FC)
                 }
                 if(is.null(ylab)){
                   ylab=bquote(~-Log[10]~p.val)
                 }
                 Plot <- EnhancedVolcano::EnhancedVolcano (Combined_DMA,
                                          lab = Combined_DMA$Metabolite,#Metabolite name
                                          x = "Log2FC",#Log2FC
                                          y = "p.val",#p-value or q-value
                                          xlab = xlab,
                                          ylab = ylab,#(~-Log[10]~adjusted~italic(P))
                                          pCutoff = pCutoff,
                                          FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                          pointSize = 3,
                                          labSize = 3,
                                          colCustom = keyvals,
                                          titleLabSize = 16,
                                          col=c("black", "grey", "grey", "purple"),#if you want to change colors
                                          colAlpha = 1,
                                          title=paste(OutputPlotName),
                                          subtitle = bquote(italic("Differential Metabolomics Analysis (DMA)")),
                                          caption = paste0("total = ", (nrow(Combined_DMA)/2), " Metabolites"),
                                          xlim =  c(min(Combined_DMA$Log2FC[is.finite(Combined_DMA$Log2FC )])-0.2,max(Combined_DMA$Log2FC[is.finite(Combined_DMA$Log2FC )])+0.2  ),
                                          ylim = c(0,(ceiling(-log10(Reduce(min,Combined_DMA$p.val))))),
                                          #drawConnectors = TRUE,
                                          #widthConnectors = 0.5,
                                          #colConnectors = "black",
                                          cutoffLineType = "dashed",
                                          cutoffLineCol = "black",
                                          cutoffLineWidth = 0.5,
                                          legendLabels=c('No changes',paste(FCcutoff,"< |Log2FC|"),paste("p.val <",pCutoff) , paste('p.val<',pCutoff,' &',FCcutoff,"< |Log2FC|")),
                                          legendPosition = 'right',
                                          legendLabSize = 12,
                                          legendIconSize = 5.0,
                                          gridlines.major = FALSE,
                                          gridlines.minor = FALSE,
                                          drawConnectors = connectors)
                 Plot <- Plot+Theme
                # ggsave(file=paste("Results_", Sys.Date(), "/Volcano_plots/", OutputPlotName,  ".", Save_as, sep=""), plot=Plot, width=12, height=9)
                
               }else{
                 print("Please select a significance test p.val or p.adj")
               }
  }else if(Multiple == TRUE & pathway==TRUE){
    print("Multiple yes pathway yes")
    InputData_1 <- data
    Condition_1<- Cond1name
    InputData_2<- data2
    Condition_2<-  Cond2name
    
    
    #1. Include a column naming the set Proteomics or RNAseq:
    InputData_1[,"comparison"]  <- as.character("InputData1")
    InputData_2[,"comparison"]  <- as.character("InputData2")
    #2. Add the colour:
    InputData_1[,"shape"]  <- 17
    InputData_2[,"shape"]  <- 64
    #3. Combine the files
    Combined_DMA <- rbind(InputData_1,InputData_2)
    #4.Prepare new colour scheme
    keyvalsshape <- ifelse(                  # shape doesnt work
      Combined_DMA$shape == 17, 17,
      ifelse(Combined_DMA$shape == 64, 64,
             4))
    keyvalsshape[is.na(keyvalsshape)] <- 4
    names(keyvalsshape)[keyvalsshape == 17] <- paste(Condition_1)
    names(keyvalsshape)[keyvalsshape == 64] <- paste(Condition_2)
    names(keyvalsshape)[keyvalsshape == 4] <- 'X'
    if(test=="p.adj"){
      # Change plot labs if the user has put the input
      if(is.null(xlab)){
        xlab=bquote(~Log[2]~ FC)
      }
      if(is.null(ylab)){
        ylab=bquote(~-Log[10]~p.adj)
      }
      DMA_PathwaysPlot_IEC <- Combined_DMA
      
      #Make a list of metabolites that we want to see on the plot:
      Labels <- subset(DMA_PathwaysPlot_IEC,DMA_PathwaysPlot_IEC$Pathway != "unknown")
      Labels <-Labels[,1]
      
      #Prepare new colour scheme:
      keyvals <- ifelse(
        DMA_PathwaysPlot_IEC$Pathway == "amino acid and conjugate", "blue",
        ifelse(DMA_PathwaysPlot_IEC$Pathway == "glycolysis and PPP", "red",
               ifelse(DMA_PathwaysPlot_IEC$Pathway == "long chain acylcarnitine", "gold4",
                      ifelse(DMA_PathwaysPlot_IEC$Pathway == "nucleotides", "seagreen",
                             ifelse(DMA_PathwaysPlot_IEC$Pathway == "purine metabolism", "darkorchid1",
                                    ifelse(DMA_PathwaysPlot_IEC$Pathway == "pyrimidine metabolism", "darkorchid4",
                                           ifelse(DMA_PathwaysPlot_IEC$Pathway == "redox homeostasis", "orange",
                                                  ifelse(DMA_PathwaysPlot_IEC$Pathway == "short chain acylcarnitine", "gold1",
                                                         ifelse(DMA_PathwaysPlot_IEC$Pathway == "medium chain carnitine", "gold3",
                                                                ifelse(DMA_PathwaysPlot_IEC$Pathway == "TCA cycle", "firebrick4",
                                                                       "gray"))))))))))
                                                                       keyvals[is.na(keyvals)] <- 'gray'
                                                                         names(keyvals)[keyvals == 'gray'] <- 'Other'
                                                                         names(keyvals)[keyvals == 'blue'] <- "Amino acid and Conjugates"
                                                                         names(keyvals)[keyvals == 'red'] <- "Glycolysis and PPP"
                                                                         names(keyvals)[keyvals == 'gold4'] <- "Long chain Acylcarnitines"
                                                                         names(keyvals)[keyvals == 'seagreen'] <- "Nucleotides"
                                                                         names(keyvals)[keyvals == 'darkorchid1'] <- "Purine metabolism"
                                                                         names(keyvals)[keyvals == 'darkorchid4'] <- "Pyrimidine metabolism"
                                                                         names(keyvals)[keyvals == 'orange'] <- "Redox homeostasis"
                                                                         names(keyvals)[keyvals == 'gold1'] <- "Short chain Acylcarnitines"
                                                                         names(keyvals)[keyvals == 'firebrick4'] <- "TCA cycle"
                                                                         names(keyvals)[keyvals == 'gold3'] <- "Medium chain Carnitines"
                                                                         
                                                                         #Plot
                                                                         #DMA_PathwaysPlot_IEC$shape <- as.numeric(DMA_PathwaysPlot_IEC$shape) ### shape doesnt work
                                                                         
                                                                         Plot<- EnhancedVolcano::EnhancedVolcano (DMA_PathwaysPlot_IEC,
                                                                                                 lab = DMA_PathwaysPlot_IEC$Metabolite,#Metabolite name
                                                                                                 selectLab =Labels,
                                                                                                 x = "Log2FC",#Log2FC
                                                                                                 y = "p.adj",#p-value or q-value
                                                                                                 xlab = xlab,
                                                                                                 ylab = ylab,#(~-Log[10]~adjusted~italic(P))
                                                                                                 pCutoff = pCutoff,
                                                                                                 FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                                                                 pointSize = 3,
                                                                                                 labSize = 3,
                                                                                                 titleLabSize = 16,
                                                                                                 colCustom = keyvals,
                                                                                                 shapeCustom = keyvalsshape,#### shape doesnt work
                                                                                                 colAlpha = 1,
                                                                                                 title= paste(OutputPlotName),
                                                                                                 subtitle = bquote(italic("Differential metabolomics analysis")),
                                                                                                 caption = paste0("total = ", nrow(DMA_PathwaysPlot_IEC), " Metabolites"),
                                                                                                 xlim =  c(min(DMA_PathwaysPlot_IEC$Log2FC[is.finite(DMA_PathwaysPlot_IEC$Log2FC )])-0.2,max(DMA_PathwaysPlot_IEC$Log2FC[is.finite(DMA_PathwaysPlot_IEC$Log2FC )])+0.2  ),
                                                                                                 ylim = c(0,(ceiling(-log10(Reduce(min,data$p.adj))))),
                                                                                                 #xlim = c(-5,10),
                                                                                                 #ylim = c(0,65),
                                                                                                 #drawConnectors = TRUE,
                                                                                                 #widthConnectors = 0.5,
                                                                                                 #colConnectors = "black",
                                                                                                 #arrowheads=FALSE,
                                                                                                 cutoffLineType = "dashed",
                                                                                                 cutoffLineCol = "black",
                                                                                                 cutoffLineWidth = 0.5,
                                                                                                 legendLabels=c('No changes',paste(FCcutoff,"< |Log2FC|"),paste("p.adj <",pCutoff) , paste('p.adj<',pCutoff,' &',FCcutoff,"< |Log2FC|")),
                                                                                                 legendPosition = 'right',
                                                                                                 legendLabSize = 8,
                                                                                                 legendIconSize =4,
                                                                                                 gridlines.major = FALSE,
                                                                                                 gridlines.minor = FALSE,
                                                                                                 drawConnectors = connectors)
                                                                         Plot <- Plot+Theme
                                                                      #   ggsave(file=paste("Results_", Sys.Date(), "/Volcano_plots/", OutputPlotName,  ".", Save_as, sep=""), plot=Plot, width=12, height=9)
                                                                       
    }else if(test=="p.val"){
      # Change plot labs if the user has put the input
      if(is.null(xlab)){
        xlab=bquote(~Log[2]~ FC)
      }
      if(is.null(ylab)){
        ylab=bquote(~-Log[10]~p.val)
      }
      DMA_PathwaysPlot_IEC <- Combined_DMA
      
      #Make a list of metabolites that we want to see on the plot:
      Labels <- subset(DMA_PathwaysPlot_IEC,DMA_PathwaysPlot_IEC$Pathway != "unknown")
      Labels <-Labels[,1]
      
      #Prepare new colour scheme:
      keyvals <- ifelse(
        DMA_PathwaysPlot_IEC$Pathway == "amino acid and conjugate", "blue",
        ifelse(DMA_PathwaysPlot_IEC$Pathway == "glycolysis and PPP", "red",
               ifelse(DMA_PathwaysPlot_IEC$Pathway == "long chain acylcarnitine", "gold4",
                      ifelse(DMA_PathwaysPlot_IEC$Pathway == "nucleotides", "seagreen",
                             ifelse(DMA_PathwaysPlot_IEC$Pathway == "purine metabolism", "darkorchid1",
                                    ifelse(DMA_PathwaysPlot_IEC$Pathway == "pyrimidine metabolism", "darkorchid4",
                                           ifelse(DMA_PathwaysPlot_IEC$Pathway == "redox homeostasis", "orange",
                                                  ifelse(DMA_PathwaysPlot_IEC$Pathway == "short chain acylcarnitine", "gold1",
                                                         ifelse(DMA_PathwaysPlot_IEC$Pathway == "medium chain carnitine", "gold3",
                                                                ifelse(DMA_PathwaysPlot_IEC$Pathway == "TCA cycle", "firebrick4",
                                                                       "gray"))))))))))
                                                                       keyvals[is.na(keyvals)] <- 'gray'
                                                                         names(keyvals)[keyvals == 'gray'] <- 'Other'
                                                                         names(keyvals)[keyvals == 'blue'] <- "Amino acid and Conjugates"
                                                                         names(keyvals)[keyvals == 'red'] <- "Glycolysis and PPP"
                                                                         names(keyvals)[keyvals == 'gold4'] <- "Long chain Acylcarnitines"
                                                                         names(keyvals)[keyvals == 'seagreen'] <- "Nucleotides"
                                                                         names(keyvals)[keyvals == 'darkorchid1'] <- "Purine metabolism"
                                                                         names(keyvals)[keyvals == 'darkorchid4'] <- "Pyrimidine metabolism"
                                                                         names(keyvals)[keyvals == 'orange'] <- "Redox homeostasis"
                                                                         names(keyvals)[keyvals == 'gold1'] <- "Short chain Acylcarnitines"
                                                                         names(keyvals)[keyvals == 'firebrick4'] <- "TCA cycle"
                                                                         names(keyvals)[keyvals == 'gold3'] <- "Medium chain Carnitines"
                                                                         
                                                                         #Plot
                                                                         #DMA_PathwaysPlot_IEC$shape <- as.numeric(DMA_PathwaysPlot_IEC$shape) ### shape doesnt work
                                                                         
                                                                         Plot<- EnhancedVolcano::EnhancedVolcano (DMA_PathwaysPlot_IEC,
                                                                                                 lab = DMA_PathwaysPlot_IEC$Metabolite,#Metabolite name
                                                                                                 selectLab =Labels,
                                                                                                 x = "Log2FC",#Log2FC
                                                                                                 y = "p.val",#p-value or q-value
                                                                                                 xlab = xlab,
                                                                                                 ylab = ylab,#(~-Log[10]~adjusted~italic(P))
                                                                                                 pCutoff = pCutoff,
                                                                                                 FCcutoff = FCcutoff,#Cut off Log2FC, automatically 2
                                                                                                 pointSize = 3,
                                                                                                 labSize = 3,
                                                                                                 titleLabSize = 16,
                                                                                                 colCustom = keyvals,
                                                                                                 shapeCustom = keyvalsshape,
                                                                                                 colAlpha = 1,
                                                                                                 title= paste(OutputPlotName),
                                                                                                 subtitle = bquote(italic("Differential metabolomics analysis")),
                                                                                                 caption = paste0("total = ", nrow(DMA_PathwaysPlot_IEC), " Metabolites"),
                                                                                                 xlim =  c(min(DMA_PathwaysPlot_IEC$Log2FC[is.finite(DMA_PathwaysPlot_IEC$Log2FC )])-0.2,max(DMA_PathwaysPlot_IEC$Log2FC[is.finite(DMA_PathwaysPlot_IEC$Log2FC )])+0.2  ),
                                                                                                 ylim = c(0,(ceiling(-log10(Reduce(min,data$p.val))))),
                                                                                                 #xlim = c(-5,10),
                                                                                                 #ylim = c(0,65),
                                                                                                 #drawConnectors = TRUE,
                                                                                                 #widthConnectors = 0.5,
                                                                                                 #colConnectors = "black",
                                                                                                 #arrowheads=FALSE,
                                                                                                 cutoffLineType = "dashed",
                                                                                                 cutoffLineCol = "black",
                                                                                                 cutoffLineWidth = 0.5,
                                                                                                 legendLabels=c('No changes',paste(FCcutoff,"< |Log2FC|"),paste("p.val <",pCutoff) , paste('p.val<',pCutoff,' &',FCcutoff,"< |Log2FC|")),
                                                                                                 legendPosition = 'right',
                                                                                                 legendLabSize = 8,
                                                                                                 legendIconSize =4,
                                                                                                 gridlines.major = FALSE,
                                                                                                 gridlines.minor = FALSE,
                                                                                                 drawConnectors = connectors)
                                                                         Plot <- Plot+Theme
                                                                        # ggsave(file=paste("Results_", Sys.Date(), "/Volcano_plots/", OutputPlotName,  ".", Save_as, sep=""), plot=Plot, width=12, height=9)
                                                                        
                                                                         
    }else{
      print("Please select a significance test p.val or p.adj")
    }
  }
  
  # # save as svg
  # svg(filename = paste("Results_", Sys.Date(), "/Volcano_plots/", OutputPlotName, ".svg", sep=""),
  #     width = 10,
  #     height = 8)
  # plot(Plot)
  # dev.off()
  ggsave(file=paste(Results_folder_plots_Volcano_folder, "/", OutputPlotName,  ".", Save_as , sep=""), plot=Plot, width=12, height=9)
  
}

###------------------------
### Test the Volcano plots
## multi false path false
#plotVolcano(data =DMA_output ,pathway= FALSE, test = "p.adj", pCutoff= 0.05, FCcutoff=1, OutputPlotName= 'Volcano yes yes padj')
## muli false path true
#plotVolcano(data =DMA_output ,pathway= TRUE, test = "p.adj", pCutoff= 0.05, FCcutoff=0.9, OutputPlotName= 'Volcano yes path true')
## test multiple true pathways false
#DMA_output2 <- DMA_output
#DMA_output2$Log2FC<- DMA_output2$Log2FC+1
#plotVolcano(data =DMA_output ,pathway= FALSE, test = "p.adj", pCutoff= 0.05, FCcutoff=1, OutputPlotName= 'Volc Multi yes paths no', data2=DMA_output2)
## test multiple true pathways true
#plotVolcano(data =DMA_output ,pathway= TRUE, test = "p.adj", pCutoff= 0.05, FCcutoff=1, OutputPlotName= 'VOlc all in55', data2=DMA_output2)
###------------------------


######################################
### ### ### Alluvial Plots ### ### ###
######################################=
#' Notes; Check the alluvial in the Metabolic Clusters. This one is a bit old but kept here to return to if nedded

MetaProVizMCAlluvial <- function(Input_data1, Input_data2,  Output_Name = "Metabolic_Clusters_Output_Condition1-versus-Condition2",  Condition1,  Condition2,
                                 pCutoff= 0.05 , FCcutoff=0.5, test = "p.adj", plot_column_names= c("class", "MetaboliteChange_Significant", "Overall_Change", "Metabolite"),
                                 safe_colorblind_palette = c("#88CCEE",  "#DDCC77","#661100",  "#332288", "#AA4499","#999933",  "#44AA99", "#882255",  "#6699CC", "#117733", "#888888","#CC6677", "#FFF", "#000"), # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
                                 plot_color_variable = "Overall_Change",  plot_color_remove_variable = "SameDirection_NoChange",
                                 Save_as = pdf ){
  
  
  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_MetabolicCluster_folder = paste(Results_folder,"/Plots_AlluvialCluster",  sep="")
  if (!dir.exists(Results_folder_plots_MetabolicCluster_folder)) {dir.create(Results_folder_plots_MetabolicCluster_folder)}
  
  
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  Save_as_var <- Save_as
  
  C1 <- Input_data1
  C1 <- na.omit(C1)
  C1$class <- paste (Condition1)
  C2 <- Input_data2
  C2 <- na.omit(C2)
  C2$class <- paste (Condition2)
  
  # 0. Overall Regulation: label the selection of metabolites that change or where at least we have a change in one of the two conditions
  if(test== "p.val"){
    C1 <- C1 %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
    C2 <- C2%>%
      mutate(MetaboliteChange_Significant1 = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                       Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                       TRUE ~ 'No_Change'))%>%
      rename("class1"="class",
             "Log2FC1"="Log2FC",
             "p.val1"="p.val",
             "p.adj1"="p.adj")
    
  }else{
    C1 <- C1 %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
    C2 <- C2%>%
      mutate(MetaboliteChange_Significant1 = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                       Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                       TRUE ~ 'No_Change'))%>%
      rename("class1"="class",
             "Log2FC1"="Log2FC",
             "p.val1"="p.val",
             "p.adj1"="p.adj")
  }
  
  
  
  MergeDF<- merge(C1, C2[,c("Metabolite","MetaboliteChange_Significant1", "class1","Log2FC1","p.val1","p.adj1")], by="Metabolite")%>%
    mutate(Change_Specific = case_when((class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ 'OppositeChange',
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ 'OppositeChange',
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ 'SameDirection_UP',
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ 'SameDirection_DOWN',
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "No_Change")& (class1== paste(Condition2)& MetaboliteChange_Significant1 == "DOWN") ~ paste("ChangeOnly", Condition2, "DOWN", sep="_"),
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "No_Change")& (class1== paste(Condition2)& MetaboliteChange_Significant1 == "UP") ~ paste("ChangeOnly", Condition2, "UP", sep="_"),
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "UP")       & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "No_Change") ~ paste("ChangeOnly", Condition1, "UP", sep="_"),
                                       (class== paste(Condition1)& MetaboliteChange_Significant == "DOWN")     & (class1== paste(Condition2)& MetaboliteChange_Significant1 == "No_Change") ~ paste("ChangeOnly", Condition1, "DOWN", sep="_"),
                                       TRUE ~ 'SameDirection_NoChange'))
  
  
  MergeDF_C1 <-MergeDF %>% select(-c( "MetaboliteChange_Significant1", "class1", "Log2FC1", "p.val1", "p.adj1")) %>%
    unite(col=UniqueID, c(Metabolite, MetaboliteChange_Significant, class), sep = "_", remove = FALSE, na.rm = FALSE)%>%
    mutate(Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                       Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP',
                                       Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_DOWN',
                                       Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                       Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~paste("ChangeOnly", Condition2, "DOWN", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") ~paste("ChangeOnly", Condition2, "UP", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~paste("Unique", Condition1,"DOWN", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~paste("Unique", Condition1, "UP", sep="_"),
                                       TRUE ~ paste("Unique", Condition1, sep="_")))
  
  MergeDF_C2 <- MergeDF %>% select(-c( "MetaboliteChange_Significant", "class", "Log2FC", "p.val", "p.adj")) %>%
    unite(col=UniqueID, c(Metabolite, MetaboliteChange_Significant1, class1), sep = "_", remove = FALSE, na.rm = FALSE)%>%
    rename("class"="class1",
           "MetaboliteChange_Significant"="MetaboliteChange_Significant1",
           "Log2FC"="Log2FC1",
           "p.val"="p.val1",
           "p.adj"="p.adj1")%>%
    mutate(Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                       Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP',
                                       Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_DOWN',
                                       Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                       Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~paste("ChangeOnly", Condition1, "DOWN", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~paste("ChangeOnly", Condition1, "UP", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~paste("Unique", Condition2,"DOWN", sep="_"),
                                       Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") ~paste("Unique", Condition2, "UP", sep="_"),
                                       TRUE ~ paste("Unique", Condition1, sep="_")))
  
  
  Alluvial_DF <- rbind(MergeDF_C1, MergeDF_C2)
  Alluvial_DF<- Alluvial_DF%>%
    mutate(Amount_Change_Specific = case_when(Change_Specific== 'OppositeChange' ~ paste((sum(Alluvial_DF$Change_Specific=="OppositeChange", na.rm=T))/2),
                                              Change_Specific== 'SameDirection_UP' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_UP", na.rm=T))/2),
                                              Change_Specific== 'SameDirection_DOWN' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_DOWN", na.rm=T))/2),
                                              Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition1, "UP", sep="_"), na.rm=T)),
                                              Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition1, "DOWN", sep="_"), na.rm=T)),
                                              Change_Specific== paste("ChangeOnly", Condition2, "UP" ,sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition2, "UP", sep="_"), na.rm=T)),
                                              Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") ~ paste(sum(Alluvial_DF$Change_Specific==paste("ChangeOnly", Condition2, "DOWN", sep="_"), na.rm=T)),
                                              Change_Specific== 'SameDirection_NoChange' ~ paste((sum(Alluvial_DF$Change_Specific=="SameDirection_NoChange", na.rm=T))/2),
                                              Change_Specific== paste("Unique", Condition1,"DOWN", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition1,"DOWN", sep="_"), na.rm=T)),
                                              Change_Specific== paste("Unique", Condition1,"UP", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition1,"UP", sep="_"), na.rm=T)),
                                              Change_Specific== paste("Unique", Condition2,"DOWN", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition2,"DOWN", sep="_"), na.rm=T)),
                                              Change_Specific== paste("Unique", Condition2,"UP", sep="_") ~paste(sum(Alluvial_DF$Change_Specific==paste("Unique", Condition2,"UP", sep="_"), na.rm=T)),
                                              TRUE ~ 'FALSE'))
  Alluvial_DF<- Alluvial_DF%>%
    mutate(Overall_Change = case_when(Change_Specific== 'OppositeChange' ~ 'OppositeChange',
                                      Change_Specific== 'SameDirection_UP'~ 'SameDirection_UP_or_DOWN',
                                      Change_Specific== 'SameDirection_DOWN' ~ 'SameDirection_UP_or_DOWN',
                                      Change_Specific== 'SameDirection_NoChange' ~ 'SameDirection_NoChange',
                                      Change_Specific== paste("ChangeOnly", Condition2, "DOWN", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition2, sep="_"),
                                      Change_Specific== paste("ChangeOnly", Condition2, "UP", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition2, sep="_"),
                                      Change_Specific== paste("Unique", Condition1, "DOWN",sep="_") ~paste("Unique", Condition1, sep="_"),
                                      Change_Specific== paste("Unique", Condition1, "UP",sep="_") ~paste("Unique", Condition1, sep="_"),
                                      Change_Specific== paste("ChangeOnly", Condition1, "DOWN", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition1, sep="_"),
                                      Change_Specific== paste("ChangeOnly", Condition1, "UP", sep="_") & MetaboliteChange_Significant == "No_Change" ~ paste("ChangeOnly", Condition1, sep="_"),
                                      Change_Specific== paste("Unique", Condition2,"DOWN", sep="_") ~paste("Unique", Condition2, sep="_"),
                                      Change_Specific== paste("Unique", Condition2,"UP", sep="_") ~paste("Unique", Condition2, sep="_"),
                                      TRUE ~ "FALSE"))
  Alluvial_DF <- Alluvial_DF %>%
    mutate(Amount_Overall_Change = case_when(Overall_Change== 'OppositeChange' ~ paste((sum(Alluvial_DF$Overall_Change=="OppositeChange", na.rm=T))/2),
                                             Overall_Change== 'SameDirection_UP_or_DOWN' ~ paste((sum(Alluvial_DF$Overall_Change=="SameDirection_UP_or_DOWN", na.rm=T))/2),
                                             Overall_Change== paste("Unique", Condition1, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("Unique", Condition1, sep="_"), na.rm=T)),
                                             Overall_Change== paste("Unique", Condition2, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("Unique", Condition2, sep="_"), na.rm=T)),
                                             Overall_Change== paste("ChangeOnly", Condition1, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("ChangeOnly", Condition1, sep="_"), na.rm=T)),
                                             Overall_Change== paste("ChangeOnly", Condition2, sep="_") ~ paste(sum(Alluvial_DF$Overall_Change==paste("ChangeOnly", Condition2, sep="_"), na.rm=T)),
                                             Overall_Change== 'SameDirection_NoChange' ~ paste((sum(Alluvial_DF$Overall_Change=="SameDirection_NoChange", na.rm=T))/2),
                                             TRUE ~ 'FALSE'))
  
  # Create Alluvial final output df
  C1.final<- Alluvial_DF %>% filter(class == Condition1)
  C1.final <- C1.final %>% select(-c("UniqueID", "class", "Change_Specific","Amount_Change_Specific", "Overall_Change", "Amount_Overall_Change" ))
  C2.final<- Alluvial_DF %>% filter(class == Condition2)
  C2.final <- C2.final %>% select(-c("UniqueID", "class" ))
  Alluvial_DF.final <- merge(C1.final, C2.final, by = "Metabolite")
  names(Alluvial_DF.final) <- gsub(".x",paste(".",substr(Condition1, 1, 3), sep = ""),names(Alluvial_DF.final))
  names(Alluvial_DF.final) <- gsub(".y",paste(".",substr(Condition2, 1, 3), sep = ""),names(Alluvial_DF.final))
  names(Alluvial_DF.final) <- gsub(x = names(Alluvial_DF.final), pattern = "MetaboliteChange_Significant", replacement =  paste("MetaboliteChange_Significant_",test,pCutoff,"logFC",FCcutoff, sep = ""))
  
  ##Write to file
  # This is not needed fot the plots
  #  writexl::write_xlsx(Alluvial_DF.final, paste("Results_", Sys.Date(), "/MetabolicCluster_plots/","Metabolic_Clusters_Output_",Condition1,"-versus-",Condition2,Output_Name,".xlsx", sep = ""))
  # write.csv(Alluvial_DF2, paste("AlluvianDF", Output, ".csv", sep="_"), row.names= TRUE)
  
  
  # 1. Regulation:
  Alluvial_DF2  <- Alluvial_DF  %>%
    mutate(MetaboliteChange = case_when(Log2FC >= FCcutoff  ~ 'UP',
                                        Log2FC <= -FCcutoff ~ 'DOWN',
                                        TRUE ~ 'No_Change'))
  
  if (test == "p.val"){
    # 2. Excluded according to p-value:
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(Excluded_by_pval = case_when(p.val <= pCutoff ~ 'NO',
                                          p.val > pCutoff ~ 'YES'))
    #3. Excluded according to P-value?
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Excluded = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                   Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                   Log2FC >= FCcutoff & p.val >= pCutoff ~ 'UP_Excluded',
                                                   Log2FC <= -FCcutoff & p.val >= pCutoff ~ 'DOWN_Excluded',
                                                   TRUE ~ 'No_Change'))
    #4. After exlusion by p-value
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.val < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.val < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
  }else if(test=="p.adj"){
    # 2. Excluded according to p-value:
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(Excluded_by_pval = case_when(p.adj <= pCutoff ~ 'NO',
                                          p.adj > pCutoff ~ 'YES'))
    #3. Excluded according to P-value?
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Excluded = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                   Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                   Log2FC >= FCcutoff & p.adj >= pCutoff ~ 'UP_Excluded',
                                                   Log2FC <= -FCcutoff & p.adj >= pCutoff ~ 'DOWN_Excluded',
                                                   TRUE ~ 'No_Change'))
    #4. After exlusion by p-value
    Alluvial_DF2  <- Alluvial_DF2  %>%
      mutate(MetaboliteChange_Significant = case_when(Log2FC >= FCcutoff & p.adj < pCutoff ~ 'UP',
                                                      Log2FC <= -FCcutoff & p.adj < pCutoff ~ 'DOWN',
                                                      TRUE ~ 'No_Change'))
  }
  
  
  #5. Frequency:
  Alluvial_DF2[,"Frequency"]  <- as.numeric("1")
  #6. Safe DF
  # Alluvial_DF3 <- Alluvial_DF2[, c(1:6,12:14,7:8,10,9,11,15)]
  # Alluvial_Plot <- Alluvial_DF3[,c(6:14,2,15)]
  
  Alluvial_Plot <- Alluvial_DF2
  
  #Make the SelectionPlot:
  if(plot_color_remove_variable %in% Alluvial_Plot[,plot_color_variable] ){
    Alluvial_Plot<- Alluvial_Plot[-which(Alluvial_Plot[,plot_color_variable]==plot_color_remove_variable),]#remove the metabolites that do not change in either of the conditions
  }
  
  
  Save_as_var(paste(Results_folder_plots_MetabolicCluster_folder,"/Metabolic_Clusters_",Condition1,"-versus-",Condition2,"_", Output_Name,  ".",Save_as, sep=""), width=12, height=9)
  par(oma=c(2,2,8,2), mar = c(2, 2, 0.1, 2)+0.1)#https://www.r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
  alluvial( Alluvial_Plot %>% select(all_of(plot_column_names)), freq=Alluvial_Plot$Frequency,
            col = case_when(Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[1] ~ safe_colorblind_palette[1],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[2] ~ safe_colorblind_palette[2],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[3] ~ safe_colorblind_palette[3],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[4] ~ safe_colorblind_palette[4],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[5] ~ safe_colorblind_palette[5],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[6] ~ safe_colorblind_palette[6],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[7] ~ safe_colorblind_palette[7],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[8] ~ safe_colorblind_palette[8],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[9] ~ safe_colorblind_palette[9],
                            Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[10] ~ safe_colorblind_palette[10],
                            TRUE ~ 'black'),
            border = case_when(Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[1] ~ safe_colorblind_palette[1],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[2] ~ safe_colorblind_palette[2],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[3] ~ safe_colorblind_palette[3],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[4] ~ safe_colorblind_palette[4],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[5] ~ safe_colorblind_palette[5],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[6] ~ safe_colorblind_palette[6],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[7] ~ safe_colorblind_palette[7],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[8] ~ safe_colorblind_palette[8],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[9] ~ safe_colorblind_palette[9],
                               Alluvial_Plot[,plot_color_variable] == unique( Alluvial_Plot[,plot_color_variable])[10] ~ safe_colorblind_palette[10],
                               
                               TRUE ~ 'black'),
            hide = Alluvial_Plot$Frequency == 0,
            cex = 0.3,
            cex.axis=0.5)
  mtext("Selection of metabolites that change in at least one of the two conditions", side=3, line=6, cex=1.2, col="black", outer=TRUE) #https://www.r-graph-gallery.com/74-margin-and-oma-cheatsheet.html
  mtext(paste("",Output_Name), side=3, line=5, cex=0.8, col="black", outer=TRUE)
  mtext(paste("Underlying comparison: ",Condition1,"-versus-",Condition2), side=2, line=0, cex=0.8, col="black", outer=TRUE)
  #mtext("Legend", side=3, line=5, adj=1.0, cex=1, col="black", outer=TRUE)
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[1])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[1]), side=3, line=6, adj=0, cex=0.6, col=safe_colorblind_palette[1], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[2])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[2]), side=3, line=5, adj=0, cex=0.6, col=safe_colorblind_palette[2], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[3])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[3]), side=3, line=4, adj=0, cex=0.6, col=safe_colorblind_palette[3], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[4])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[4]), side=3, line=3, adj=0, cex=0.6, col=safe_colorblind_palette[4], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[5])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[5]), side=3, line=2, adj=0, cex=0.6, col=safe_colorblind_palette[5], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[6])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[6]), side=3, line=1, adj=0, cex=0.6, col=safe_colorblind_palette[6], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[7])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[7]), side=3, line=0, adj=0, cex=0.6, col=safe_colorblind_palette[7], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[8])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[8]), side=3, line=7, adj=1, cex=0.6, col=safe_colorblind_palette[8], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[9])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[9]), side=3, line=6, adj=1, cex=0.6, col=safe_colorblind_palette[9], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[10])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[10]), side=3, line=5, adj=1, cex=0.6, col=safe_colorblind_palette[10], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[11])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[11]), side=3, line=4, adj=1, cex=0.6, col=safe_colorblind_palette[11], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[12])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[12]), side=3, line=3, adj=1, cex=0.6, col=safe_colorblind_palette[12], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[13])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[13]), side=3, line=2, adj=1, cex=0.6, col=safe_colorblind_palette[13], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[14])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[14]), side=3, line=1, adj=1, cex=0.6, col=safe_colorblind_palette[14], outer=TRUE)
  }
  if( is.na(unique(Alluvial_Plot[,plot_color_variable])[15])==FALSE)  {
    mtext(paste(unique( Alluvial_Plot[,plot_color_variable])[15]), side=3, line=7, adj=1, cex=0.6, col=safe_colorblind_palette[15], outer=TRUE)
  }
  dev.off()# Close the pdf file
}



##########--------------------------
## use  function
#Alluvial_Plot <- plotAlluvial(Input1=dataNor, Input2=dataHyp, Condition1="Normoxia", Condition2="Hypoxia",test = "p.val", OutputPlotName= "Normaxia vs Hypoxia", Comparison="KO versus WT (DMEM)")
#MetaProVizplotMetabolicCluster(Input_data1 = DMA_output,  Input_data2 = DMA_output2, Condition1 = "Rot vs ctlr",
#                              Condition2 = "3NPA vs ctrl",    pCutoff = 0.05 ,FCcutoff = 0.5, test = "p.val", Output_Name = "lala",
#                              plot_column_names = c("class", "MetaboliteChange_Significant","Pathway","Metabolite"),
#                              plot_color_variable = "Pathway", plot_color_remove_variable = "unknown")

#######----------------------------

#####################################
### ### ### Lolipop Plots ### ### ###
#####################################
#' Should probaly add to be able to change the size and color variable to whatever we want
#' Issues with save as pdf, and some in svg for the Together plot
#' for none the color and size show the same thing
#library(tidyverse)
#library(cowplot)
#library(showtext)

MetaProVizLolipop <- function(Input_data, pCutoff= 0.05 , FCcutoff=0.5, test = "p.adj",OutputPlotName= "lolipop plot", plot_pathways = "none",# or "Individual" or "Together
                              Theme=theme_classic(), Save_as = svg ){
  
  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_Lolipop_folder = paste(Results_folder,"/Lolipop_Plots",  sep="")
  if (!dir.exists(Results_folder_plots_Lolipop_folder)) {dir.create(Results_folder_plots_Lolipop_folder)}
  
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  #Save_as= deparse(substitute(Save_as))
  
  
  # Remove rows with NAs
  Input_data<- Input_data %>% drop_na()
  #Select metabolites for the cut offs selected
  loli.data <- Input_data %>% mutate(names=Metabolite) %>% filter( abs(Log2FC) >=FCcutoff)
  
  
  if(plot_pathways == "Individual"){
    
    Pathway_Names <- unique(Input_data$Pathway)
    for (i in Pathway_Names){
      print(i)
      loli.data_path_indi <- loli.data %>% filter(Pathway==i)
      
      lolipop_plot <- ggplot(loli.data_path_indi , aes(x = Log2FC, y = names)) +
        geom_segment(aes(x = 0, xend = Log2FC, y = names, yend = names)) +
        geom_point(aes(colour = p.adj, size = p.adj ))   +
        scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +
        scale_colour_gradient(low = "red", high = "blue", limits = c(0, max(loli.data[,test]))) +
        ggtitle(label = paste(i," Pathway"), subtitle = paste("Metabolites with > |",FCcutoff,"| logfold change")) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Metabolites")+Theme
      #ggsave(filename = "Loli_plot2.pdf", plot = last_plot(), width=10, height=8)
      
      ggsave(file=paste(Results_folder_plots_Lolipop_folder, "/",i,"_", OutputPlotName, ".",Save_as, sep=""), plot=lolipop_plot, width=10, height=10)
      
      # svg(filename = paste("Results_", Sys.Date(), "/Lolipop_plots/",i,"_", OutputPlotName, ".svg", sep=""),
      #     width = 10,
      #     height = 8)
      # plot(lolipop_plot)
      # dev.off()
    }}else if(plot_pathways == "Together"){
      
      loli.data <- loli.data %>%
        arrange(Pathway, Metabolite)
      
      loli.data_avg <- loli.data %>%
        arrange(Pathway, Metabolite) %>%
        mutate(Metab_name = row_number()) %>%
        group_by(Pathway) %>%
        mutate(
          avg = mean(Log2FC)
        ) %>%
        ungroup() %>%
        mutate(Pathway = factor(Pathway))
      
      
      loli_lines <-   loli.data_avg %>%
        arrange(Pathway, Metabolite) %>%
        group_by(Pathway) %>%
        summarize(
          start_x = min(Metab_name) -0.5,
          end_x = max(Metab_name) + 0.5,
          y = 0#unique(avg)
        ) %>%
        pivot_longer(
          cols = c(start_x, end_x),
          names_to = "type",
          values_to = "x"
        ) %>%
        mutate(
          x_group = if_else(type == "start_x", x + .1, x - .1),
          x_group = if_else(type == "start_x" & x == min(x), x_group - .1, x_group),
          x_group = if_else(type == "end_x" & x == max(x), x_group + .1, x_group) )
      
      #rm(p2)
      p2 <- loli.data_avg %>%
        ggplot(aes(Metab_name, Log2FC)) + # names in aes ro Metab_name
        geom_hline(
          data = tibble(y = -5:5),
          aes(yintercept = y),
          color = "grey82",
          size = .5 )
      
      p2 <- p2 + geom_segment(
        aes(
          xend = Metab_name,          # names
          yend = 0,#avg,
          color = Pathway,
          #color = after_scale(colorspace::lighten(color, .2))
        ))
      
      p2 <- p2 + # geom_line( data = loli_lines, aes(x, y),  color = "grey40"  ) +
        geom_line(
          data = loli_lines,
          aes( x_group, y,
               color = Pathway,
               #  color = after_scale(colorspace::darken(color, .2))
          ), size = 2.5) +  geom_point(aes(size = p.adj, color = Pathway)
          )
      
      p2<- p2 + coord_flip()
      p2<-p2+ theme(axis.text.x=element_text())
      
      lab_pos_metab <- loli.data_avg %>% filter(Log2FC>0) %>% select(Metabolite, Metab_name, Log2FC)
      p2<- p2+ annotate("text", x = lab_pos_metab$Metab_name, y = lab_pos_metab$Log2FC+1.5, label = lab_pos_metab$Metabolite, size = 3)
      
      lab_neg_metab <- loli.data_avg %>% filter(Log2FC<0) %>% select(Metabolite, Metab_name, Log2FC)
      p2<- p2+ annotate("text", x = lab_neg_metab$Metab_name, y = lab_neg_metab$Log2FC-1.5, label = lab_neg_metab$Metabolite, size = 3)
      
      p2 <- p2+ annotate("text", x = max(lab_neg_metab$Metab_name)+ 3, y = 0, label = "Significantly changed metabolites and their pathways", size = 8)
      
      p2 <- p2+Theme
      # p2 + xlab("Metabolites") +  ggtitle("Titile is missing") # desnt work
      
      ggsave(file=paste(Results_folder_plots_Lolipop_folder, "/","Together", OutputPlotName, ".",Save_as, sep=""), plot=p2, width=20, height=20)
      # dev.off()
      # svg(filename = paste("Results_", Sys.Date(), "/Lolipop_plots/","Together", OutputPlotName, ".svg", sep=""),
      #     width = 14,
      #     height = 10)
      # plot(p2)
      # dev.off()
      
    }else if(plot_pathways == "none"){
      lolipop_plot <- ggplot(loli.data , aes(x = Log2FC, y = names)) +
        geom_segment(aes(x = 0, xend = Log2FC, y = names, yend = names)) +
        geom_point(aes(colour = p.adj, size = p.adj ))   +
        scale_size_continuous(range = c(1,5))+# , trans = 'reverse') +
        scale_colour_gradient(low = "red", high = "blue", limits = c(0, max(loli.data[,test]))) +
        ggtitle(paste("Metabolites with > |",FCcutoff,"| logfold change")) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Metabolites")+Theme
      #ggsave(filename = "Loli_plot2.pdf", plot = last_plot(), width=10, height=8)
      
      ggsave(file=paste(Results_folder_plots_Lolipop_folder, "/", OutputPlotName, ".",Save_as, sep=""), plot=lolipop_plot, width=10, height=10)
      
      # svg(filename = paste("Results_", Sys.Date(), "/Lolipop_plots/", OutputPlotName, ".svg", sep=""),
      #     width = 10,
      #     height = 8)
      # plot(lolipop_plot)
      # dev.off()
    }
  
}

###########----------------------
# Use function
#plotLolipop(Input_data = DMA_output, plot_pathways = "Together") # or Individual ot Together
######---------------------------


####################################
### ### ### Violin Plots ### ### ###
#####################################
#' @param Input_data Preprocessed data
#' @param Experimental_design Experimental_design
#' @param OutputPlotName File name for Output_plots=Together . When Output_plots = Individual the name is the name of each metabolite
#' @param Output_plots Option to save plots Individually or Together
#' @param Selected_Conditions Select which conditions will be plotted. IF nothing is selected then all conditions are plotted.
#' @param Selected_Comparisons Select between which comparison to take statistics.Works only Together with the Selected_Conditions. To select Comparisons add in a list
# the positions of the Conditions in the Selected_Conditions i.e. c(1,2) to take at-test between the 1st and 2nd Condition.
#
MetaProVizViolin <- function(Input_data,Experimental_design, OutputPlotName= "Violin", Output_plots = "Individual", # or "Together"
                             Selected_Conditions = NULL, Selected_Comparisons = NULL, Theme=theme_classic() , Save_as = svg){
  
  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_Violin_folder = paste(Results_folder,"/Violin_Plots",  sep="")
  if (!dir.exists(Results_folder_plots_Violin_folder)) {dir.create(Results_folder_plots_Violin_folder)}
  
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  #Save_as= deparse(substitute(Save_as))
  
  
  Metabolite_Names <- colnames(Input_data)
  
  # make a list for plotting all plots Together
  violin_plot_list <- list()
  k=1
  
  for (i in Metabolite_Names){
    
    
    violinplotdata <- Input_data %>%  select(i) %>%                         # Get mean & standard deviation by group
      group_by(Conditions=Experimental_design$Conditions)
    names(violinplotdata) <- c("Intensity", "Conditions")
    
    if (is.null(Selected_Conditions) == "FALSE"){
      violinplotdata <- violinplotdata %>% filter(Conditions %in% Selected_Conditions)
    }
    
    violinplot <- ggplot(violinplotdata, aes(x=Conditions, y=Intensity)) +
      geom_violin(fill="skyblue",width = 1) +
      geom_jitter(shape=16, position=position_jitter(0.2), alpha=1)
    
    l=1
    for (m in Selected_Comparisons){
      m <- unlist(m)
      violinplot <- violinplot+  stat_compare_means(data=violinplotdata, comparisons = list( c( Selected_Conditions[m][1],  Selected_Conditions[m][2])), method = "t.test", paired=TRUE, vjust = l)
      l=l+1
    }
    
    
    violinplot <- violinplot + Theme
    violinplot <- violinplot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
    violinplot <- violinplot + ggtitle(paste(i))
    violinplot
    
    if(Output_plots=="Individual"){
      
      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))
      
      ggsave(file=paste(Results_folder_plots_Violin_folder, "/", i, ".",Save_as, sep=""), plot=violinplot, width=10, height=10)
      
      # svg(filename = paste("Results_", Sys.Date(), "/Violin_plots/",i, ".", Save_as, sep=""),
      #     width = 10,
      #     height = 8)
      # plot(violinplot)
      # dev.off()
    } else if(Output_plots=="Together"){
      
      plot(violinplot)
      # save plot
      violin_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
    
  }
  
  pdf(file= paste(Results_folder_plots_Violin_folder,"/", OutputPlotName,".pdf", sep = ""), onefile = TRUE ) # or multivariate quality control chart
  for (plot in violin_plot_list){
    replayPlot(plot)
  }
  dev.off()
  
}

####---------------------
#plotViolin(Input_data = preprocessing_output$data_processed, Experimental_design = preprocessing_output$Experimental_design, Output_plots = "Individual" ,
#           Selected_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3)) )
#####----------------------




#################################
### ### ### Bargraphs ### ### ###
#################################
#; selected comparisons doesnt work

MetaProVizBarplots <- function(Input_data,Experimental_design, OutputPlotName= "Barplot", Output_plots = "Individual", # or "Together"
                               Selected_Conditions = NULL, Selected_Comparisons = NULL,Theme=theme_classic(), Save_as=svg){
  
  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_Barplots_folder = paste(Results_folder,"/Barplots_Plots",  sep="")
  if (!dir.exists(Results_folder_plots_Barplots_folder)) {dir.create(Results_folder_plots_Barplots_folder)}
  
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  #Save_as= deparse(substitute(Save_as))
  
  
  Metabolite_Names <- colnames(Input_data)
  
  # make a list for plotting all plots Together
  outlier_plot_list <- list()
  k=1
  
  for (i in Metabolite_Names){
    
    
    barplotdataMeans <- Input_data %>%  select(i) %>%                         # Get mean & standard deviation by group
      group_by(Conditions=Experimental_design$Conditions) %>%
      summarise_at(vars(i), list(mean = mean, sd = sd)) %>%
      as.data.frame()
    
    barplotdata <- Input_data %>%  select(i) %>%  group_by(Conditions=Experimental_design$Conditions)  %>%
      as.data.frame()
    names(barplotdata) <- c("Intensity", "Conditions")
    
    if (is.null(Selected_Conditions) == "FALSE"){
      barplotdataMeans <- barplotdataMeans %>% filter(Conditions %in% Selected_Conditions)
      barplotdata <- barplotdata %>% filter(Conditions %in% Selected_Conditions)
    }
    
    
    barplot<- ggplot(barplotdataMeans) +
      geom_bar( aes(x=Conditions, y=mean), stat="identity", fill="skyblue", alpha=0.7,width = 0.7) +
      geom_errorbar( aes(x=Conditions, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=0.5)
    
    # This is for selected compatrisons and it doesnt work
    #    l=1
    #     for (m in Selected_Comparisons){
    #       m <- unlist(m)
    #       stats <- compare_means(value ~ Intensity, group.by = "Conditions", data = barplotdata, method = "t.test")
    #       barplot <- barplot+  stat_compare_means(data=barplotdata, comparisons = list( c( Selected_Conditions[m][1],  Selected_Conditions[m][2])), method = "t.test", paired=TRUE)
    #       l=l+1
    #     }
    
    barplot <- barplot + Theme
    barplot <- barplot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
    barplot <- barplot + ggtitle(paste(i))
    
    
    
    if(Output_plots=="Individual"){
      
      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))
      
      ggsave(file=paste(Results_folder_plots_Barplots_folder, "/",i, ".",Save_as, sep=""), plot=barplot, width=10, height=8)
      
      # svg(filename = paste("Results_", Sys.Date(), "/Bar_plots/",i, ".svg", sep=""),
      #     width = 10,
      #     height = 8)
      # plot(barplot)
      # dev.off()
    } else if(Output_plots=="Together"){
      
      plot(barplot)
      # save plot
      outlier_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
    
  }
  
  pdf(file= paste(Results_folder_plots_Barplots_folder,"/", OutputPlotName,".pdf", sep = ""), onefile = TRUE ) # or multivariate quality control chart
  for (plot in outlier_plot_list){
    replayPlot(plot)
  }
  dev.off()
  
}

###########---------------------
#plotBargraphs(Input_data = preprocessing_output$data_processed, add_statistics = DMA_output, Experimental_design = preprocessing_output$Experimental_design )
# plotBarplots(Input_data = preprocessing_output$data_processed, Experimental_design = preprocessing_output$Experimental_design, Output_plots = "Together" )
#############--------------------






################################
### ### ### Boxplots ### ### ###
################################

plotBoxplots <- function(Input_data,
                         Conditions = "Conditions",
                         Experimental_design,
                        OutputPlotName= "Boxplot",
                         Output_plots = "Together", #"Individual" or "Together"
                         Selected_Conditions = NULL,
                         Selected_Comparisons = NULL,
                        STATtest = "t.test", # "wilcox.test"
                         Theme = theme_classic(),
                        Paired = FALSE,
                         Save_as = "svg"){
  
  ## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))
  
  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  RequiredPackages <- c("tidyverse", "ggplot2", "ggpubr")
  new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  suppressMessages(library(tidyverse))
  ##------------- Set the Conditions -----------##

  # Input_data = dplyr::select_if(Processed_data, is.numeric)
  # Experimental_design = Experimental_design
  # Output_plots = "Together"
  # OutputPlotName = "t.test"
  # Selected_Conditions = c("P", "H")
  # Conditions = "Pair"
  # Paired = TRUE
  # Selected_Comparisons = list(c(1,2))
  if (Conditions != "Conditions") {
    colnames(Experimental_design)[colnames(Experimental_design) == "Conditions"] = "Conditions2"
    colnames(Experimental_design)[colnames(Experimental_design) == Conditions] = "Conditions"
   }
  ## ------------ Check Input files ----------- ##
  #1. Input_data and Conditions
  if(any(duplicated(row.names(Input_data)))==TRUE){
    stop("Duplicated row.names of Input_data, whilst row.names must be unique")
  } else if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Experimental_design to obtain Condition1 and Condition2.")
  } else{
    #is this necessary?
    Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    #Test_num <- select_if(as.data.frame(Input_data), is.numeric)
    if((any(Test_num) ==  FALSE) ==  TRUE){
      stop("Input_data needs to be of class numeric")
    } else{
      Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
        stop("row.names Input_data need to match row.names Experimental_design.")
      } else(
        data <- Input_data
      )
    }
    Experimental_design <- Experimental_design
  }
  
  Output_plots_options <- c("Individual", "Together")
  if (Output_plots %in% Output_plots_options == FALSE){
    stop("Check Input the Plot_pathways option is incorrect. The Allowed options are the following: ",paste(Output_plots_options,collapse = ", "),"." )
  }
  if("Conditions" %in% colnames(Experimental_design)==FALSE){
    stop("There is no column named `Conditions` in Input_data.")
  }
  if(is.null(Selected_Conditions)==FALSE){
    #Conditions <- "CM4N_WT"
    for (Conditions in Selected_Conditions){
      if(Conditions %in% Experimental_design$Conditions==FALSE){
        stop("Check Input. The Selected_Conditions were not found in the Conditions Column.")
      }
    }
  }
  Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }
  
  ## ------------ Create Output folders ----------- ##
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  name <- paste0("Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  ### Create Volcano plots folder in  result directory ###
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} # Make Results folder
  Results_folder_plots_Boxplots_folder = file.path(Results_folder, "Boxplot")
  if (!dir.exists(Results_folder_plots_Boxplots_folder)) {dir.create(Results_folder_plots_Boxplots_folder)}  # check and create folder
  Results_folder_plots_Boxplots_folder_Analysis = file.path(Results_folder_plots_Boxplots_folder, Output_plots)
  if (!dir.exists(Results_folder_plots_Boxplots_folder_Analysis)) {dir.create(Results_folder_plots_Boxplots_folder_Analysis)}
  STATS.Test <- ""
  
  
  Metabolite_Names <- colnames(data)
  data <- merge(Experimental_design,data, by="row.names")
  data <- column_to_rownames(data, "Row.names")
  data <- data %>% arrange(factor(row.names(data), levels = row.names(Experimental_design)))
  # make a list for plotting all plots Together
  box_plot_list <- list()
  k=1
  #i<- "1-methylnicotinamide"
  for (i in Metabolite_Names){
    # data %>% select(i)
    boxplotdata <- data %>%  select(i, "Conditions") %>%                        # Get mean & standard deviation by group
      group_by(Conditions)
    names(boxplotdata) <- c("Intensity", "Conditions")
    
    if (is.null(Selected_Conditions) == "FALSE"){
      boxplotdata <- boxplotdata %>% filter(Conditions %in% Selected_Conditions)
    }
    
    #create a data.frame to be used with ggpubr
    if (Paired == TRUE){
      extract <- data %>%  select(i, "Conditions", "Patient") %>% as.data.frame()
     names(extract) <- c("Intensity", "Conditions", "ID")
     extractP <- extract[which(extract$Conditions == "P"), ] 
     extractH <- extract[which(extract$Conditions == "H"), ] 
     extract <- rbind(extractH,extractP)
    } else {
      extract <- data %>%  select(i, "Conditions") %>% as.data.frame()
      names(extract) <- c("Intensity", "Conditions")
    }
    #names(barplotdataMeans)[2] <- "Intensity"
    if(is.null(Selected_Comparisons)== TRUE){
      # Dimitris, those need to be clarified
      # names(barplotdataMeans)[2] <- "Intensity"
      # a <- max(barplotdataMeans$Intensity)
      boxplot <- ggplot(extract, aes(x=Conditions, y=Intensity, color = Conditions, palette = "jco")) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.4, size = 0.5) +
        ggplot2::geom_point(size = 2, alpha=0.7, colour = "darkred", shape = 16, position = position_jitter(width = 0.2, height  = 0.1, seed = 1234))+
        #ggplot2::geom_line(aes(group = factor(extract$ID)),colour="black", linetype=5, alpha = 0.7, position = position_jitter(width = 0.2, height  = 0.1, seed = 1234)) +
        #ggplot2::geom_jitter(shape=16, width = 0.2, height = 0.1, position = position_jitter(seed = 1234), alpha=0.7)+
        ggplot2::xlab("Conditions")+
        ggplot2::ylab("Mean Intensity")
      
    }else if (Paired == TRUE){
      #Theo, this line needs to be substituted by ggplot, see above.
      boxplot <- ggpubr::ggpaired(extract, x = "Conditions", y = "Intensity",
               color = "Conditions", line.color = "grey", line.size = 0.4,
                position = jitter)+
        ylab("Mean Intensity")+
        #geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.7)+
        ggpubr::stat_compare_means(paired = TRUE)
      #STATS.Test <-
        boxplot$layers
      }else{
      # names(barplotdataMeans)[2] <- "Intensity"
      # a <- max(barplotdataMeans$Intensity)
      boxplot <- ggplot(extract, aes(x=Conditions, y=Intensity)) +
        geom_boxplot(outlier.shape = NA, fill="skyblue") +
        geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6)+
        xlab("Conditions")+
        ylab("Mean Intensity")+
        ggpubr::stat_compare_means(comparisons = Selected_Comparisons,
                                   method = STATtest,
                                   label = "p.format", paired = FALSE, hide.ns = TRUE, position = position_dodge(0.9), vjust = 0.25, show.legend = FALSE) +
        theme(legend.position = "right")+xlab("Conditions")+ ylab("Mean Intensity")
      
      STATS.Test<- boxplot$layers[[3]]$stat_params$test
      }
    
    
    
    
    
    boxplot <- boxplot + Theme
    boxplot <- boxplot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
    if (Paired == TRUE){
      boxplot <- boxplot + ggtitle(paste(i, "p.val"))
    }else{
    boxplot <- boxplot + ggtitle(paste(i, ",",STATS.Test,"p.val"))}
    
    
    if(Output_plots=="Individual"){
      
      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))
      
      if(OutputPlotName ==""){
        ggsave(file=paste(Results_folder_plots_Boxplots_folder_Analysis, "/",i, ".",Save_as, sep=""), plot=boxplot, width=10, height=8)
      }else{
        ggsave(file=paste(Results_folder_plots_Boxplots_folder_Analysis, "/",OutputPlotName,"_",i, ".",Save_as, sep=""), plot=boxplot, width=10, height=8)
      }
      
    } else if(Output_plots=="Together"){
      
      plot(boxplot)
      # save plot
      box_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
  }
  if(Output_plots=="Together"){
    if(OutputPlotName ==""){
      pdf(file= paste(Results_folder_plots_Boxplots_folder_Analysis,"/Boxplots.pdf", sep = ""), onefile = TRUE )
    }else{
      pdf(file= paste(Results_folder_plots_Boxplots_folder_Analysis,"/Boxplots_", OutputPlotName,".pdf", sep = ""), onefile = TRUE )
    }
    for (plot in box_plot_list){
      replayPlot(plot)
    }
    dev.off()
  }
}

####---------------------

# plotBoxplots(Input_data = preprocessing_output$data_processed, Experimental_design = preprocessing_output$Experimental_design, Output_plots = "Individual",
#              ,Selected_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3)) )
#####----


###################################
### ### ### Super plots ### ### ###
###################################
#; This works only if you have biological replicates in the experimental design


plotSuperplots <- function(Input_data,Experimental_design,select_conditions, OutputPlotName= "Superplot", Output_plots = "Individual", Selected_Conditions = NULL,
                           Selected_Comparisons = NULL, Theme=theme_classic()){
  
  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  Results_folder = paste(getwd(), "/Results_",Sys.Date(),  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_Superplots_folder = paste(Results_folder,"/Superplots",  sep="")
  if (!dir.exists(Results_folder_plots_Superplots_folder)) {dir.create(Results_folder_plots_Superplots_folder)}
  
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  #Save_as= deparse(substitute(Save_as))
  
  #select the consitions
  Metabolite_Names <- colnames(Input_data)
  
  super_plot_list <- list()
  for (i in Metabolite_Names){
    
    cond_selected <- Input_data %>% select(i)
    names(cond_selected) <- "metabolite"
    cond_selected$Conditions <-  Experimental_design$Conditions
    cond_selected$Biological_Replicates <- Experimental_design$Biological_Replicates
    
    if (is.null(Selected_Conditions) == "FALSE"){
      cond_selected <- cond_selected %>% filter(Conditions %in% Selected_Conditions)
    }
    
    
    ReplicateAverages <- cond_selected %>%
      group_by(Conditions, Biological_Replicates) %>% summarise_each(list(mean))
    
    superplot<- ggplot(cond_selected, aes(x=Conditions,y=metabolite,color=factor(Biological_Replicates))) +
      geom_beeswarm(cex=1) + scale_colour_brewer(palette = "Set1") +
      geom_beeswarm(data=ReplicateAverages, size=4)
    
    l=1
    for (k in Selected_Comparisons){
      k <- unlist(k)
      superplot <- superplot+  stat_compare_means(data=ReplicateAverages, comparisons = list( c( Selected_Conditions[k][1],  Selected_Conditions[k][2])), method = "t.test", paired=TRUE, vjust = l)
      l=l+1
    }
    
    superplot<-  superplot+ theme(legend.position="right")
    superplot <- superplot + Theme
    superplot <- superplot + theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
    superplot <- superplot + ggtitle(paste(i)) +ylab(NULL)
    
    if(Output_plots=="Individual"){
      
      i <- (gsub("/","_",i))#remove "/" cause this can not be safed in a PDF name
      i <- (gsub(":","_",i))
      
      ggsave(file=paste(Results_folder_plots_Superplots_folder, "/",i, ".",Save_as, sep=""), plot=superplot, width=10, height=8)
      
      # svg(filename = paste("Results_", Sys.Date(), "/Super_plots/",i, ".svg", sep=""),
      #     width = 10,
      #     height = 8)
      # plot(superplot)
      # dev.off()
    } else if(Output_plots=="Together"){
      
      plot(superplot)
      # save plot
      super_plot_list[[k]] <- recordPlot()
      dev.off()
      k=k+1
    }
    
  }
  
  pdf(file= paste(Results_folder_plots_Superplots_folder,"/", OutputPlotName,".pdf", sep = ""), onefile = TRUE ) # or multivariate quality control chart
  for (plot in super_plot_list){
    replayPlot(plot)
  }
  dev.off()
  
}



########-------------------------
#plotSuperplots(Input_data = preprocessing_output$data_processed, Experimental_design = preprocessing_output$Experimental_design,
#               Output_plots = "Individual",Selected_Conditions = c("Control", "Rot", "3NPA"), Selected_Comparisons = list(c(1,2), c(1,3)) )


###############################
### ### ### Heatmap ### ### ###
###############################

plotHeatmaps<- function(Input_data,
						Experimental_design,
						Clustering_Condition = "Conditions",
						Clustering_method = "single",
						OutputPlotName= "Heatmap",
						Save_as=svg,		#Save_as="svg"
						kMEAN = c(5,10,15),
						SCALE = "row"){
  library("pheatmap")
  set.seed(12345)
  ####################################################
  # This searches for a Results directory within the current working directory and if its not found it creates a new one
  SysDATE <- Sys.Date()
  Results_folder = paste(getwd(), "/Results_",SysDATE,  sep="")
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)}
  ### Create Volcano plots folder in  result directory ###
  Results_folder_plots_Heatmaps_folder = paste(Results_folder,"/Heatmaps",  sep="")
  if (!dir.exists(Results_folder_plots_Heatmaps_folder)) {dir.create(Results_folder_plots_Heatmaps_folder)}
  
  #####################################################
  ### ### ### make output plot save_as name ### ### ###
  #Save_as= deparse(substitute(Save_as))

  my_annot<- NULL
  for (i in Clustering_Condition){
  my_annot[i] <- Experimental_design %>% select(i) %>% as.data.frame()
  }
  
  my_annot<- as.data.frame(my_annot)
  rownames(my_annot) <- rownames(Experimental_design)
  #Input_data = Processed_data
  Input_data <- as.data.frame(merge(dplyr::select_if(Input_data, is.numeric), Experimental_design,  by = "row.names"))
  row.names(Input_data) <- Input_data$Row.names
  Input_data <- Input_data[,-1]
  
  
  HEATMAP <- pheatmap(as.matrix(t(dplyr::select_if(Input_data, is.numeric))),
              clustering_method = Clustering_method,
              scale = "row",
              annotation_col = my_annot)
	ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", OutputPlotName, ".",Save_as, sep=""), plot=HEATMAP, width=10, height=12)
	
  #k=5
  #kMEAN=5
  #SCALE=column
	if (is.empty(kMEAN) == FALSE){
 for (k in kMEAN){
   out <-pheatmap(as.matrix(t(dplyr::select_if(Input_data, is.numeric))),
                  clustering_method =  "complete",
                  scale = SCALE,
                  kmeans_k = k,
                  clustering_distance_rows = "correlation",
                  annotation_col = my_annot)
   Results_folder_plots_Cluster_Analysis_folder= file.path(Results_folder_plots_Heatmaps_folder,paste0("MetaProVizplots_lineplots_k-means",k))
   if (!dir.exists(Results_folder_plots_Cluster_Analysis_folder)) {dir.create(Results_folder_plots_Cluster_Analysis_folder)}
   ggsave(file=paste0(Results_folder_plots_Heatmaps_folder,"/", OutputPlotName, "_kmeans-",k,".",Save_as), plot=out, width=10, height=12)
   Metabolite_clusters <- out[["kmeans"]][["cluster"]] %>% as.data.frame()
   names(Metabolite_clusters) <- "Clusters"
   
   Cluster_Analysis <- merge(t(Input_data), Metabolite_clusters, by = 'row.names' )
   names(Cluster_Analysis)[1] <- "Metabolite"
   Cluster_Analysis_selectec <- Cluster_Analysis 
   filter<-dplyr::filter
   writexl::write_xlsx(Cluster_Analysis_selectec, paste(Results_folder_plots_Cluster_Analysis_folder,"/","Clustering_k=",k,"_t(Input_data).xlsx", sep=""),col_names = TRUE)
   #cluster = 1
   #for Temperature in 
   for(cluster in unique(Cluster_Analysis$Clusters)){
     cluster1 <- Cluster_Analysis_selectec %>% filter(Clusters == cluster) %>% select(Metabolite)
     econ <- as.vector(cluster1$Metabolite)
     a <- Input_data %>% select(econ)
     b <- cbind(a, Experimental_design[which(Experimental_design$Temperature=="20oC"),])
     for (tissue in unique(Experimental_design$Tissue)){
       
       c <- b %>% filter(Tissue==tissue)
       d <- c %>% select(econ)
       d <- d %>% group_by(c$Species, c$Time) %>% summarise_at(vars(econ), list(name = mean))
       names(d)[1] <-"Species"
       names(d)[2]<- "Time"
       
       if (length(econ)==1){
         names(d)[3]<- econ[1]
         
       }
       
       #  library(MASS) 
       #library(reshape2) 
       d <- as.data.frame(d)
       e <- reshape2::melt(d, id = c("Species","Time")) 
       e$variable <-  gsub("_name", "",e$variable)
       #detach("package:reshape2", unload=TRUE)
       
       cluster1 <- ggplot(e, aes(x=Time, y=log(value,2), colour=variable, #shape = variable,
                                 group=interaction(Species, variable))) + 
         geom_point() + geom_line(linewidth=1.5) + 
         geom_text(data = subset(e, Time == "60mins"), aes(label = variable, x = Inf, y = log(value)), hjust = -.1) +
         # scale_colour_discrete(guide = 'none')  +    
         theme(plot.margin = unit(c(1,0.5,1,2), "lines"))   + theme(legend.position="right") + facet_grid(. ~Species)
       #cluster1
      
       ggsave(file=paste(Results_folder_plots_Cluster_Analysis_folder, "/Clustering_k=4_Time_Total_NMRat_Quad_logfc_important_Quad","_cluster_",cluster, ".svg", sep=""), plot=cluster1, width=22, height=8)
     }
   }
 }
}
  # selarate genes in # groups and get which genes are in which group
  # hc <- pl$tree_row
  # lbl <- cutree(hc, 5) # you'll need to change '5' to the number of gene-groups you're interested in
  # which(lbl==1) 
  
  ggsave(file=paste(Results_folder_plots_Heatmaps_folder,"/", OutputPlotName, ".",Save_as, sep=""), plot=heatmap, width=10, height=12)

  }

# https://www.biostars.org/p/287512/

#plotHeatmaps(Input_data = Quad_Dat, Experimental_design=Quad_Exp, Clustering_Consition = c("Conditions","Species"), OutputPlotName= "Heatmap_Quad")


