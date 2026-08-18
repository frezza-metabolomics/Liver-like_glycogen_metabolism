

Boxplot <- function(Result_Folder,
                    Input_data,
                    Input_SettingsFile,
                    Input_SettingsInfo,
                    OutputName=NULL,
                    FolderName = NULL ,
                    color=NULL
){
  
  library(ggpubr)
  library(tidyverse)
  if(is.null(color)==TRUE){
    boxdata <- merge(Input_SettingsFile %>% select(Conditions), Input_data, by=0) %>% column_to_rownames(., "Row.names")
  }else{
    boxdata <- merge(Input_SettingsFile %>% select(Conditions, color), Input_data, by=0) %>% column_to_rownames(., "Row.names")
  }
  

  
  Results_folder_Boxplot_folder = file.path(Result_Folder,  if (is.null(FolderName)) {paste("Boxplot") }else{ paste("Boxplot",FolderName,sep = "_" )}) # assign preprocessing folder
  
  if (!dir.exists(Results_folder_Boxplot_folder)) {dir.create(Results_folder_Boxplot_folder)}  # make preprocessing folder
  
  plot_list = NULL
  plot_list <-vector(mode = "list", length = ncol(boxdata)-1)
  names(plot_list)
  plot_list <- split(plot_list, names(boxdata[,-1]))
  message("Making Boxplots")
  
  if(is.null(color)==TRUE){
    for (i in 2:ncol(boxdata)){ # i=4
      #plotdata <- boxdata %>% select(c(i,1))
      #boxdata <- boxdata %>% filter(Conditions  %in% Input_SettingsInfo)
      #plotdata[,1] <- log(plotdata[,1],2)
      set.seed(1234)
      boxplot <- ggplot(boxdata[c(1,i)], aes(x = Conditions, y = boxdata[i][,1] )) +
        geom_boxplot(outlier.shape = NA, # Remove outlier points
                     alpha = 0.7, 
                     fill = "lightblue") +
        labs(title = paste0(colnames(boxdata[i])[1]),
             subtitle = "",
             x = "Conditions",
             y = "Normalised Intensity") +
        theme_minimal() +
        geom_jitter(#data = boxdata[i], aes(x = Conditions, y = boxdata[i][,1]),
          shape = 16, position = position_jitter(0.2), color = "black", size = 1.7) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Larger x-axis labels
              axis.text.y = element_text(size = 14),  # Larger y-axis labels
              axis.title.x = element_text(size = 16),  # Larger x-axis title
              axis.title.y = element_text(size = 16),  # Larger y-axis title
              plot.title = element_text(size = 18)) + # Larger plot title
        scale_y_continuous(limits = c(0, NA))#+
       # stat_compare_means(aes(label = ..p.signif..), method = "anova", label = "p.format", size = 5, label.x = 1.5) 
      plot_list[[colnames(boxdata[i])[1]]] <- boxplot
      plot(boxplot)
      plot_list[[colnames(boxdata[i])[1]]] <- recordPlot()
    }
    pdf(file= file.path(Results_folder_Boxplot_folder,paste0(OutputName,".pdf")), onefile = TRUE )
    for (plot in plot_list){
      replayPlot(plot)
    }
    dev.off()
    dev.off()
    
    
  }else{
    
    for (i in 3:ncol(boxdata)){
      plotdata <- boxdata %>% select(c(i,1,2))
      #plotdata[,1] <- log(plotdata[,1],2)
      set.seed(1234)
      
      # Install and load ggpubr package if not already installed
      #if (!require("ggpubr")) install.packages("ggpubr")
      
      # Create the boxplot with statistical annotations
      boxplot <- ggplot(boxdata[i], aes(x = Conditions, y = boxdata[i][,1], fill = Conditions)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        labs(title = paste0(colnames(boxdata[i])[1], " by Condition"),
             x = "Condition",
             y = "Normalised Intensity") +
        theme_minimal() +
        geom_jitter(aes(x = Conditions, y = boxdata[i][,1]),
                    shape = 16, position = position_jitter(0.2), color = "black", size = 1.7) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Larger x-axis labels
              axis.text.y = element_text(size = 14),  # Larger y-axis labels
              axis.title.x = element_text(size = 16),  # Larger x-axis title
              axis.title.y = element_text(size = 16),  # Larger y-axis title
              plot.title = element_text(size = 18),    # Larger plot title
              strip.text = element_text(size = 16)) +  # Larger facet label text
        # scale_y_continuous(limits = c(0, NA)) +
        facet_wrap(~ Data, scales = "free_x") +
        stat_compare_means(aes(label = ..p.signif..), method = "anova", label = "p.signif", size = 5, label.x = 1.5)
      
      #ggsave(filename =  paste0(Results_folder_Boxplot_folder, "/",colnames(plotdata)[1],"_.svg"), plot = boxplot, device = "svg", width = 10,  height = 10)
      plot(boxplot)
      plot_list[[colnames(boxdata[i])[1]]] <- recordPlot()
    }
    pdf(file= file.path(Results_folder_Boxplot_folder,paste0(OutputName,".pdf")), onefile = TRUE )
    for (plot in plot_list){
      replayPlot(plot)
    }
    dev.off()
    dev.off()
  }
}
