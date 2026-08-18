# Analysis of Naked Mole Rat project

#' Here we have run 3 batches
#' 1. AB01 : 200 samples
#' Mouse/ Naked Mole Rat
#' All under anoxia
#' Brain/ Kidney/ Heart/ Liver/ Muscle
#' 0/ 5/ 10/ 30/ 60 minutes time
#' 1/ 2/ 3/ 4/ replicates
#' 
#' 2. AB02: 153 samples
#' Mouse/ Naked Mole Rat
#' Normoxia/ Hypoxia10/ Hypoxia5(only for Nakes Mole Rat)
#' Brain/ Kidney/ Heart/ Liver/ Muscle
#' 1/ 2/ 3/ 4/ 5/ replicates
#' 
#' #' 3. AB06: 10 samples
#' Mouse/ Naked Mole Rat
#' All Heart
#' 0/ 30 mins time
#' NO inhibitor/ Amylase Inhibitor/ Bafilomycin/ Glycogen phos inhibitor (all for time 30 min)
#' 1 replicate
#' 

# #FOR WINDOWS USERS ONLY
# convert_clip <- function() {
#   gsub("\\\\", "/", readClipboard())
# }
# #copy the file path from the windows explorer
# pathtoFile <- convert_clip()
#pathtoFile <- "C:/Users/Dimitrios/Desktop/04. AB01 Analysis"

pathtoFile <- getwd()

### ### ### Load the data ### ### ###
source("scripts/00_Load_Data.R")
LoadData_result <-   LoadData(file.path(pathtoFile,"20220406_Amanda_AB01_pHILIC_E240.xlsx"), Sheet = 2)
outliers <- c("AB01_069", "AB01_056")
LoadData_result <- lapply(LoadData_result, function(x) x %>% filter(!rownames(x) %in% outliers))

################################################################################
########################## Analysis start ######################################
################################################################################

# make the experimental design as it should be
Experimental_design <- LoadData_result$Experimental_design

Experimental_design$Conditions <- str_replace(Experimental_design$Conditions, "_0mins", "_00mins")
Experimental_design$Conditions <- str_replace(Experimental_design$Conditions, "_5mins", "_05mins")
Experimental_design$Conditions <- str_replace(Experimental_design$Conditions, "20mins", "30mins")

Experimental_design <- cbind(Experimental_design, separate(Experimental_design, col=Conditions,sep="_", into = c('Species', 'Tissue', "Treatment", "Time") ) )
Experimental_design$Replicate <- substring(Experimental_design$Species, 3,3)
Experimental_design$Species <- substring(Experimental_design$Species, 1,2)
Experimental_design$Conditions_Rep <- paste(Experimental_design$Species, Experimental_design$Tissue, Experimental_design$Time, sep = "_")

names(Experimental_design)[1] <- "Full_name"
names(Experimental_design)[7] <- "Conditions"


#####################################
### ### ### Preprocessing ### ### ###
#####################################
source("scripts/01_MetaProVizPreprocessing.R")

Preprocessing_result <- MetaProVizPreprocessing(Input_data = LoadData_result$Data,
                                                Experimental_design = Experimental_design,
                                                Normalization = "TIC",
                                                PQN_selected_persentage=100,
                                                Standards = c("valine-d8 (IS)", "hippurate-d5 (IS)"))

Processed_data <- Preprocessing_result$Processed_data

source("scripts/04_MetaProVizVisualization.R")
#################################################################
### ### ###  PCA plots of all the tissues separately  ### ### ###
# Get Brain data
LoadData_result <- lapply(LoadData_result, function(x) x %>% filter(!rownames(x) %in% outliers))

#liver PCA
Liver_Exp <- Experimental_design %>% filter(Tissue =="Liver")     
Liver_Dat <- Processed_data %>% filter(row.names(Processed_data) %in% row.names(Liver_Exp))     
row.names(Liver_Dat) == row.names(Liver_Exp)
# Remove zero variance variables
#Liver_Dat <- Liver_Dat[, sapply(Liver_Dat, var) != 0]
MetaProVizPCA(data= Liver_Dat, Design= Liver_Exp, Color = "Time", Shape = "Species", Show_Loadings = FALSE,  Scaling = TRUE, 
              OutputPlotName= 'Liver_Colour_Species_Shape_Time')

#Heart_PCA
Heart_Exp <- Experimental_design %>% filter(Tissue =="Heart")     
Heart_Dat <-Processed_data %>% filter(row.names(Processed_data) %in% row.names(Heart_Exp))        

# Remove zero variance variables
Heart_Dat <- Heart_Dat %>%
  select(where(~ !is.numeric(.) || var(., na.rm = TRUE) != 0))
MetaProVizPCA(data = Heart_Dat, Design = Heart_Exp, Color = "Time", Shape = "Species", Show_Loadings = FALSE,  Scaling = TRUE, 
              OutputPlotName = 'Heart_Colour_Species_Shape_Time')

# Heart and Liver PCA
library(ggConvexHull)
#heart and liver
Heart_Liver_Exp <- Experimental_design %>% filter(Tissue =="Heart" | Tissue =="Liver")     
Heart_Liver_Dat <-Processed_data %>% filter(row.names(Processed_data) %in% row.names(Heart_Liver_Exp))      
row.names(Heart_Liver_Dat) == row.names(Heart_Liver_Exp)
#Heart_Liver_Dat <- Heart_Liver_Dat[, sapply(Heart_Liver_Dat, var) != 0]
MetaProVizPCA(data = Heart_Liver_Dat, Design = Heart_Liver_Exp, Color = "Tissue", Shape = "Time",
              Fill = "Species", Show_Loadings = FALSE,  Scaling = TRUE, 
              OutputPlotName = 'Heart_Liver_Colour_Tissue_Shape_Time_Shape_Species')

#######################################################################
# Heatmaps
# Individual Tissues
plotHeatmaps(Input_data = Heart_Dat, Experimental_design= Heart_Exp, Clustering_Condition = c("Time","Species"), OutputPlotName= "Heatmap_Heart")
plotHeatmaps(Input_data = Liver_Dat, Experimental_design= Liver_Exp, Clustering_Condition = c("Time","Species"), OutputPlotName= "Heatmap_Liver")

###############################################################################
# boxplots
Heart_Dat_num <- Heart_Dat %>%
  select(where(is.numeric))

plotBoxplots(Input_data = Heart_Dat_num,Experimental_design = Heart_Exp, OutputPlotName= "Boxplots_Heart",# out_plots = "together", # or "together"
             Selected_Conditions = NULL, Selected_Comparisons = NULL, Theme=theme_classic(), Save_as="svg" )


#############################################################################
# clustering with heatmap
#Heart_PCA
Heart_Exp <- Experimental_design %>% filter(Tissue =="Heart")     
Heart_Dat <-Processed_data %>% filter(row.names(Processed_data) %in% row.names(Heart_Exp))  
Heart_Dat <- Heart_Dat[,-c(1:8)]


##WE get the Heatmap data
library(pheatmap)

test_data <- Heart_Dat 
test_SettingsFile<- Heart_Exp

test_data <- test_data[, sapply(test_data, function(x) var(x, na.rm = TRUE) > 0)]

# Function to scale a column from -1 to 1
scale_to_minus1_to_1 <- function(x) {
  (2 * (x - min(x)) / (max(x) - min(x))) - 1
}

# Apply the scaling function to all columns
test_scaled <- as.data.frame(lapply(test_data, scale_to_minus1_to_1))
rownames(test_scaled) <- rownames(test_data)
colnames(test_scaled) <- colnames(test_data)

#test_scaled <-scale(test_data)
#test_scaled <-test_data
Heart_Exp$Conditions <- factor(Heart_Exp$Conditions, levels = unique(Heart_Exp$Conditions))

test_SettingsFile_ordered <- Heart_Exp %>% arrange(Conditions)
test_scaled_ordered <- test_scaled[rownames(test_SettingsFile_ordered),]


library(cluster)
library(factoextra)
library(dplyr)

set.seed(12345)
# Function to compute silhouette width for different values of k
compute_silhouette <- function(scaled_data, k_max) {
  silhouette_scores <- numeric(k_max - 1)
  
  for (k in 2:k_max) {
    kmeans_result <- kmeans(scaled_data, centers = k, nstart = 25)
    sil <- silhouette(kmeans_result$cluster, dist(scaled_data, method = "euclidean"))
    silhouette_scores[k - 1] <- mean(sil[, 3])  # Extract mean silhouette width
  }
  
  # Plot Silhouette scores
  plot(2:k_max, silhouette_scores, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", ylab = "Average Silhouette Width",
       main = "Silhouette Method for Optimal k Selection")
  points(14, silhouette_scores[14 - 2 + 1], col = "red", pch = 19, cex = 1.5)
  
  
  # Return the optimal k
  optimal_k <- which.max(silhouette_scores) + 1
  return(optimal_k)
}

# Define maximum k to test
k_max <- 25  # You can adjust this based on your dataset

# Run the silhouette method to find optimal k
optimal_k <- compute_silhouette(test_scaled_ordered, k_max)

cat("Optimal number of clusters based on silhouette method:", optimal_k, "\n")



set.seed(12345)
out <- pheatmap(t(test_scaled_ordered),
                clustering_method =  "complete", # "complete"
                # scale = "row",
                kmeans_k = 14,
                cluster_cols = FALSE,
                clustering_distance_rows = "correlation",
                # clustering_distance_cols =  "none",
                main = " Mouse and HG",
                annotation_col = test_SettingsFile_ordered %>% select(2,5)
)
ggsave(file=paste("Clustering_k=14.xlsx.svg", sep=""), plot=out, width=12, height=8)



Metabolite_clusters <- out[["kmeans"]][["cluster"]] %>% as.data.frame()
names(Metabolite_clusters) <- "Clusters"
Metabolite_clusters <- rownames_to_column(Metabolite_clusters, "Feature")


writexl::write_xlsx(Metabolite_clusters, paste("Clustering_k=14.xlsx", sep=""),col_names = TRUE)


library(dplyr)
library(tidyr)
library(pheatmap)

Cluster_heatmaps_folder <- file.path(
  "Cluster_heatmaps"
)

if (!dir.exists(Cluster_heatmaps_folder)) {
  dir.create(Cluster_heatmaps_folder, recursive = TRUE)
}

for(cluster in unique(Metabolite_clusters$Clusters)) {
  
  econ <- Metabolite_clusters %>%
    filter(Clusters == cluster) %>%
    pull(Feature)
  
  cluster_data <- Heart_Dat %>%
    select(all_of(econ)) %>%
    mutate(
      Species = Heart_Exp$Species,
      Time = Heart_Exp$Time
    )
  
  cluster_mean <- cluster_data %>%
    group_by(Species, Time) %>%
    summarise(
      across(all_of(econ), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  heatmap_data <- cluster_mean %>%
    mutate(
      Group = paste(Species, Time, sep = "_")
    ) %>%
    select(-Species, -Time) %>%
    pivot_longer(
      cols = -Group,
      names_to = "Metabolite",
      values_to = "Value"
    ) %>%
    pivot_wider(
      names_from = Group,
      values_from = Value
    ) %>%
    as.data.frame()
  
  rownames(heatmap_data) <- heatmap_data$Metabolite
  heatmap_data$Metabolite <- NULL
  
  # desired column order
  wanted_order <- c(
    "MM_00mins", "MM_05mins", "MM_10mins", "MM_30mins", "MM_60mins",
    "HG_00mins", "HG_05mins", "HG_10mins", "HG_30mins", "HG_60mins"
  )
  
  heatmap_data <- heatmap_data[, wanted_order, drop = FALSE]
  
  # Heatmap
  svg(
    filename = file.path(
      Cluster_heatmaps_folder,
      paste0("Heart_Cluster_", cluster, ".svg")
    ),
    width = 10,
    height = max(5, nrow(heatmap_data) * 0.25)
  )
  
  pheatmap(
    heatmap_data,
    scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "grey60",
    main = paste("Cluster", cluster)
  )
  
  dev.off()
}


##############################################################################
# Differential abundance analysis


############################################
# 1a) Naked Mole rat -vs- Mouse time zero all tissues all timepoints
# Prepare the data
Processed_data_filtered <- Processed_data#[-which(rownames(Processed_data) == "AB01_069" | rownames(Processed_data) =="AB01_056"),]
Experimental_design_filtered <- Experimental_design#[-which(rownames(Processed_data) == "AB01_069" | rownames(Processed_data) =="AB01_056"),]


DMA_data <-  Processed_data_filtered
DMA_Exp_design_Species <- Experimental_design_filtered %>% select(Species) %>% rename("Conditions" = "Species")

source("scripts/02_MetaProVizDifferentialMetaboliteAnalysis.R")
DMA_species_all_tissue_res <-MetaProVizDMA(Input_data =DMA_data, Experimental_design= DMA_Exp_design_Species,
                                           Condition1 = "MM", Condition2 = "HG", OutputName= 'all_tissues_together_all_timepoints_together')

### Make experimetal design with species and time information
DMA_Exp_design_Species_Time <- DMA_Exp_design_Species
DMA_Exp_design_Species_Time[1] <- paste(Experimental_design_filtered$Species, Experimental_design_filtered$Time, sep="_")

# DMA Naked Mole rat -vs- Mouse time zero all tissues  time 00mins
DMA_species_all_tissue_time_00mins_res <- MetaProVizDMA(Input_data =DMA_data, Experimental_design= DMA_Exp_design_Species_Time,
                                                        Condition1 = "MM_00mins", Condition2 = "HG_00mins", OutputName= 'all_tissues_together_time_00mins')

# DMA Naked Mole rat -vs- Mouse time zero all tissues  time 05mins
DMA_species_all_tissue_time_05mins_res <- MetaProVizDMA(Input_data =DMA_data, Experimental_design= DMA_Exp_design_Species_Time,
                                                        Condition1 = "MM_05mins", Condition2 = "HG_05mins", OutputName = 'all_tissues_together_time_05mins')

# DMA Naked Mole rat -vs- Mouse time zero all tissues  time 10mins
DMA_species_all_tissue_time_10mins_res <- MetaProVizDMA(Input_data =DMA_data, Experimental_design= DMA_Exp_design_Species_Time,
                                                        Condition1 = "MM_10mins", Condition2 = "HG_10mins", OutputName = 'all_tissues_together_time_10mins')

# DMA Naked Mole rat -vs- Mouse time zero all tissues  time 30mins
DMA_species_all_tissue_time_30mins_res <- MetaProVizDMA(Input_data =DMA_data, Experimental_design= DMA_Exp_design_Species_Time,
                                                        Condition1 = "MM_30mins", Condition2 = "HG_30mins", OutputName = 'all_tissues_together_time_30mins')

# DMA Naked Mole rat -vs- Mouse time zero all tissues  time 60mins
DMA_species_all_tissue_time_60mins_res <- MetaProVizDMA(Input_data =DMA_data, Experimental_design= DMA_Exp_design_Species_Time,
                                                        Condition1 = "MM_60mins", Condition2 = "HG_60mins", OutputName = 'all_tissues_together_time_60mins')




###############################################################################
###------------------------------ END --------------------------------------###



