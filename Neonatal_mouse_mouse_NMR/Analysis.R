
library(dplyr)
file <- read.csv("AB01_AB08_Saccharides_Valine-d8.csv", sep = ",", na.strings = NULL)
samples <- read.csv("samples.txt", sep = "\t", na.strings = NULL, header = FALSE)
length(rownames(file))
length(rownames(file[1]))
rownames(file) <- file$X 
file <- file[-1]
file <- cbind(file, Batch=c(rep(1,32),rep(2,200)))
file <- cbind(file, Tissue=c(rep("Heart",16),rep("Liver",16),rep(c(rep("Heart",5),rep("Liver",5),rep("Brain",5),rep("Kidney",5),rep("Quad",5)),8)))
file <- cbind(file, Spieces=c(rep("Mouse",32),rep("Mouse",100),rep("MoleRat",100)))

rownames(samples) <- samples$V1
samples <- samples[-1]
colnames(samples) <- c("Tissue", "Time")
rownames(file) <-  rownames(samples)
merGED <- merge(file, samples, 
      by = c('row.names', "Tissue"), all = TRUE,
      sort = FALSE) 
rownames(merGED) <- merGED$Row.names
merGED <- merGED[-1]
merGED$Time <- gsub("min","",merGED$Time)
merGED <- merGED %>% relocate(Tissue, .after=Time)
merGED$Tetrasaccharide<- as.numeric(merGED$Tetrasaccharide+1)
write.csv(merGED, "merGEDraw.csv")
write.csv(samples, "samples.csv")
write.csv(file, "file.csv")



#install.packages("ggstatsplot")
# library(ggstatsplot)
# library(tidyverse)
GATHER <- as.data.frame(reshape::melt(merGED[-4], id = c("Tissue", "Spieces", "Batch", "Time")))
GATHERVal <- as.data.frame(reshape::melt(merGED[-c(1:3)], id = c("Tissue", "Spieces", "Batch", "Time")))
# ggbetweenstats(data = GATHER,
#   x = Tissue,
#   y = value       
# )
# ggbetweenstats(data = GATHER,
# x = variable,
# y = value)
library(ggplot2)
library(viridis)
summary(GATHER)
GATHERlog <- cbind(GATHER[1:5], log2(GATHER[6]))
as_fct_innatural <- function(x, ordered= TRUE) factor(x, stringr::str_sort(unique(x), numeric = TRUE), ordered = ordered)
vars(as_fct_innatural(Time))
ggplot(GATHERlog, aes(x=variable, y=value)) + 
  geom_violin(trim=FALSE) +
  geom_jitter(aes(x=variable, y=value, color = Time, shape= Spieces))  + 
  scale_color_viridis(discrete = TRUE, option = "D")+
  facet_grid(Tissue ~ as_fct_innatural(Time))



#B
#Normalization by Median (Batch-wise):
#Calculate the median Valine-d8 intensity for each batch and divide each metabolite intensity within that batch by the corresponding median. This can be helpful if batch effects are suspected.

set.seed(1234)
data_norm_batch_median <- subset(merGED, Tissue == c("Heart","Liver")) %>%
  dplyr::group_by(Batch) %>%
  dplyr::mutate_at(vars(Disaccharide:Trisaccharide), ~ .x / median(valine.d8)) %>% as.data.frame()
rownames(data_norm_batch_median)
rownames(data_norm_batch_median) <- rownames(subset(merGED, Tissue == c("Heart","Liver")))
GATHERmedian <- as.data.frame(reshape::melt(as.data.frame(data_norm_batch_median)[-4], id = c("Tissue", "Spieces", "Batch", "Time")))

set.seed(1234); ggplot(subset(GATHERmedian, Tissue == c("Liver","Heart")), aes(x=variable, y=value)) + 
  geom_violin(trim=FALSE) +
  geom_jitter(aes(x=variable, y=value, color = factor(Batch), shape= factor(Tissue)))#+
  

subsetED <- subset(merGED, Tissue == c("Heart","Liver"))


#rAW
merged_Heart_Liver <- subset(merGED, Tissue == c("Heart","Liver"))
GATHER_Heart_Liver <- as.data.frame(reshape::melt(merged_Heart_Liver[-4], id = c("Tissue", "Spieces", "Batch", "Time")))

#Median by batch
data_norm_batch_median <- merged_Heart_Liver %>%
  dplyr::group_by(Batch) %>%
  dplyr::mutate_at(vars(Disaccharide:Trisaccharide), ~ .x / median(valine.d8)) %>% as.data.frame()
rownames(data_norm_batch_median)
rownames(data_norm_batch_median) <- rownames(merged_Heart_Liver)
GATHERmedian <- as.data.frame(reshape::melt(data_norm_batch_median[-4], id = c("Tissue", "Spieces", "Batch", "Time")))


data_norm_batch_median$Spieces[grepl("^AB08", rownames(data_norm_batch_median))] <- "Neonatal"
data_norm_batch_median[, 1:3] <- log2(data_norm_batch_median[, 1:3])

ggpubr::ggline(data_norm_batch_median[-4],
               x = "Time", y = c("Disaccharide","Tetrasaccharide","Trisaccharide"), 
               add = c("mean", "jitter"),
               color = "Spieces", palette = "jco",
               facet.by = "Tissue",  order = c("0","5","10", "20", "30", "60"),
               xlab = NULL,
               ylab = "Ratio normalised Intensity (A.U.)")#




