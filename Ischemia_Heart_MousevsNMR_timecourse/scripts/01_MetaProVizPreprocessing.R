## ---------------------------
##
## Script name: MetaProVizPreprocessing
##
## Purpose of script: Metabolomics (raw ion counts) pre-processing, normalisation and outlier detection
##
## Author: Dimitrios Prymidis and Christina Schmidt
##
## Date Created: 2022-10-28
##
## Copyright (c) Dimitrios Prymidis and Christina Schmidt
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

#' MetaProVizPreprocessing
#'
#' Applies 80%-filtering rule, total-ion count normalisation, missing value imputation and HotellingT2 outlier detection
#'
#' @param Input Data matrix which contains unique sample identifiers as row names and metabolite numerical values in columns with metabolite identifiers as column names. Use NA for metabolites that were not detected.
#' @param Experimental_design Data matrix which contains information about the samples, which will be combined with your input data based on the unique sample identifiers used as rownames. Column "Conditions" with information about the sample conditions (e.g. "N" and "T" or "Normal" and "Tumor"), can be used for feature filtering and colour coding in the PCA. Column "AnalyticalReplicate" including numerical values, defines technical repetitions of measurements, which will be summarised. Column "BiologicalReplicates" including numerical values. Please use the following names: "Conditions", "Biological_Replicates", "Analytical_Replicates"
#' @param Feature_Filtering \emph{Optional: }If set to "None" then no feature filtering is performed. If set to Standard then it applies the 80%-filtering rule (Bijlsma S. et al., 2006) on the metabolite features on the whole dataset. If is set to "Modified",filtering is done based on the different conditions, thus a column named "Conditions" must be provided in the Experimental_design including the individual conditions you want to apply the filtering to (Yang, J et al., 2015).\strong{Default=TRUE} \strong{Default=Modified}
#' @param Feature_Filt_Value \emph{Optional: } Percentage of feature filtering (Bijlsma S. et al., 2006).\strong{Default=0.8}
#' @param Normalization \emph{Optional: } If TIC, Total Ion Count normalization is performed. \strong{Default=TRUE}
#' @param HotellinsConfidence \emph{Optional: } Defines the Confidence of Outlier identification in HotellingT2 test. Must be numeric.\strong{Default = 0.99}
#' @param ExportQCPlots \emph{Optional: } Select whether the quality control (QC) plots will be exported. \strong{Default = TRUE}
#' @param Filter_Potential_Outliers \emph{Optional: } If TRUE, potential outliers between confidence of 0.95 and 0.99 are filtered out.\strong{Default=FALSE}
#' @param CoRe \emph{Optional: } If TRUE, a consumption-release experiment has been performed and the CoRe value will be calculated. Please consider providing a Normalisation factor column called "CoRe_norm_factor" in your "Experimental_design" DF, where the column "Conditions" matches. Th normalisation factor must be a numerical value obtained from growth rate that has been obtained from a growth curve or growth factor that was obtained by the ratio of cell count/protein quantification at the start point to cell count/protein quantification at the end point.. Additionally control media samples have to be available in the "Input" DF and defined as "blank" samples in the "Conditions" column in the "Experimental_design" DF, e.g. "blank_1", "blank_2". \strong{Default=FALSE}
#' @param Save_as \emph{Optional: } Select the file type of output plots. Options are svg, png, pdf, jpeg, tiff, bmp. \strong{Default=svg}
#'
#' @keywords 80% filtering rule, Missing Value Imputation, Total Ion Count normalization, PCA, HotellingT2, multivariate quality control charts,
#' @export


###################################################
### ### ### Metabolomics pre-processing ### ### ###
###################################################
#Input_data=test
#Experimental_design = With_pooled
MetaProVizPreprocessing <- function(Input_data=data,
						Experimental_design = design,
						Feature_Filtering = "Modified",
						Feature_Filt_Value = 0.8,
						Normalization = "TIC",
						PQN_selected_persentage = 100,
						Filter_Potential_Outliers = FALSE,
						#theo#the hottelinsConfidence
						HotellinsConfidence = 0.99,
                        ExportQCPlots = TRUE,
						CoRe = FALSE,
						OutlierLoop = 2,
						npcs=NULL,
						Standards = c("valine-d8", "hippurate-d5"),
						Save_as = "svg"
						){

## ------------ Setup and installs ----------- ##
  RequiredPackages <- c("tidyverse", # general scripting
                        "factoextra", # visualize PCA
                        "qcc", # for hotelling plots
                        "ggplot2", # For visualization PCA
                        "hash", # Dictionary in R for making column of outliers
                        "inflection")# For finding inflection point/ Elbow knee /PCA component selection # https://cran.r-project.org/web/packages/inflection/inflection.pdf # https://deliverypdf.ssrn.com/delivery.php?ID = 454026098004123081018105104090015093000085002012023032095093077109069092095000114006057018122039107109012089110120018031068078025094036037013095100070100076109026029024044005068010070117123085122016083112098002109001027028000024115096122101001083084026&EXT = pdf&INDEX = TRUE # https://arxiv.org/abs/1206.5478
# Load libraries
new.packages <- RequiredPackages[!(RequiredPackages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(library(tidyverse)) # general scripting
#tidyverse::tidyverse_packages()
#library(ggplot2) # For visualization PCA
# mutate_all <- dplyr::mutate_all
# select_if <- dplyr::select_if
# `%>%` <- magrittr::`%>%`
# sort <- base::sort
# ggplot <- ggplot2::ggplot
# geom_point <- ggplot2::geom_point
# aes <- ggplot2::aes
# theme_classic <- ggplot2::theme_classic
# geom_hline <- ggplot2::geom_hline
# geom_vline <- ggplot2::geom_vline
# geom_text <- ggplot2::geom_text
# ggtitle <- ggplot2::ggtitle
# scale_x_continuous <- ggplot2::scale_x_continuous
# scale_y_continuous <- ggplot2::scale_y_continuous
# theme <- ggplot2::theme
# geom_line <- ggplot2::geom_line
# scale_linetype_discrete <- ggplot2::scale_linetype_discrete
# ggsave <- ggplot2::ggsave
# annotate <- ggplot2::annotate
# ggplot_build<- ggplot2::ggplot_build
# 
# library(factoextra) # visualize PCA
# library(qcc) # for hotelling plots
# library(hash) # Dictionary in r. for making column of outliers
# library(inflection) # For finding inflection point/ Elbow knee /PCA component selection # https://cran.r-project.org/web/packages/inflection/inflection.pdf # https://deliverypdf.ssrn.com/delivery.php?ID=454026098004123081018105104090015093000085002012023032095093077109069092095000114006057018122039107109012089110120018031068078025094036037013095100070100076109026029024044005068010070117123085122016083112098002109001027028000024115096122101001083084026&EXT=pdf&INDEX=TRUE # https://arxiv.org/abs/1206.5478

## ------------------ Run ------------------- ##

#######################################################################################
### ### ### Check Input Information and add Experimental_design information ### ### ###

  #1.  Input data
  if(any(duplicated(row.names(Input_data))) ==  TRUE){# Is the "Input_data" has unique IDs as row names and numeric values in columns?
	stop("Duplicated row.names of Input_data, whilst row.names must be unique")
	} else{
	Test_num <- apply(Input_data, 2, function(x) is.numeric(x))
    if((any(Test_num) ==  FALSE) ==  TRUE){
			stop("Input_data needs to be of class numeric")
		} else {
		Test_match <- merge(Experimental_design, Input_data, by.x = "row.names",by.y = "row.names", all =  FALSE) # Do the unique IDs of the "Input_data" match the row names of the "Experimental_design"?
      if(nrow(Test_match) ==  0){
			stop("row.names Input_data need to match row.names Experimental_design.")
		} else(
        Input_data <- Input_data
			)
		}
	}

  #2. Conditions
  if ( "Conditions" %in% colnames(Experimental_design)){   # Parse Condition and Replicate information
	Conditions <- Experimental_design$Conditions
	}else{
	Conditions <- NULL
	}
if("Biological_Replicates" %in% colnames(Experimental_design)){
	Biological_Replicates <- Experimental_design$Biological_Replicates
	}else{
	Biological_Replicates <-NULL
	}
if("Analytical_Replicates" %in% colnames(Experimental_design)){
	Analytical_Replicates <- Experimental_design$Analytical_Replicates
	}else{
	Analytical_Replicates <-NULL
	}

  #3. Core parameters
  if (CoRe ==  TRUE){   # parse CoRe normalisation factor
	message("For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.")
    if ("CoRe_norm_factor" %in% colnames(Experimental_design)){
		CoRe_norm_factor <- Experimental_design$CoRe_norm_factor
	}else{
		warning("No growth rate or growth factor provided for normalising the CoRe result, hence CoRe_norm_factor set to 1 for each sample")
		CoRe_norm_factor <- as.numeric(rep(1,length(Experimental_design$Conditions)))
	}
    if (length(CoRe_norm_factor) !=  length(Experimental_design$Conditions)){ # Check is the length of normalization factor and conditions in 1 to 1
	stop("The CoRe_norm_factor length is different from the amount of samples. Please input a vector with a value for each sample. Blanks should take a value of 1.")
	}
    if(length(grep("blank", Conditions)) < 1){     # Check for blank samples
      stop("No blank samples were provided in the 'Conditions' in the Experimental design'. For a CoRe experiment control media samples without cells have to be measured and be added in the 'Conditions'
           column labeled as 'blank' (see @param section). Please make sure that you used the correct labelling or whether you need CoRE = FALSE for your analysis")
	}
}

#4. General parameters
  if(Feature_Filtering %in% c("Standard","Modified", "none") == FALSE ){
    stop("Check input. The selected Feature_Filtering option is not valid. Please select one of the folowwing: ",paste(Feature_Filtering_options,collapse = ", "),"." )
  }
  if( is.numeric(Feature_Filt_Value) == FALSE |Feature_Filt_Value > 1 | Feature_Filt_Value < 0){
    stop("Check input. The selected Filtering value should be numeric and between 0 and 1.")
  }
    if( is.numeric(HotellinsConfidence)== FALSE |HotellinsConfidence > 1 | HotellinsConfidence < 0){
    stop("Check input. The selected Filtering value should be numeric and between 0 and 1.")
  }
    if(is.logical(ExportQCPlots) == FALSE){
    stop("Check input. The ExportQCPlots value should be either =TRUE if QC plots are to be exported or =FALSE if not.")
  }
   if(is.logical(CoRe) == FALSE){
    stop("Check input. The CoRe value should be either =TRUE for preprocessing of Consuption/Release experiment or =FALSE if not.")
  }
   Save_as_options <- c("svg","pdf", "jpeg", "tiff", "png", "bmp", "wmf","eps", "ps", "tex" )
  if(Save_as %in% Save_as_options == FALSE){
    stop("Check input. The selected Save_as option is not valid. Please select one of the folowwing: ",paste(Save_as_options,collapse = ", "),"." )
  }
  
Input_data <- as.matrix(mutate_all(as.data.frame(Input_data), function(x) as.numeric(as.character(x))))

####################################################
### ### ### Create output folders ### ### ###
# This searches for a folder called "Results" within the current working directory and if its not found it creates one
name <- paste0("Results_",Sys.Date())
  WorkD <- getwd()
  Results_folder <- file.path(WorkD, name)
  #Make Results_folder
  if (!dir.exists(Results_folder)) {dir.create(Results_folder)} 
# This searches for a folder called "Preprocessing" within the "Results" folder in the current working directory and if its not found it creates one
Results_folder_Preprocessing_folder = file.path(Results_folder, "Preprocessing")
if (!dir.exists(Results_folder_Preprocessing_folder)) {dir.create(Results_folder_Preprocessing_folder)}  # check and create folder
# Create Outlier_Detection directory
Results_folder_Preprocessing_Outlier_detection_folder = file.path(Results_folder_Preprocessing_folder, "Outlier_detection")
if (!dir.exists(Results_folder_Preprocessing_Outlier_detection_folder)) {dir.create(Results_folder_Preprocessing_Outlier_detection_folder)}
# Create Quality_Control_PCA directory
if (ExportQCPlots ==  TRUE){   # Create Quality_Control_PCA directory
   Results_folder_Preprocessing_folder_Quality_Control_PCA_folder = file.path(Results_folder_Preprocessing_folder, "Quality_Control_PCA")
   if (!dir.exists(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder)) {dir.create(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder)}
}

#########################################
### ### ### Feature filtering ### ### ###

message("Feature filtering is performed to reduce missing values that can bias the analysis and cause methods to underperform, which leads to low precision in the statistical analysis. REF: Steuer et. al. (2007), Methods Mol Biol. 358:105-26., doi:10.1007/978-1-59745-244-1_7.")

if (Feature_Filtering == "Modified"){
	message("Here we apply the modified 80%-filtering rule that takes the class information (Column `Conditions`) into account, which additionally reduces the effect of missing values. REF: Yang et. al., (2015), doi: 10.3389/fmolb.2015.00004)")
		message(paste("filtering value selected:", Feature_Filt_Value))
	Input_data <- as.data.frame(Input_data)
	unique_conditions <- unique(Conditions) # saves the different conditions
	if(is.null(unique_conditions) == TRUE){
		stop("Condition information is missing from the Experimental design.")
	}
    if(length(unique_conditions) == 1){
      stop("To perform the Modified feature filtering there have to be at least 2 different Conditions in the `Condition` column in the Experimental design. Consider using the Standard feature filtering option.")
	}   

	miss <- c()
	message("***Performing modified feature filtering***")
   # for (i in unique_conditions){
	split_Input <- split(Input_data, Conditions) # splits data frame into a list of dataframes by condition
	# Select metabolites to be filtered for different conditions
	for (m in split_Input){
		for(i in 1:ncol(m)) {
          if(length(which(is.na(m[,i]))) > (1-Feature_Filt_Value)*nrow(m)) ## Check complete.case instead of is.na. it is faster and you dont have to use which
		miss <- append(miss,i)
			}
		}
    #}
  
    if(length(miss) ==  0){ #remove metabolites if any are found
	message("There where no metabolites exluded")
	filtered_matrix <- Input_data
	# save filtering result
	feat_file_res <- "There where no metabolites exluded"
      write.table(feat_file_res,row.names = FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
	} else {
      message(paste(length(unique(miss)) ," metabolites where removed."))
	message(unique(colnames(Input_data)[miss]))
	filtered_matrix <- Input_data[,-miss]
	# save filtering output
      write.table(unique(colnames(Input_data)[miss]),row.names = FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep =  ""))
	}
}
if (Feature_Filtering== "Standard"){
	message("Here we apply the so-called 80%-filtering rule is used, which removes metabolites with missing values in more than 80% of samples. REF: Smilde et. al. (2005), Anal. Chem. 77, 6729–6736., doi:10.1021/ac051080y")
	Feature_Filt_Value <- as.numeric(Feature_Filt_Value)
	# Check if filtering value input is correct
		if(Feature_Filt_Value > 0 & Feature_Filt_Value < 1){
		message(paste("filtering value selected:", Feature_Filt_Value))
		} else{
		stop("Check input. The filtering value should be between 0 and 1")
		}
	split_Input<- Input_data
	# Select metabolites to be filtered for one condition
	miss <- c()
	message("***Performing standard feature filtering***")
	for(i in 1:ncol(split_Input)) {
		if(length(which(is.na(split_Input[,i]))) > (Feature_Filt_Value)*nrow(split_Input))
		miss <- append(miss,i)
	}
	#remove metabolites if any are found
	if(length(miss) ==0){
		message("There where no metabolites exluded")
		filtered_matrix <- Input_data
		feat_file_res <- "There where no metabolites exluded"
		write.table(feat_file_res,row.names = FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep = ""))
	} else {
		message(paste( length(unique(miss)) ,"metabolites where removed:"))
		message(unique(colnames(Input_data)[miss]))
		filtered_matrix <- Input_data[,-miss]
		# save filtering output
		write.table(unique(colnames(Input_data)[miss]),row.names = FALSE, file = paste(Results_folder_Preprocessing_folder,"/Filtered_metabolites","_",Feature_Filt_Value,"%_",Feature_Filtering,".csv",sep = ""))
	}
}
  if (Feature_Filtering == "None"){
    warning("No feature filtering is selected.")
    filtered_matrix <- as.data.frame(Input_data)
}

filtered_matrix <- as.matrix(mutate_all(as.data.frame(filtered_matrix), function(x) as.numeric(as.character(x))))
################################################
### ### ###Zero variance metabolites ### ### ###
filtered_matrix <- filtered_matrix %>% as.data.frame()
zero_var_metabs <- filtered_matrix[, sapply(filtered_matrix, var,na.rm=TRUE) == 0] %>% colnames()
zero_var_metabs <-paste(zero_var_metabs, collapse = ', ')
if(length(zero_var_metabs)>1){
	message(paste("There are metabolits with Zero variance. These are: ",zero_var_metabs ))
}
# Remove zero variance metabolites
zero_var_removed <- filtered_matrix[, sapply(filtered_matrix, var,na.rm=TRUE) != 0]
filtered_matrix <- zero_var_removed
################################################
### ### ### Missing value Imputation ### ### ###
message("Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `half minimum value`. REF: Wei et. al., (2018), Reports, 8, 663, doi:https://doi.org/10.1038/s41598-017-19120-0")
NA_removed_matrix <- replace(filtered_matrix, filtered_matrix %in% NA, ((min(filtered_matrix, na.rm = TRUE))/2))
# Check if a metabolite has very small changes
#which(apply(NA_removed_matrix,2,sd)<1)

#######################################################
### ### ### Normalisation ### ### ###
if (Normalization == "TIC"){ #Normalise for total ion counts
	message("Total Ion Count (TIC) normalization is used to reduce the variation from non-biological sources, while maintaining the biological variation. REF: Wulff et. al., (2018), Advances in Bioscience and Biotechnology, 9, 339-351, doi:https://doi.org/10.4236/abb.2018.98022")
	RowSums <- rowSums(NA_removed_matrix)
	data_raw_no_zero <- NA_removed_matrix
	Median_RowSums <- median(RowSums) #This will built the median
	Data_TIC_Pre <- apply(NA_removed_matrix, 2, function(i) i/RowSums) #This is dividing the ion intensity by the total ion count
	Data_TIC <- Data_TIC_Pre*Median_RowSums #Multiplies with the median metabolite intensity
	Data_TIC<- as.data.frame(Data_TIC)
	}else if (Normalization == "PQN"){ #Normalise using the PQN method
	source("mswsd_resamp_publi.R")
	  #this is a trick to avoid the Error in select(., any_of(vars)) : unused argument (any_of(vars))
	#select <- dplyr::select
	#select_if <- dplyr::select_if
	NA_removed_matrix<- select_if(as.data.frame(NA_removed_matrix), is.numeric) %>%
								select(-any_of(Standards))
	temp <- NA_removed_matrix
	unbal_reg(t(temp))
	data_raw_no_zero <- temp
	data_raw_no_zero[data_raw_no_zero==0] <- 1
	resamp_mswsd(t(data_raw_no_zero))
	selected_percentage = PQN_selected_persentage
	data_norm_PQN <- norm_unbal(t(data_raw_no_zero),selected_percentage,"PQN") #selected percentages is set to 100% by default. change it in the function call if nescessary
	Data_TIC <- as.data.frame(t(data_norm_PQN))
	} else {# No normalization is applyed
	Data_TIC <- as.data.frame(NA_removed_matrix)
	warning("***NO normalization has been performed***")
	message("Total Ion Count (TIC) normalization is used to reduce the variation from non-biological sources, while maintaining the biological variation. REF: Wulff et. al., (2018), Advances in Bioscience and Biotechnology, 9, 339-351, doi:https://doi.org/10.4236/abb.2018.98022")
}

# Process CoRe
if (CoRe == TRUE){
    
	blankMeans <- colMeans( Data_TIC[grep("blank", Conditions),])
	blankSd <- as.data.frame( apply(Data_TIC[grep("blank", Conditions),], 2, sd) )
	blank_df <- as.data.frame(data.frame(blankMeans, blankSd))
	names(blank_df) <- c("blankMeans", "blankSd")
	blank_df$CV <- blank_df$blankSd / blank_df$blankMeans
	# CV is the sd/mean and it sa measure of variability. above 1 is high and below is ok.
	CV <- (sum(blank_df[blank_df$blankMeans != 0,]$CV > 1)/length(blank_df[blank_df$blankMeans != 0,]$CV))*100
	message(paste0(CV, " of variables have very high variability in the blank samples"))
	#Subtract from each sample the blank mean
	Data_TIC_CoRe <- as.data.frame(t( apply(t(Data_TIC),2, function(i) i-blank_df$blankMeans)))
	message("CoRe data are normalised using CoRe_norm_factor")
	Data_TIC <- apply(Data_TIC_CoRe, 2, function(i) i*CoRe_norm_factor)
	if (var(CoRe_norm_factor) == 0){
		warning("The growth rate or growth factor for normalising the CoRe result, is the same for all samples") 
	}
}

  data_norm <- Data_TIC %>% as.data.frame()


#####################################################
### ### ### Sample outlier identification ### ### ###

message("Identification of outlier samples is performed using Hotellin's T2 test to define sample outliers in a mathematical way (Confidence = 0.99 ~ p.val < 0.01) REF: Hotelling, H. (1931), Annals of Mathematical Statistics. 2 (3), 360–378, doi:https://doi.org/10.1214/aoms/1177732979.")
message(paste("HotellinsConfidence value selected:", HotellinsConfidence))

Outlier_filtering_loop = OutlierLoop
sample_outliers <- list()
scree_plot_list <- list()
outlier_plot_list <- list()
metabolite_zero_var_total_list <- list() #Theo#this is new
zero_var_metab_warning = FALSE
  k =  1
  a =  1
  for (loop in 1:Outlier_filtering_loop){   
    #################################################
    ### ### ### Zero variance metabolites ### ### ###
    metabolite_var <- as.data.frame( apply(data_norm, 2, var) %>% t()) # calculate each metabolites variance
    metabolite_zero_var_list <- list( colnames(metabolite_var)[which(metabolite_var[1,]==0)]) # takes the names of metabollites with zero variance and puts them in list

    if(sum(metabolite_var[1,]==0)==0){
      metabolite_zero_var_total_list[loop] <- 0
    } else if(sum(metabolite_var[1,]==0)>0){
    metabolite_zero_var_total_list[loop] <- metabolite_zero_var_list
    zero_var_metab_warning = TRUE # This is used later to print and save the zero variance metabolites if any are found.
    }

    for (metab in metabolite_zero_var_list){  # Remove the metabolites with zero variance from the data to do PCA
      data_norm <- data_norm %>% select(-all_of(metab))
    }

	### ### PCA ### ###
	PCA.res <- prcomp(data_norm, center = TRUE, scale. = TRUE) # Do PCA
	outlier_PCA_data <- data_norm
	outlier_PCA_data$Conditions <- Conditions
	dtp <- data.frame('Conditions' = outlier_PCA_data$Conditions, PCA.res$x[,c(1,2)])
	pca_outlier <- ggplot(data = dtp) +
				geom_point(aes(x = PC1, y = PC2, colour = Conditions), size = 4, alpha = 0.8, show.legend = TRUE) +
      ggtitle(paste("PCA outlier test filtering round ",loop))+
				theme_classic()+
				geom_hline(yintercept = 0, colour = "black", linewidth = 0.1)+
				geom_vline(xintercept = 0, colour = "black", linewidth = 0.1)+
				geom_text(aes(x = PC1, y = PC2, label = rownames(outlier_PCA_data)),hjust = 0.3, vjust = -0.5,size = 3,alpha = 0.6 )+
	scale_x_continuous(paste("PC1 ",summary(PCA.res)$importance[2,][[1]]*100,"%")) +
	scale_y_continuous(paste("PC2 ",summary(PCA.res)$importance[2,][[2]]*100,"%"))
	plot(pca_outlier)
	outlier_plot_list[[k]] <- recordPlot()
	dev.off()
	k = k+1
	### ### Scree plot ### ### what is the assumption for the knee (% of variance expained)
	# get Scree plot values for inflection point calculation
	inflect_df <- as.data.frame(c(1:length(PCA.res$sdev)))
	colnames(inflect_df) <- "x"
	inflect_df$y <- summary(PCA.res)$importance[2,]
	inflect_df$Cumulative <- summary(PCA.res)$importance[3,]
	#make cumulative variation labels for plot
	screeplot_cumul <- format(round(inflect_df$Cumulative[1:20]*100, 1), nsmall = 1)
	# Calculate the knee and select optimal number of components
	knee = inflection::uik(inflect_df$x,inflect_df$y)
	if (is.null(npcs)){
	npcs = knee -1 #Note: we subtract 1 components from the knee cause the root of the knee is the PC that does not add something. npcs=30
	}
###############define the knee based on PAC plots###################
	# Make a scree plot with the selected component cut-off for HotellingT2 test
	screeplot <- factoextra::fviz_screeplot(PCA.res, main = paste("PCA Explained variance plot filtering round ",loop, sep = ""),
										addlabels = TRUE,
										ncp = 20,
										geom = c("bar", "line"),
										barfill = "grey",
										barcolor = "grey",
										linecolor = "black",linetype = 1) + theme_classic()+ geom_vline(xintercept = npcs+0.5, linetype = 2, colour = "red") +
										annotate("text", x = c(1:20),y = -0.8,label = screeplot_cumul,col = "black", size = 3)
	plot(screeplot)
	# save plot
	outlier_plot_list[[k]] <- recordPlot()
	dev.off()
	k = k+1
	### ### HotellingT2 test for outliers ### ###
	data_hot <- as.matrix(PCA.res$x[,1:npcs])
	message("***Checking for outliers***")
    hotelling_qcc <- qcc::mqcc(data_hot, type = "T2.single",labels = rownames(data_hot),confidence.level = HotellinsConfidence, title = paste("Outlier filtering via HotellingT2 test filtering round ",loop,", with ",HotellinsConfidence, "% Confidence",  sep = ""), plot = FALSE)
	### Multivariate quality control chart for Hotteling T2
	HotellingT2plot_data <- as.data.frame(hotelling_qcc$statistics)
	HotellingT2plot_data <- rownames_to_column(HotellingT2plot_data, "Samples")## alternative to cbind
	#HotellingT2plot_data <- cbind(row.names(HotellingT2plot_data),HotellingT2plot_data) # an alternative is to use the rownames_to_column function
	colnames(HotellingT2plot_data) <- c("Samples", "Group summary statisctics")
	outlier <- HotellingT2plot_data %>% filter(HotellingT2plot_data$`Group summary statisctics`>hotelling_qcc$limits[2])
	limits <- as.data.frame(hotelling_qcc$limits)
	legend <- colnames(HotellingT2plot_data[2])
	LegendTitle = "Limits"

	HotellingT2plot <- ggplot(HotellingT2plot_data, aes(x = Samples, y = `Group summary statisctics`, group = 1, fill = ))
	HotellingT2plot <- HotellingT2plot + 
					geom_point(aes(x = Samples,y = `Group summary statisctics`), colour = 'blue', size = 2) + 
					geom_point(data = outlier, aes(x = Samples,y = `Group summary statisctics`), colour = 'red',size = 3) + 
					geom_line(linetype = 2)

	#draw the horizontal lines corresponding to the LCL,UCL
	HotellingT2plot <- HotellingT2plot + geom_hline(aes(yintercept = limits[,1]), colour = "black", data = limits,  show.legend = F) + 
	geom_hline(aes(yintercept = limits[,2], linetype = "UCL"), colour = "red", data = limits, show.legend = T) +
	#only the LCl and UCL to be shown in y axis
	scale_y_continuous(breaks = sort(c(ggplot_build(HotellingT2plot)$layout$panel_ranges[[1]]$y.major_source, c(limits[,1],limits[,2]))))

    HotellingT2plot <- HotellingT2plot + theme_classic()
    HotellingT2plot <- HotellingT2plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    HotellingT2plot <- HotellingT2plot + ggtitle(paste("Outlier filtering via Hotelling ", hotelling_qcc$type ," test filtering round ",loop,", with ", 100 * hotelling_qcc$confidence.level,"% Confidence"))
    #plot1 <- plot1 + scale_linetype_manual(name = LegendTitle,values = "dashed") # this line instead of the next for bashed red line
    HotellingT2plot <- HotellingT2plot + scale_linetype_discrete(name = LegendTitle,)
    HotellingT2plot <- HotellingT2plot + theme(plot.title = element_text(size = 10, face = "bold")) +
      theme(axis.text = element_text(size = 12))

    plot(HotellingT2plot)
    outlier_plot_list[[k]] <- recordPlot()
    dev.off()
    k = k+1

    ### Save the outlier detection plots in the outlier detection folder
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "/PCA_OD_round_" ,a ,".", Save_as, sep = ""),
           plot = pca_outlier, width = 10,height = 8)
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "//Scree_plot_OD_round_" ,a ,".",Save_as, sep = ""),
           plot = screeplot, width = 10,height = 8)
    ggsave(filename = paste(Results_folder_Preprocessing_Outlier_detection_folder, "/Hotelling_OD_round_" ,a ,".", Save_as, sep = ""),
           plot = HotellingT2plot, width = 10,height = 8)
    a = a+1

    # Here for the outliers we use confidence of 0.999 and p.val < 0.01. We could use the potential outliers for pval 0.05 and confidence 95
    if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 0){ # loop for outliers until no outlier is detected
      data_norm <- data_norm  # filter the selected outliers from the data
      break
    } else if (length(hotelling_qcc[["violations"]][["beyond.limits"]]) == 1){
      data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]
      Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]

      if(is.null(Biological_Replicates)!=TRUE){
        Biological_Replicates <- Biological_Replicates[-hotelling_qcc[["violations"]][["beyond.limits"]]]
      }
      if(is.null(Analytical_Replicates)!=TRUE){
        Analytical_Replicates <- Analytical_Replicates[-hotelling_qcc[["violations"]][["beyond.limits"]]]
      }

      # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
      hotelling_qcc[["violations"]][["beyond.limits"]][1] <-   rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][1]]
      sample_outliers[loop] <- list(hotelling_qcc[["violations"]][["beyond.limits"]])

    } else {
      data_norm <- data_norm[-hotelling_qcc[["violations"]][["beyond.limits"]],]
      Conditions <- Conditions[-hotelling_qcc[["violations"]][["beyond.limits"]]]
      if(is.null(Biological_Replicates)!=TRUE){
        Biological_Replicates <- Biological_Replicates[-hotelling_qcc[["violations"]][["beyond.limits"]]]
      }
      if(is.null(Analytical_Replicates)!=TRUE){
        Analytical_Replicates <- Analytical_Replicates[-hotelling_qcc[["violations"]][["beyond.limits"]]]
      }

      # Change the names of outliers in mqcc . Instead of saving the order number it saves the name
      sm_out <- c() # list of outliers samples
      for (i in 1:length(hotelling_qcc[["violations"]][["beyond.limits"]])){
        sm_out <-  append(sm_out, rownames(data_hot)[hotelling_qcc[["violations"]][["beyond.limits"]][i]]  )
      }
      sample_outliers[loop] <- list(sm_out )
    }
  }

  # Save outlier result
 # pdf(file= paste("Results_", Sys.Date(),"/","Preprocessing", "/Outlier_testing.pdf", sep = ""), onefile = TRUE ) # or multivariate quality control chart
  Results_folder_Preprocessing_Outlier_detection_folder_out<- file.path(Results_folder_Preprocessing_Outlier_detection_folder, "Outlier_testing.pdf")
  pdf(file = Results_folder_Preprocessing_Outlier_detection_folder_out, onefile = TRUE )
  for (plot in outlier_plot_list) {
    replayPlot(plot)
  }
  dev.off()

  
  # Print Outlier detection results about samples and metabolites
  if(length(sample_outliers) > 0){   # Print outlier samples
    message("There are possible outlier samples in the data") #This was a warning
    for (i in 1:length(sample_outliers)  ){
      message("Filtering round ",i ," Outlier Samples:", paste( head(sample_outliers[[i]]) ," "))
    }
  }else{  message("No sample outliers were found")}

  ### Potential Outlier testing.
  # When we are done with the filtering of the outliers. We do one more round with lower thresholds for potential outliers.
  PCA.res<- prcomp(data_norm, center = TRUE, scale. = TRUE) # Do PCA
  inflect_df<- as.data.frame(c(1:length(PCA.res$sdev)))
  colnames(inflect_df)<- "x"
  inflect_df$y<- summary(PCA.res)$importance[2,]
  knee=inflection::uik(inflect_df$x,inflect_df$y)
  npcs = knee -1
  data_hot<- as.matrix(PCA.res$x[,1:npcs])

  mqcc95 <- qcc::mqcc(data_hot, type="T2.single",labels=rownames(data_hot), confidence.level = 0.95, plot = FALSE) # these appear as outliers in the 0.95 confidence of the mqcc plot
  potent_Outlier<- as.vector(mqcc95[["violations"]][["beyond.limits"]])
  if(length(potent_Outlier)==1){
    n_potent_outlier <-  rownames(data_hot)[mqcc95[["violations"]][["beyond.limits"]]]
    warning(paste(length(potent_Outlier), "Potential outlier found. Sample", n_potent_outlier))
  }else if(length(potent_Outlier)>1){
    n_potent_outlier <-  rownames(data_hot)[mqcc95[["violations"]][["beyond.limits"]]]
    paste(c(length(potent_Outlier), " Potential outliers found. Samples", n_potent_outlier), collapse=" ")
  } else{  message("No potential outliers were found")}

  if(length(potent_Outlier)>0){
    # The user can filter the potentially outlier samples
    if (Filter_Potential_Outliers== TRUE){
      message("*** Filtering potential outliers ***")
      data_norm<- data_norm[-potent_Outlier,]
      Conditions <- Conditions[-potent_Outlier]
      if(is.null(Biological_Replicates)!=TRUE){
        Biological_Replicates <- Biological_Replicates[-potent_Outlier]
      }
      if(is.null(Analytical_Replicates)!=TRUE){
        Analytical_Replicates <- Analytical_Replicates[-potent_Outlier]
      }
    }else{
      message("*** No filtering for potential sample outliers performed ***")
    }
  }
  data_norm_filtered <- data_norm
  data_norm_filtered <- as.data.frame(as.matrix(mutate_all(as.data.frame(data_norm_filtered), function(x) as.numeric(as.character(x)))))

  # IF in the output we use the filtered exp design the this line is ok . If not then remove this line
  Experimental_design_filtered<- Experimental_design%>%
    filter( rownames(Experimental_design) %in% rownames(data_norm_filtered) )

  ### make outlier columns and add them to output dataframe
  total_outliers<- hash::hash() # make a dictionary
  # Create columns with outliers to merge to output dataframe
  if(length(sample_outliers)>0){
    for (i in 1:length(sample_outliers)  ){
      total_outliers[[paste("Outlier_filtering_round_",i, sep = "")]] <- sample_outliers[i]
    }
  }

  if(length(potent_Outlier)>0){
    total_outliers[["Potential_Outliers"]] <- list(n_potent_outlier)
  }

  data_norm_filtered_full<- as.data.frame(Data_TIC)
  
    
  # add outlier information to the full output dataframe
  if(length(total_outliers)>0){
    data_norm_filtered_full$Outliers <- "no"
    for (i in 1:length(total_outliers)){
      for (k in 1:length( hash::values(total_outliers)[i] ) ){
        data_norm_filtered_full[as.character(hash::values(total_outliers)[[i]]) , "Outliers"] <- hash::keys(total_outliers)[i]
      }
    }
  }else{
    data_norm_filtered_full$Outliers <- "no"
  }
  
  
  #Put Outlier columns in the front
  data_norm_filtered_full <- data_norm_filtered_full %>% relocate(Outliers)


  # add the design in the output df (merge by rownames/sample names)
  data_norm_filtered_full <- merge(Experimental_design, data_norm_filtered_full,  by = 0)
  rownames(data_norm_filtered_full) <- data_norm_filtered_full$Row.names
  data_norm_filtered_full$Row.names <- c()

  ################################################
  
  ### ### ### Quality Control (QC) PCA ### ### ###

  QC_PCA_data<- data_norm_filtered
  QC_PCA_data$Conditions <- Conditions
  if(is.null(Biological_Replicates)!=TRUE){
    QC_PCA_data$Biological_Replicates <-  as.character(Biological_Replicates)
  }
  pca.obj <- prcomp(data_norm_filtered, center = TRUE, scale. = TRUE)
  dtp <- data.frame('Conditions' = QC_PCA_data$Conditions, pca.obj$x[,1:2])
  
  ## PCA conditions
  pca_QC <- ggplot(data = dtp) +
    geom_point(aes(x = PC1, y = PC2, colour = Conditions), size = 4, alpha=0.8) +  
    ggtitle("Quality COntrol PCA COndition clustering check")+
    theme_classic()+
    geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
    geom_vline(xintercept=0,  color = "black", linewidth=0.1)+
    geom_text(aes(x = PC1, y = PC2, label = rownames(QC_PCA_data)),hjust=0.3, vjust=-0.5,size=3,alpha=0.6 )+
    scale_x_continuous(paste("PC1 ",summary(pca.obj)$importance[2,][[1]]*100,"%")) +
    scale_y_continuous(paste("PC2 ",summary(pca.obj)$importance[2,][[2]]*100,"%"))
  #create a qUALITY cONTROL OF THE CORRECTED MEANS
  mdatRAW <- reshape2::melt(t(log(as.matrix(data_raw_no_zero),2)))  ## convert to long format
  mdat <- reshape2::melt(t(log(as.matrix(data_norm_filtered),2)))  ## convert to long format
  YlinE <- log2(median(Biobase::rowMedians(as.matrix(data_norm))))
  YlinERAW <- log2(median(Biobase::rowMedians(as.matrix(data_norm))))
  if (Normalization == "PQN"){
    NORM = paste( "with selected %", PQN_selected_persentage, sep = " ")
  }else{
    NORM = "" 
  }
  plotQC_RAW <- ggplot2::ggplot(mdatRAW,aes(x=factor(Var2),y=value, fill = Var2))+
    ggplot2::geom_violin(size = 0.1)+
    geom_boxplot(aes(x = factor(Var2), y = value),outlier.colour = "darkred",
                 outlier.size = 0.9,width = 0.2, size = 0.3)+
    ggplot2::guides(fill=guide_legend(title="samples")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = YlinERAW), size = 0.2) +
    theme_classic()+
    theme(axis.text.x = element_text(size = rel(0.9), angle=45, hjust = 1, vjust = 1)) +
    ggtitle(paste("QC of raw values", sep = " "))+
    theme(plot.title = element_text(size = 10, face = "bold")) +
    xlab("") +
    ylab("log2 of ")
  
  plotQC_NORM <- ggplot2::ggplot(mdat,aes(x=factor(Var2),y=value, fill = Var2))+
                        ggplot2::geom_violin(size = 0.1)+
                        geom_boxplot(aes(x = factor(Var2), y = value),outlier.colour = "darkred",
                                     outlier.size = 0.9,width = 0.2, size = 0.3)+
                        ggplot2::guides(fill=guide_legend(title="samples")) +
                        ggplot2::geom_hline(ggplot2::aes(yintercept = YlinE), size = 0.2) +
                        theme_classic()+
                        theme(axis.text.x = element_text(size = rel(0.9), angle=45, hjust = 1, vjust = 1)) +
                        ggtitle(paste("QC of Normaliation Using", Normalization ,NORM, sep = " "))+
                        theme(plot.title = element_text(size = 10, face = "bold")) +
                        xlab("") +
                        ylab("log2 of ")
  
    #geom_jitter(aes(x = reorder(factor(Var2), desc(value)), y = value)) +
    arrange <- ggpubr::ggarrange(plotQC_RAW,  plotQC_NORM, ncol = 1, nrow = 2) 
    WIDTH = (length(unique(mdat$Var2))+4)*0.6
    HEIGTH = max(mdat$value)*0.5
  ggsave(filename = paste0(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder, "/PCA_Condition_Clustering.",Save_as), plot=pca_QC, width =10,  height = 8)
  ggsave(filename = paste0(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder, "/Median_correction.",Save_as), plot=arrange,  device ="svg")
  
  
  if(is.null(Biological_Replicates)!=TRUE){
    ### ### QC PCA color for replicates
    pca_QC_repl <- ggplot(data = dtp) +
      geom_point(aes(x = PC1, y = PC2, colour = Conditions, shape =as.factor(Biological_Replicates)), size = 4, alpha=0.8) +  
      ggtitle("Quality Control PCA replicate spread check")+
      theme_classic()+
      geom_hline(yintercept=0,  color = "black", linewidth=0.1)+
      geom_vline(xintercept=0,  color = "black", linewidth=0.1)+
      geom_text(aes(x = PC1, y = PC2, label = rownames(QC_PCA_data)),hjust=0.3, vjust=-0.5,size=3,alpha=0.6 )+
      scale_x_continuous(paste("PC1 ",summary(pca.obj)$importance[2,][[1]]*100,"%")) +
      scale_y_continuous(paste("PC2 ",summary(pca.obj)$importance[2,][[2]]*100,"%"))
    message()
    ggsave(filename = paste0(Results_folder_Preprocessing_folder_Quality_Control_PCA_folder, "/PCA_replicate_distribution.",Save_as), plot=pca_QC_repl, width =10,  height = 8)
    
  }

  ###################################################################
  ### ### ### Merge analytical replicates (if they exist) ### ### ###

  if(is.null(Analytical_Replicates)!=TRUE){
    data_norm_filtered_sumed<- merge(Experimental_design%>% select(Biological_Replicates,Analytical_Replicates, Conditions),data_norm_filtered, by=0)
    rownames(data_norm_filtered_sumed) <- data_norm_filtered_sumed$Row.names
    data_norm_filtered_sumed$Row.names <- c()

    # Problem here what names and information from design to keep in the data frame
    data_norm_filtered_sumed<-as.data.frame( data_norm_filtered_sumed %>%
                                               group_by(Biological_Replicates, Conditions) %>%
                                               summarise_all("mean"))
    
    data_norm_unfiltered_sumed<- merge(Experimental_design%>% select(Biological_Replicates,Analytical_Replicates, Conditions),Data_TIC, by=0)
    rownames(data_norm_unfiltered_sumed) <- data_norm_unfiltered_sumed$Row.names
    data_norm_unfiltered_sumed$Row.names <- c()
    data_norm_unfiltered_sumed<-as.data.frame(data_norm_unfiltered_sumed %>%
                                                group_by(Biological_Replicates, Conditions) %>%
                                                summarise_all("mean"))
    
  }

  
  ####################################
  ### ### ###  Make output ### ### ###

  output_list <- list()  #Here we make a list in which we will save the output

  # In the next line we use the experimental_design. there is also the experimental_design_filtered which has the same sample with the data_processed_filtered
  # return dictionary with the above
  if(is.null(Analytical_Replicates)!=TRUE){
    # merge by row names (by=0 or by="row.names")
    preprocessing_output_list <- list(Experimental_design=Experimental_design, Raw_data=as.data.frame(Input_data), Processed_data= data_norm_filtered_full,
                                      Processed_data_matrix_OutliersRemoved =data_norm_filtered, Processed_summed_data =data_norm_filtered_sumed,
                                      Processed_unfiltered_summed_data = data_norm_unfiltered_sumed)
  }else{
    preprocessing_output_list <- list(Experimental_design=Experimental_design, Raw_data=as.data.frame(Input_data), Processed_data= data_norm_filtered_full,
                                      Processed_data_matrix_OutliersRemoved =data_norm_filtered)
  }


  ##Write to file
  preprocessing_output_list_out <- lapply(preprocessing_output_list, function(x) rownames_to_column(x, "Sample_ID")) #  # use this line to make a sample_ID column in each dataframe
  writexl::write_xlsx(preprocessing_output_list_out, paste(Results_folder_Preprocessing_folder, "/Preprocessing_output.xlsx", sep = ""))#,showNA=TRUE)

  # Return the result
  message("Done")
  return(preprocessing_output_list )
}




