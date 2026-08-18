# frezzaMetaboAnalysis

Here are all the steps for the Feature Selection pipeline for metabolomics analysis in Frezza lab.


The Starting Script is called 'BiomarkerDiscoveryAnalysis.Rmd' and is located in the frezzaMetaboAnalysis/FeatureSelection/Scripts folder.


#-1 Compound Discoverer
Before any analysis, Compound Discoverer is used to extract the features and raw Intensity values from the raw data. Add how to run CD.
 The initial CD Analysis is saved in the '3. XXXX CD Data/XXXX'  (XXXX is the ProjectCode (ie. DS02)).

Exporting CD result as table:
A) On the resulting table in Compound Discoverer on the top left there is a square (above the feature numbering, left of the column names of the table) which selects the columns shown on the table called 'Chose fields'. Select: Area, Formula, m/z, Name, Peak Rating?, Reference Ion, RT [min]. 
B) Filter Compounds: In the filtering area select: Reference Ion contains ']+' for positive or does not contain ']+' for negative.
C) Export as Excel with name: 'Compounds_pos.xlsx' or 'Compounds_neg.xlsx' inside the ProjectCode folder inside the CD folder inside of the ProjectFolder (ie. '3. XXXX CD Data/XXXX').

These 2 files are the main input for the pipeline.


#0. Load libraries and directories
In the first block of code we set the ProjectFolder and the workingDirectory. The ProjectFolder is the main folder of the project in the Mass spec Data on the network (ie. ProjectFolder <- "Z:/1. Mass Spec Projects/FREZZA/Debarati SANYAL/2024_Debarati_DS02") and the workingDirectory is the main Github folder (ie. setwd("C:/Users/Dimitrios/Documents/GitHub")).


#1. CD to Binner
This chunk takes the CD output and processes it in order to be used in Binner software. We use the PreBinnerAnalysis script.
Inputs:
ProjectFolder = ProjectFolder, is the projectFolder from #0.
FolderName = "preBinnerAnalysis", The output folder name inside the 5. XXXX Analysed Data folder, as this is the analysis done before Binner the name is standardized to preBinnerAnalysis.
InputName = NULL, The standard CD Analysis output should be in the folder '3. XXXX CD Data/XXXX'. In case we have another CD Analysis in the folder (ie. XXXX_test) and would like to use it as input we set the InputName to the additional name (ie. InputName = "test", the underscore is automatically added).


PreBinnerAnalysis
For each file in the '3. XXXX CD Data/XXXX' folder (ie. Compounds_XXXX_pos/neg.xlsx)
take Name, mz, rt, areas. 
Make featureName = mz_RT and message if there are any duplications. 
In pos dataset if we find the feature with min(abs(mz - 126.13645)) and take the RT. If the RT-5.67 <0.2 we set the name of that feature to "Valine-d8".
We remove any blank samples.
We load the SSF and get the metadata (The xlsx files that exist in the '1. RP01 SSFs' folder).
We load the runOrder (The csv files that exist in the '2. RP01 Raw Data' folder).


We check the data for zero variance features and remove them.
Then we run PoolEstimation.
	We take the pool samples and calculate, the CV for each feature and the percentage of NAs. We remove features with a single value and those with CV> 30% as unstable features, and save the names of metabs as lists in csv files.
Then we use the filtered data, remove all features that have at least 1 missing value and do  PCA as PCA doesnt accept NAs. For PCA iF we have 0 metabs with no missing values then an empty plot is created, if we have 1 feature then we do violin plot and if we have more, we do normal pca. We also make a histogram and a violin plot of CVs from the Pool samples. Using the sample order we order the samples as they were run in the mass spec, calculate the median intensity of the sample and make barplot of that across the sample order. With this we visually check the signal drift. We also do the same but with grouping based on the conditions and calculate the group means and add them as horizontal lines in the plot. Finally, we do the signal drift plot but only for valined-8 in the positive dataset.
	After Pool estimation we filter features using the 80% modified rule. We take the data ignoring the pool samples. we remove the features that have more NAs that 20% in each group. Meaning that we keep the features that have more than 80% of values in at least one sample group.
	Then Missing value imputation is performed using the group mean. For each feature we calculate the NAs. IF all values are NA or if there is only 1 value, we set all values to half minimum of the whole dataset. If we have more than 1 value we calculate the group mean and set all NAs equal to the group mean.
	Next is normalisation. First we set all zero values to one (there shouldnt be any). Then we use the PQN normalisation. We plot the sample median intensity of the raw and normalised data as well as PCA.
	Next is Outlier Detection. We check for zero variance features and remove them. We do PCA and get the components. Find the inflection point (knee), select only the wanted PCA components and visualise with a scree plot. Then we do the hotellings test on the selected components and a plot. IF any outliers we found they are removed and the processed is repeated. After outlier removal we check the features for zero variance and do PCA.

Question: In this stage we remove any samples as outliers. Should we or not? 

The above creates 2 files in the Results_preBinnerAnalysis foder. Compounds_XXXX_pos/neg_processed.csv

#2. Binner
Binner is used to remove redundant features from the data. Binner is a solfware outside R. To run it, open it and as input files go to ProjectFolder, then go to the subfolder: 5. XXXX Analysed Data/Results_preBinnerAnalysis and use Binner with the Compounds_XXXX_pos/neg_processed.csv
This will create the Compounds_XXXX_pos/neg_processed_Report.xls in the same folder.

In binner we uncheck the log transformation. (checking it or not does not change results (tested and visually checked) as they are already normalised.)

Now we remove redundant features and combine the 2 datasets to create the clean dataset for analysis using the "PostBinnerDataPreparation.R" script. The script for each mode (pos/neg) automatically searches the '3. XXXX CD Data/XXXX_InputName' for the data and loads them. Then, searches the Results_preBinnerAnalysis folder for xlsx files (the result from Binner), loads them, identifies the redundant features (the marked yellow rows), saves the feature names (mz_RT) as csv and removes them from the data. In pos mode valine is added as a name. Finally the pos and neg modes are merged and the result is saved as Data_clean_XXXX in the Results_porBinnerAnalysis folder.

Here we have to add a function for duplicated features (features found in both modes) and select one.


#3. Data Preprocessing.
Now that we have cleaned the raw data we preprocess them. First we lead the clean data, the SSF and the sample run order as did before. Make the nesesary changes and corrections to the SSF and check for zero variance features. We run Pool Estimation and afterwards Valide-d8 is removed, Feature filtering, Missing Value Impoutation, Normalisation and Hotelling same as before. Additionaly PCAs are created for all combinations of the first 5 components and colored by all columns in SSF one at a time. 


