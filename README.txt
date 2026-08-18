# Metabolomic analysis of ischemia in mouse and naked mole-rat

This repository contains the R code used for the metabolomics analyses presented in *Liver-like glycogen metabolism supports glycolysis in naked mole-rat heart during ischemia*.

The analyses include LC-MS metabolomics data from mouse and naked mole-rat tissues following ischemia, as well as neonatal mouse samples.

## Analyses

The repository contains the code used for the following analyses:

1. **Heart ischemia time course – mouse and naked mole-rat (Figure 1)**  
   Metabolomic analysis of heart tissue collected at different time points following ischemia.

2. **Ischemic risk and adjacent zones (Figure 1)**  
   Metabolomic analysis comparing ischemic tissue regions.

3. **Neonatal mouse comparison (Figure 3)**  
   Analysis of neonatal mouse samples and comparison with the previously acquired mouse and naked mole-rat metabolomics data.

## LC-MS data processing

Raw metabolite intensity tables were generated using TraceFinder or Compound Discoverer and processed in R using MetaProViz.

For the analyses corresponding to Figures 1B–1H, feature filtering was performed using the 80% filtering rule per condition, followed by half-minimum missing value imputation, total ion count (TIC) normalisation and outlier detection based on Hotelling's T2.

For Figures 1K–1Q, the same processing workflow was applied, except that metabolite intensities were obtained using Compound Discoverer and probabilistic quotient normalisation (PQN) was used instead of TIC normalisation.

For Figures 3I–J, neonatal mouse metabolite intensities generated using TraceFinder were batch corrected against the previously acquired dataset using internal-standard median correction. Each saccharide intensity was divided by the median valine-d8 intensity of its respective analytical batch.

## Reproducibility

The R scripts required to reproduce the analyses can be found in the corresponding folders:

- `Figure1 Ischemia_Heart_MousevsNMR_timecourse/Ischemia_mouse_NMR.R`
- `Figure1 Ischemia_risk_and_adj_zone/2024_Maria_MSJ01/frezzaMetaboAnalysis/FeatureSelection/Analysis.Rmd`
- `Figure3 Neonatal_mouse_mouse_NMR/Analysis.R`

The scripts include data preprocessing, normalisation, PCA, clustering and generation of the corresponding visualisations.

PCA was performed using `prcomp`. K-means clustering and selection of the number of clusters using the silhouette method were performed using the `factoextra` and `cluster` packages. Heatmaps were generated using `pheatmap`.

Additional data manipulation and visualisation were performed using `tidyverse`, `reshape2` and `readxl`.

## Data

The metabolomics data associated with this study are deposited on Zenodo under DOI: `10.5281/zenodo.21979323`.

The input tables required to reproduce the R analyses are provided in the corresponding analysis folders.

LC-MS data acquisition was performed by the Frezza Laboratory, CECAD Research Center, Faculty of Medicine and University Hospital Cologne, Germany.

## Software

The analyses were performed in R. Major packages include:

- MetaProViz
- tidyverse
- factoextra
- cluster
- pheatmap
- reshape2
- readxl
