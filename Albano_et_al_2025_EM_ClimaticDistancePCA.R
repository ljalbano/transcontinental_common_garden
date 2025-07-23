# Data analysis for Albano et al., submitted to Ecological Monographs
# Title: Adaptation to environmental variation in the native and introduced ranges of a cosmopolitan plant

# This script contains the necessary information to recreate the principal component analysis used to calculate weighted climatic distance

# Set the working directory
setwd("")

# Load all required packages
library(ggplot2) # for plotting principal component analyses
library(factoextra) # for generating and plotting principal component analyses
library(ggpubr) # for arranging multi-panel figures
library(writexl) # for writing excel output files of PCA results

# Prior to this PCA construction:

# Extract data for 6 bioclimatic variables of interest for each T. repens population and common garden location (from WorldClim 2.1)
    # 1 = Mean Annual Temperature
    # 5 = Max Temperature of Warmest Month
    # 6 = Min Temperature of Coldest Month
    # 12 = Annual Precipitation
    # 13 = Precipitation of Wettest Month
    # 14 = Precipitation of Driest Month
# Combine into a single dataframe and z-score standardize all values for PCA construction

#***************************** 
# PCA CONSTRUCTION AND EXPORT
#*****************************

# Load and attach z-score standardized bioclimatic variables
# CHECK FILENAMES
Bioclim = read.csv("Albano_et_al_2025_EM_WorldClim_z-standardized.csv", header = TRUE, sep = ",")
attach(Bioclim)

# Run principal component analysis on CG_Bio dataset (z-score adjusted)
climatic_distance_pca <- princomp(Bioclim)

# Visualize the PC axes and datapoints using the following plots:

# PCA-Biplot containing all datapoints and vectors for each BIO variable
fviz_pca(climatic_distance_pca)
# Some collinearity within temperature-based and precipitation-based variables

# Scree plot of 6 dimensions (used to select dimensions 1, 2, and 3 based on the amount of variance they explain)
fviz_eig(climatic_distance_pca)
# Can then extract the actual eigenvalues and variance percentages for each dimension
get_eigenvalue(climatic_distance_pca)

# Separately view the datapoints and vectors on PCA plots (easier viewing)
fviz_pca_ind(climatic_distance_pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
fviz_pca_var(climatic_distance_pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# create three separate biplots for the three combinations of three PC axes chosen to remain in the PCA to view the contributions of each
biplot.1.2<-fviz_pca_var(climatic_distance_pca, axes = c(1,2), repel = TRUE, label = "all", title = "")
biplot.1.3<-fviz_pca_var(climatic_distance_pca, axes = c(1,3), repel = TRUE, label = "all", title = "")
biplot.2.3<-fviz_pca_var(climatic_distance_pca, axes = c(2,3), repel = TRUE, label = "all", title = "")

# Combine biplots into Figure S1
figureS1<-ggarrange(biplot.1.2, biplot.1.3, biplot.2.3,
                   align = "hv",
                   ncol = 3, nrow = 1,
                   labels=c("A","B","C"),
                   font.label=list(face="plain"),
                   widths = c(1, 1))
figureS1

# Extract loadings
PCA_scores_df <- as.data.frame(climatic_distance_pca$scores)  # convert to a plain matrix

PCA_scores <- princomp(PCA_scores_df)
# Visualize the PC axes and datapoints using the following plots and the scores dataframe
fviz_pca(PCA_scores)
fviz_eig(PCA_scores)
fviz_pca_ind(PCA_scores, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

# Get PC coordinates and write in excel file
climatic_distance_pca_ind <- get_pca_ind(climatic_distance_pca)
climatic_distance_pca_coord <- climatic_distance_pca_ind$coord
write.table(climatic_distance_pca_coord, file = "Albano_et_al_2025_EM_Climatic_Distance_PCA_Scores.csv", sep = ",", quote = FALSE, row.names = F)
# This output table will be used to calculate weighted climatic distance for the first 3 PC axes (done outside of R; see Methods)
# That variable will be contained as WeightedCD in the experimental data .csv file
