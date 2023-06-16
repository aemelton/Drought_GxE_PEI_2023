### AE Melton, 2020
# Perform a principal component analysis on data collected using the "GetRasterData.R" script
# Identify clusters of environmental space among geographically isolated populations

#
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(ggfortify)
library(ggpubr)
library(pca3d)
library(forcats)
library(dplyr)
library(ggmap)


# Read in data
df <- read.csv("Outputs/env.csv")

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

avg_sil <- function(k) {
  km.res <- kmeans(df, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(df))
  mean(ss[, 3])
}

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)
pdf("Silhouette.pdf") 
plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")
#
dev.off()
#
fviz_nbclust(df, kmeans, method = "silhouette") # Optimal k is 2
#

#
km.res <- kmeans(df, centers = 2, nstart = 1234) # Using kvalue from previous function

# Dimension reduction using PCA
res.pca <- prcomp(df) 
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(km.res$cluster)
# Add Species groups from the original data sett
ind.coord$Species <- df$Species
# Data inspection
head(ind.coord)
write.csv(ind.coord, "ind_coord.csv", row.names = F)
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)
write.csv(x = eigenvalue, file = "eigen_values.csv", row.names = F)
#