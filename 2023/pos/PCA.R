library(readr)
library(factoextra)
library(tidyverse) #used for data science. The eight core packages inside this library are: ggplot2 (data visualisation), dplyr (data manipulation), tidyr, readr, purrr, tibble, stringr, and forcats
library(KODAMA) # to use the normalisation function
library(ggrepel) #mainly used to repel overlapping text labels in ggplots
library(vegan) #popular library for analysing ecological diversity and for multivariate analysis of community data. Here, we use it for PCoA
library(IRdisplay) #better display of output cells in Jupyter Notebooks running iwth IRKernel. Library not needed when running the script in RStudio
library(svglite) #to save the plots in support vector graphics (svg) format
library(factoextra) #for extracting and visualizing outputs of multivariate analyses such as PCA, k-means
library(ggsci) #provides color palettes for ggplot2 that can be used for scientific journals
library(matrixStats) #contains highly optimized functions to perform statistics on matrix data
library(cowplot) #efficient functions to arrange several plots
library(ComplexHeatmap) #for visualising heatmaps
library(dendextend) # for getting dendograms
library(NbClust) # for finding the optimum no.of clusters to be used for a clustering method
library(tidyverse)


setwd('/home/anastasiia/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Multivariate')
getwd()

metadata<- read_csv("~/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Multivariate/07052023_metadata_pos.csv")
metadata <- metadata[-31:-32,] #exclude blank
scaled_table <- read_csv("~/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Multivariate/CLR_Scaled.csv")
colnames(scaled_table)[1] <- "filename"
data_merge <- metadata %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
   column_to_rownames("filename")

res.pca <- data_merge %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
#fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Diet, addEllipses = TRUE,repel = TRUE, legend.title = "Diet",
             #palette = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"))
fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("green4", "lightgreen","royalblue","skyblue2","red3", "#FAA0A1", "yellow", "orange"),
             geom.ind = "point", pointsize = 2)

#deficient only
metadata_sorted <- metadata [order(metadata $ATTRIBUTE_iron_level), ]
metadata_d <- metadata_sorted [-17:-30,]
data_merge_d <- metadata_d %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d <- data_merge_d %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d)
fviz_pca_ind(res.pca_d, col.ind = data_merge_d$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3","orange", "royalblue","green4"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")
view(pc_loadings_d)
pc_loadings_d[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(pc_loadings_d$'feature', '_', 3)
fviz_pca_var(res.pca_d, select.var = list(contrib = 50), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#deficient only-no dKO
metadata_d2 <- metadata_d [1:13,]
data_merge_d2 <- metadata_d2 %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_d2 %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_d2$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue","green4"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

pc_loadings_d[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(pc_loadings_d$'feature', '_', 3)
pc_loadings_d <- as.matrix(pc_loadings_d[, c(1,20,18,19)])

write.csv(pc_loadings_d, "pc_loadings_d.csv",row.names = TRUE)
fviz_pca_var(res.pca_d, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


####################PERMANOVA############################
distm <- dist(data_merge_d, method = "euclidean")# compute distance

adonres <- adonis2(distm ~ data_merge_d[,colnames(data_merge_d) == 'ATTRIBUTE_group'])
View(adonres)




