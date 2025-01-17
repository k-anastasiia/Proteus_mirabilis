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


setwd('~/Desktop/Proteus_mirabilis_2024/Positive_mode/Multivariate')
getwd()

metadata<- read_csv("Pos_Proteus_mirabilis_2024_metadata_no_QC.csv")
metadata <- metadata[-1:-6,] #exclude blank
metadata <- metadata[-13,] #exclude
scaled_table <- read_csv("CLR_Scaled_pos.csv")
scaled_table <- scaled_table [-45,] #exclude

colnames(scaled_table)[1] <- "filename"
data_merge <- metadata %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
   column_to_rownames("filename")

res.pca <- data_merge %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
#fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Diet, addEllipses = TRUE,repel = TRUE, legend.title = "Diet",
             #palette = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"))
fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("green4","royalblue","red3", "grey", "orange"),
             geom.ind = "point", pointsize = 2)



pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
#view(pc_loadings)
pc_loadings[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(pc_loadings$'feature', '_', 3)
fviz_pca_var(res.pca, select.var = list(contrib = 50), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#dKO vs pbt vs nrp
metadata_d2 <- metadata[-37:-48,]
data_merge_d2 <- metadata_d2 %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_d2 %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_d2$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue","green4"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d2$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

#pc_loadings_d[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(pc_loadings_d$'feature', '_', 3)
#pc_loadings_d <- as.matrix(pc_loadings_d[, c(1,20,18,19)])

#write.csv(pc_loadings_d, "pc_loadings_d.csv",row.names = TRUE)
fviz_pca_var(res.pca_d2, select.var = list(contrib = 50), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


####################PERMANOVA############################
distm <- dist(data_merge_d2, method = "euclidean")# compute distance

adonres <- adonis2(distm ~ data_merge_d2[,colnames(data_merge_d2) == 'ATTRIBUTE_group'])
View(adonres)





#BioREPLICATE1
metadata_b1 <- metadata[(c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46)),]

data_merge_b1 <- metadata_b1 %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca <- data_merge_b1%>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
#fviz_pca_ind(res.pca, col.ind = data_merge$ATTRIBUTE_Diet, addEllipses = TRUE,repel = TRUE, legend.title = "Diet",
#palette = c("Normal"="green4", "Deficient"="red3","0.6% Fe"="royalblue1"))
fviz_pca_ind(res.pca, col.ind = data_merge_b1$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("green4","royalblue","red3", "grey"),
             geom.ind = "point", pointsize = 2)


pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
pc_loadings[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(pc_loadings$'feature', '_', 3)
fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#dKO vs pbt vs nrp
metadata_b1_d2 <- metadata_b1[1:12,]
data_merge_b1_d2 <- metadata_b1_d2 %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_b1_d2 %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_b1_d2$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue","green4"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d2$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

fviz_pca_var(res.pca_d2, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


# pbt vs nrp
metadata_b1_d <- metadata_b1[-1:-4,]
metadata_b1_d <- metadata_b1_d[-9:-12,]
data_merge_b1_d <- metadata_b1_d %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_b1_d %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_b1_d$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d2$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

fviz_pca_var(res.pca_d2, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

####################PERMANOVA############################
distm <- dist(data_merge__b1d2, method = "euclidean")# compute distance

adonres <- adonis2(distm ~ data_merge_b1_d2[,colnames(data_merge_b1_d2) == 'ATTRIBUTE_group'])
View(adonres)







#BioREPLICATE2
metadata_b2 <- metadata[(c(2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47)),]

data_merge_b2 <- metadata_b2 %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca <- data_merge_b2%>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
fviz_pca_ind(res.pca, col.ind = data_merge_b2$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("green4","royalblue","red3", "grey"),
             geom.ind = "point", pointsize = 2)


pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
pc_loadings[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(pc_loadings$'feature', '_', 3)
fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#dKO vs pbt vs nrp
metadata_b2_d2 <- metadata_b2[1:12,]
data_merge_b2_d2 <- metadata_b2_d2 %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_b2_d2 %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_b2_d2$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue","green4"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d2$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

fviz_pca_var(res.pca_d2, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


# pbt vs nrp
metadata_b2_d <- metadata_b2[-1:-4,]
metadata_b2_d <- metadata_b2_d[-9:-12,]
data_merge_b2_d <- metadata_b2_d %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_b2_d %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_b2_d$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d2$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

fviz_pca_var(res.pca_d2, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)



#BioREPLICATE3
metadata_b3 <- metadata[(c(3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48)),]

data_merge_b3 <- metadata_b3 %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca <- data_merge_b3%>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)
fviz_pca_ind(res.pca, col.ind = data_merge_b3$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("green4","royalblue","red3", "grey"),
             geom.ind = "point", pointsize = 2)


pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
pc_loadings[c('row_ID', 'm/z', 'RT')] <- str_split_fixed(pc_loadings$'feature', '_', 3)
fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#dKO vs pbt vs nrp
metadata_b3_d2 <- metadata_b2[1:12,]
data_merge_b3_d2 <- metadata_b3_d2 %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_b3_d2 %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_b3_d2$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue","green4"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d2$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

fviz_pca_var(res.pca_d2, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


# pbt vs nrp
metadata_b3_d <- metadata_b3[-1:-4,]
metadata_b3_d <- metadata_b3_d[-9:-12,]
data_merge_b3_d <- metadata_b3_d %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_b3_d %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_b3_d$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d2$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

fviz_pca_var(res.pca_d2, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# dko vs nrp
metadata_b3_dz <- metadata_b3[-13:-16,]
metadata_b3_dz <- metadata_b3_dz[-5:-8,]
data_merge_b3_dz <- metadata_b3_dz %>% dplyr::select("filename", "ATTRIBUTE_group") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

res.pca_d2 <- data_merge_b3_dz %>% dplyr::select(-ATTRIBUTE_group) %>% prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca_d2)
fviz_pca_ind(res.pca_d2, col.ind = data_merge_b3_dz$ATTRIBUTE_group, addEllipses = TRUE, repel = TRUE, 
             legend.title = "Group", palette = c("red3", "royalblue"),
             geom.ind = "point", pointsize = 2)

pc_loadings_d <- res.pca_d2$rotation
pc_loadings_d <- as_tibble(pc_loadings_d, rownames = "feature")

fviz_pca_var(res.pca_d2, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


