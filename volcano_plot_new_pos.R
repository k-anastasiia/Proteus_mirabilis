#read in libraries
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

setwd('/home/anastasiia/Desktop/proteus mirabilis/new_method_pos/Univariate/')
getwd()


# Read the data into R
dset<-read.csv('Normalised_Quant_table.csv', header = TRUE)
colnames(dset) <- gsub("X", "", colnames(dset))

#def only
colnames(dset)[1] <- "feature-ID"
#dset_t <- dset %>% column_to_rownames ('filename') %>% t() %>% as.data.frame %>% rownames_to_column("feature-ID")
data<- as.matrix(dset[, c(2,4,5,6,7,12,17,18,19,20,21,22,24,25)])

featureIDs <- dset[,1]
view(featureIDs)
# Check the loaded dataset
# Dimension of the dataset
dim(data) 
head(data)


# Separate the two conditions into two smaller data frames
pbt = data[,2:5]
#pbt <- as.matrix(pbt)
#view(pbt)

nrp = data[,7:10]
#nrp<- as.numeric(nrp)
#view(nrp)

wt= data[,c(1,6,12,13,14)]
#wt <- as.numeric(wt)
#view(wt)

dKO <- data[,11]
#dKO<- as.numeric(dKO)
#view(dKO)

# Compute the means of the samples of each condition
wt.mean = apply(wt, 1, mean)
pbt.mean = apply(pbt, 1, mean)
nrp.mean = apply(nrp, 1, mean)
#dKO.mean = apply(dKO, 1, mean)

# Just get the maximum of all the means for plotting
limit = max(wt.mean, pbt.mean, nrp.mean)
means <- cbind.data.frame(wt.mean, pbt.mean, nrp.mean)

###############################################################
######### WT VS PBT ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = wt.mean / pbt.mean
log2fold <- log2(fold)
#add to results dataframe
DEres.df$logFC <- log2fold

# Generate histogram of the fold differences
foldhist <- ggplot(DEres.df, aes(x=logFC)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Fold Change Distribution") +
  labs(x= "Log fold change in abundance between groups", y = "Count") +
  theme_minimal() 
foldhist

#compute stat sig - t-test
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics
for(i in 1 : nrow(data)) { # For each gene : 
  x = wt[i,] # WT of gene number i
  y = pbt[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y, na.rm=TRUE, paired=FALSE, exact=FALSE)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}
# add results to results dataframe
DEres.df$pvalue <- pvalue

# Generate -log10(pvalues)
pvallogged <- -log10(pvalue)
# add results to results dataframe
DEres.df$logpval <- pvallogged
# Generate histogram of -log10(pvalues)
pval_hist <- ggplot(DEres.df, aes(x=logpval)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("-log10(pvalue) Distribution") +
  labs(x= "Significance levels occuring in D.E. results", y = "Count") +
  xlim(0,4) +
  theme_minimal() 
pval_hist

# Volcano Plots:
# put the biological significance (fold changes)
# and statistical significance (p-value) in one plot
# Generate the volcano plot
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "Volcano Plot of Features"){
  # Inputs: ##############################################
  
  # results.df: dataframe containing columns: 
  #             1) "featureIDs": feature IDs
  #             2) "logFC": foldchange values from DE analysis
  #             3) "logpval": log-transformed p-values from DE analysis
  # fold cutoff: significance threshold to render visually on plot;
  #               denotes fold difference between mutant and wildtype
  #               also referred to as "biologial signal"
  # p_value_cutoff: significance threshold to render visually on plot;
  #                 denotes statistical significance
  ########################################################
  # create factor denoting differential expression status
  stats.df <- results.df %>% mutate(diffexpressed = 
                                      ifelse(logFC < -(fold_cutoff) & 
                                               logpval >= -log10(pvalue_cutoff),
                                             "Significant, downregulated", 
                                             ifelse(logFC > fold_cutoff & 
                                                      logpval >= -log10(pvalue_cutoff), 
                                                    "Significant, upregulated", "Normally expressed")))
  
  stats.df$diffexpressed <- as.factor(stats.df$diffexpressed)
  # Base plot
  plot <- ggplot(stats.df, aes(x = logFC, y = logpval, col = diffexpressed, text = featureIDs)) +
    geom_point() +
    ggtitle(title) + 
    ylab("-log10(pvalue)") +
    geom_vline(xintercept=c(-fold_cutoff, fold_cutoff), col="red") +
    geom_hline(yintercept=-log10(pvalue_cutoff), col="darkturquoise") +
    scale_color_manual(values=c("black", "blue", "red"), 
                       name = "Differential expression\nstatus") +
    theme_minimal()
  return(plot)
  
}

library("plotly")
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)

install.packages('svglite')
ggsave(file="test.svg", plot=plot, width=10, height=8)

write.csv(DEres.df, "~/Desktop/test_R_statistics_for_volanco_wilcoxon_wtko_neg.csv", row.names=FALSE)





###############################################################
######### WT VS NRP ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = wt.mean / nrp.mean
log2fold <- log2(fold)
#add to results dataframe
DEres.df$logFC <- log2fold

# Generate histogram of the fold differences
foldhist <- ggplot(DEres.df, aes(x=logFC)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Fold Change Distribution") +
  labs(x= "Log fold change in abundance between groups", y = "Count") +
  theme_minimal() 
foldhist

#compute stat sig - t-test
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics
for(i in 1 : nrow(data)) { # For each gene : 
  x = wt[i,] # WT of gene number i
  y = nrp[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y, na.rm=TRUE, paired=FALSE, exact=FALSE)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}
# add results to results dataframe
DEres.df$pvalue <- pvalue

# Generate -log10(pvalues)
pvallogged <- -log10(pvalue)
# add results to results dataframe
DEres.df$logpval <- pvallogged
# Generate histogram of -log10(pvalues)
pval_hist <- ggplot(DEres.df, aes(x=logpval)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("-log10(pvalue) Distribution") +
  labs(x= "Significance levels occuring in D.E. results", y = "Count") +
  xlim(0,4) +
  theme_minimal() 
pval_hist

# Volcano Plots:
# put the biological significance (fold changes)
# and statistical significance (p-value) in one plot
# Generate the volcano plot
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "Volcano Plot of Features"){
  # Inputs: ##############################################
  
  # results.df: dataframe containing columns: 
  #             1) "featureIDs": feature IDs
  #             2) "logFC": foldchange values from DE analysis
  #             3) "logpval": log-transformed p-values from DE analysis
  # fold cutoff: significance threshold to render visually on plot;
  #               denotes fold difference between mutant and wildtype
  #               also referred to as "biologial signal"
  # p_value_cutoff: significance threshold to render visually on plot;
  #                 denotes statistical significance
  ########################################################
  # create factor denoting differential expression status
  stats.df <- results.df %>% mutate(diffexpressed = 
                                      ifelse(logFC < -(fold_cutoff) & 
                                               logpval >= -log10(pvalue_cutoff),
                                             "Significant, downregulated", 
                                             ifelse(logFC > fold_cutoff & 
                                                      logpval >= -log10(pvalue_cutoff), 
                                                    "Significant, upregulated", "Normally expressed")))
  
  stats.df$diffexpressed <- as.factor(stats.df$diffexpressed)
  # Base plot
  plot <- ggplot(stats.df, aes(x = logFC, y = logpval, col = diffexpressed, text = featureIDs)) +
    geom_point() +
    ggtitle(title) + 
    ylab("-log10(pvalue)") +
    geom_vline(xintercept=c(-fold_cutoff, fold_cutoff), col="red") +
    geom_hline(yintercept=-log10(pvalue_cutoff), col="darkturquoise") +
    scale_color_manual(values=c("black", "blue", "red"), 
                       name = "Differential expression\nstatus") +
    theme_minimal()
  return(plot)
  
}

library("plotly")
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)

install.packages('svglite')
ggsave(file="test.svg", plot=plot, width=10, height=8)

write.csv(DEres.df, "~/Desktop/test_R_statistics_for_volanco_wilcoxon_wtko_neg.csv", row.names=FALSE)




##############END




















################TEST FOR NORMALITY#############################
library("dplyr")
library("tidyverse")
#my_data3 <- read_csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797.csv") %>% column_to_rownames("group")
my_data2 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_2.csv")
my_data2 <- my_data2 %>% column_to_rownames(., var = 'group')
sapply(my_data2, class) 
my_data2$norm_peak_area <- as.numeric(as.character(my_data2$norm_peak_area))
my_matrix2 <- data.matrix(my_data2, rownames.force = NA)
hist(my_matrix2)
shapiro.test(my_matrix2)

#A7
my_data_A7 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_A7.csv")
my_data_A7 <- my_data_A7 %>% column_to_rownames(., var = 'group')
my_matrixA7 <- data.matrix(my_data_A7, rownames.force = NA)
hist(my_matrixA7)
shapiro.test(my_matrixA7)

#KO
my_data_A7 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_KO.csv")
my_data_A7 <- my_data_A7 %>% column_to_rownames(., var = 'group')
my_matrixA7 <- data.matrix(my_data_A7, rownames.force = NA)
hist(my_matrixA7)
shapiro.test(my_matrixA7)

#WT
my_data_A7 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_WT.csv")
my_data_A7 <- my_data_A7 %>% column_to_rownames(., var = 'group')
my_matrixA7 <- data.matrix(my_data_A7, rownames.force = NA)
hist(my_matrixA7)
shapiro.test(my_matrixA7)

#WTlim
my_data_A7 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797_WT2.csv")
my_data_A7 <- my_data_A7 %>% column_to_rownames(., var = 'group')
my_matrixA7 <- data.matrix(my_data_A7, rownames.force = NA)
hist(my_matrixA7)
shapiro.test(my_matrixA7)


############################################################################
#####box plot with wilcoxon and krustal-wallis stats for feature m/z 797####
library("ggpubr")
library("datatable")
my_data2 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797.csv")
p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10")
p

my_comparisons <- list( c("A7", "KO"), c("KO", "WT"), c("A7", "WT"), c("WT", "WT_limited") )
fig <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  stat_compare_means(label.y = .0012)
fig
ggsave(file="wilcox_boxplot.svg", plot=fig, width=10, height=8)

############################################################################
#####box plot with anova and stats for feature m/z 797####

library("ggpubr")
library("data.table")
my_data2 <- read.csv("~/Desktop/Processed_Norm_Quant_2032022_neg_mz_797.csv")
p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10")
p

my_comparisons <- list( c("A7", "KO"), c("KO", "WT"), c("A7", "WT"), c("WT", "WT_limited") )
fig <- p + stat_compare_means(comparisons = my_comparisons, method = "t.test") + stat_compare_means(method = "anova", label.y = .0012)
fig

ggsave(file="anova_boxplot.svg", plot=fig, width=10, height=8)
