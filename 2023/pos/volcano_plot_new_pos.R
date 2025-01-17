#read in libraries
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

setwd('/home/anastasiia/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Univariate/')
getwd()


# Read the data into R
dset<-read.csv('Normalised_Quant_table.csv', header = TRUE)
colnames(dset) <- gsub("X", "", colnames(dset))

#def only
colnames(dset)[1] <- "feature-ID"
#dset_t <- dset %>% column_to_rownames ('filename') %>% t() %>% as.data.frame %>% rownames_to_column("feature-ID")
data<- as.matrix(dset[, c(2,4,5,6,7,12,17,18,19,20,21,22,23,24,26,27)])

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

wt= data[,c(1,6,14,15,16)]
#wt <- as.numeric(wt)
#view(wt)

dKO <- data[,c(11,12,13)]
#dKO<- as.numeric(dKO)
#view(dKO)

# Compute the means of the samples of each condition
wt.mean = apply(wt, 1, mean)
pbt.mean = apply(pbt, 1, mean)
nrp.mean = apply(nrp, 1, mean)
dKO.mean = apply(dKO, 1, mean)

# Just get the maximum of all the means for plotting
limit = max(wt.mean, pbt.mean, nrp.mean, dKO.mean)
means <- cbind.data.frame(wt.mean, pbt.mean, nrp.mean, dKO.mean)

###############################################################
######### WT VS Pbt mut ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = pbt.mean / wt.mean
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
                          "WT VS Pbt mut"){
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

ggsave(file="test.svg", plot=plot, width=10, height=8)

write.csv(DEres.df, "~/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Univariate/pbt_wt.csv", row.names=FALSE)



###############################################################
######### WT VS NRP mut ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = nrp.mean / wt.mean
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
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "WT VS Nrp mut"){
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
plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)


ggsave(file="test.svg", plot=plot, width=10, height=8)

write.csv(DEres.df, "~/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Univariate/nrp_wt.csv.csv", row.names=FALSE)


###############################################################
######### PBT VS NRP mut ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = nrp.mean / pbt.mean
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
  x = pbt[i,] # WT of gene number i
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
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "Pbt VS Nrp mut"){
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


ggsave(file="test.svg", plot=plot, width=10, height=8)

write.csv(DEres.df, "~/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Univariate/pbt_nrp.csv", row.names=FALSE)





###############################################################
######### PBT VS dKO mut ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = pbt.mean / dKO.mean
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
  x = pbt[i,] # WT of gene number i
  y = dKO[i,] # KO of gene number i
  
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
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "Pbt VS dKO"){
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


plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)

ggsave(file="test.svg", plot=plot, width=10, height=8)

write.csv(DEres.df, "~/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Univariate/pbt_dKO.csv", row.names=FALSE)




###############################################################
######### NRP VS dKO mut ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = nrp.mean / dKO.mean
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
  x = nrp[i,] # WT of gene number i
  y = dKO[i,] # KO of gene number i
  
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
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "Nrp VS dKO"){
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


plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)

ggsave(file="test.svg", plot=plot, width=10, height=8)

write.csv(DEres.df, "~/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Univariate/nrp_dKO.csv", row.names=FALSE)


###############################################################
######### WT VS dKO mut ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = wt.mean / dKO.mean
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
  y = dKO[i,] # KO of gene number i
  
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
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "wt VS dKO"){
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


plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)

ggsave(file="test.svg", plot=plot, width=10, height=8)

write.csv(DEres.df, "~/Desktop/proteus mirabilis/new_method_pos_2-99_ALL/Univariate/wt_dKO.csv", row.names=FALSE)







