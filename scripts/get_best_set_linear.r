#!/usr/bin/env Rscript

#Load required packages

rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(DescTools)
library(pROC)
library(caret)
library(ggthemes)
library(ggcorrplot)

##############################
#######Read in the files

out_PRSet = args[1]
out_PRSice_on_Set = args[2]
out_lassosum_on_set = args[3]
out_PRScs_on_set = args[4]
out_LDpred2_on_set = args[5]
out_lasso2_on_set = args[6]
out_comparison = args[7]

test_prefix = args[8]
test_prefix_rsID = args[9]
training_prefix = args[10]
training_prefix_rsID = args[11]

cov_test = args[12]
cov_training = args[13]
pheno_file = args[14]

covariates = args[15]
cov_file = args[16]

cross_validation = args[17]

##############################
#######Define functions

#### PC correction

PC_correction <- function(data, score, score_res) {
  tip <- paste(score, "~", PCs, sep="")
  alpha <- lm(tip, data = data)
  beta <- residuals(alpha, "response")
  data[[score_res]] <- beta
  return(data)
}

#### Standardize

standardize <- function(data_training,data_test, score, score_sc) {
  mean_training <- mean(data_training[[score]])
  sd_training <- sd(data_training[[score]])
  data_test[[score_sc]] <- (data_test[[score]] - mean_training)/sd_training
  return(data_test)
}

#### Figures


##### Multiset_barplot
multiset_barplot <- function(data, Parameters, R2, logP, out){
    plot <- ggplot(data = data, aes_string(x = R2, y = paste0("reorder(", Parameters, ", -", R2, ")"))) +
        geom_bar(stat = "identity", aes(fill = logP), show.legend = TRUE) +  # Use logP for color fill
        scale_fill_gradient(low = "blue", high = "yellow", name = "-log10(P-value)") +  # Set color gradient
        labs(
        x = expression(R^2 ~ "(Variance Explained)"),
        y = "Gene Sets") +
        theme_few(base_size = 18) +
        theme(axis.text.y = element_text(size = 16),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
             )
    ggsave(out, plot, width=8, height=6)
}


##### Correlation plot
correlation_plot <- function(data, set_columns, set_labels, output_file) {
  # Compute correlation matrix for pathway scores
  cor_matrix <- cor(data[, set_columns], use = "pairwise.complete.obs", method = "pearson")
  
  # Compute p-values for significance filtering
  p_matrix <- matrix(NA, ncol = ncol(cor_matrix), nrow = nrow(cor_matrix))
  colnames(p_matrix) <- colnames(cor_matrix)
  rownames(p_matrix) <- rownames(cor_matrix)
  for (i in 1:ncol(cor_matrix)) {
    for (j in 1:nrow(cor_matrix)) {
      if (i != j) {  # Avoid self-correlation
        test <- cor.test(data[[set_columns[i]]], data[[set_columns[j]]], method = "pearson")
        p_matrix[i, j] <- test$p.value
      }
    }
  }

  # Bonferroni correction for multiple testing
  num_tests <- (ncol(cor_matrix) * (ncol(cor_matrix) - 1)) / 2  # Total number of pairwise tests
  bonferroni_threshold <- 0.05 / num_tests  # Adjusted significance threshold

  # Mask non-significant correlations (p > 0.05)
  cor_matrix[p_matrix > bonferroni_threshold] <- NA
  
  # Rename rows & columns in correlation matrix for better readability
  colnames(cor_matrix) <- set_labels
  rownames(cor_matrix) <- set_labels

  # Generate the correlation plot
  svg(output_file, width = 37, height = 35)
  print(ggcorrplot(cor_matrix, method = "square", type = "lower", # method = "circle" also works
             lab = TRUE, lab_size = 16, 
             colors = c("darkred", "grey89", "steelblue3"),
             #title = "Correlation Plot of PRS Pathways",
             show.diag = TRUE) +
    theme_few(base_size = 32) +
    labs(x = NULL, y = NULL) +  # Remove axis labels
    theme(axis.text.x = element_text(size = 62, angle = 45, hjust = 1),  # Adjust X-axis text size
          axis.text.y = element_text(size = 62), # Adjust Y-axis text size
          legend.text = element_text(size = 60),
          legend.title = element_text(size = 66))+
    guides(
      fill = guide_colorbar(
        barwidth = 4,     # thickness of color bar
        barheight = 40,   # length of color bar
        title.theme = element_text(size = 62),
        label.theme = element_text(size = 60))))
  dev.off()
}


##################
##General files###
##################

PCs_test <- read.table(cov_test, header=T)
PCs_training <- read.table(cov_training, header=T)

PC_only <- select(PCs_test, -FID, -IID)
PCs <- paste(colnames(PC_only), collapse = '+')

pheno <- read.table(pheno_file, header=T)
colnames(pheno) <- c("FID", "IID", "PHENO")
pheno <- pheno[!is.na(pheno$PHENO),]

if (cov_file != "" && file.exists(cov_file)) {
  cov_data <- read.table(cov_file, header = T)
  if(covariates != ""){
    cov_analysis <- "TRUE"
    covariate_list <- strsplit(covariates, ",")[[1]]
  } else {
    cat("No covariates given. \n")
    cov_analysis <- "FALSE"
  }
  
} else {
  cat("No covariate file provided.\n")
  cov_analysis <- "FALSE"
}

file_path_out_comparison <- file.path(out_comparison, "Regression_results_best_per_tool") # Read the comparison file
comparison <- read.table(file_path_out_comparison, header=TRUE, sep="\t", stringsAsFactors=FALSE)

compare_tools = TRUE

##################
#####PRSet########
##################

if (file.exists(paste0(out_PRSet, "/", test_prefix, ".all_score"))) {

  cat("Start get best set... PRSet...\n")

  #### Scores

  PRS_test_PRSet <- read.table(paste0(out_PRSet, "/", test_prefix, ".all_score"), header=T)
  PRS_training_PRSet <- read.table(paste0(out_PRSet, "/", training_prefix, ".all_score"), header=T)
  # Keep only the first occurrence of 'FID' and 'IID', and remove subsequent duplicates
  PRS_test_PRSet <- PRS_test_PRSet[, !grepl("^(FID|IID)\\.[0-9]+$", colnames(PRS_test_PRSet))]
  PRS_training_PRSet <- PRS_training_PRSet[, !grepl("^(FID|IID)\\.[0-9]+$", colnames(PRS_training_PRSet))]


  ##############################
  #######PC correction

  #### Get Pt

  #get_Pt <- select(PRS_test_PRSet, -FID, -IID)
  get_Pt <- select(PRS_test_PRSet, -starts_with("FID"), -starts_with("IID"))
  Pt <- colnames(get_Pt)
  Pt_res <- paste(Pt, "_res", sep="")

  #### Merge files with PCs

  PRS_test_PRSet_PC <- merge(PRS_test_PRSet, PCs_test, by=c("IID", "FID"))
  PRS_training_PRSet_PC <- merge(PRS_training_PRSet, PCs_training, by=c("IID", "FID"))

  #### Merge test and training for PC correction

  PRS_test_PRSet_PC$Group <- "test"
  PRS_training_PRSet_PC$Group <- "training"
  merged_PRSet <- rbind(PRS_test_PRSet_PC, PRS_training_PRSet_PC)
  ##### Account for NA columns
  na_col_PRSet <- colSums(is.na(merged_PRSet))
  all_na_cols_PRSet <- which(na_col_PRSet == nrow(merged_PRSet))
  if(length(all_na_cols_PRSet) > 0) {
    merged_PRSet <- merged_PRSet[, -all_na_cols_PRSet]
  } else {
    # If all_na_cols is empty, just use the original dataframe
    merged_PRSet <- merged_PRSet
  }
  ##### Account for rows that have all NA
  na_rows_PRSet <- apply(merged_PRSet, 1, function(row) any(is.na(row)))
  merged_PRSet <- merged_PRSet[!na_rows_PRSet,]

  #### PC correction

  for (i in 1:length(Pt)) {
    score <- Pt[i]
    score_res <- Pt_res[i]
    merged_PRSet <- PC_correction(data=merged_PRSet, score=score, score_res=score_res)
  }

  #### Split again to then do standardization
  PRS_test_PRSet_PC <- merged_PRSet[(merged_PRSet$Group=="test"),]
  PRS_test_PRSet_PC <- PRS_test_PRSet_PC[, !colnames(PRS_test_PRSet_PC) %in% "Group"]
  PRS_training_PRSet_PC <- merged_PRSet[(merged_PRSet$Group=="training"),]
  PRS_training_PRSet_PC <- PRS_training_PRSet_PC[, !colnames(PRS_training_PRSet_PC) %in% "Group"]


  ##############################
  #######Standardization

  #### Standardize

  Pt_sc <- c(Pt, Pt_res)
  Pt_sc_2 <- paste(Pt_sc, "_sc", sep="")

  for (i in 1:length(Pt_sc)) {
    score <- Pt_sc[i]
    score_sc <- Pt_sc_2[i]
    PRS_test_PRSet_PC <- standardize(data_training=PRS_training_PRSet_PC, data_test=PRS_test_PRSet_PC, score=score, score_sc=score_sc)
    PRS_training_PRSet_PC <- standardize(data_training=PRS_training_PRSet_PC, data_test=PRS_training_PRSet_PC, score=score, score_sc=score_sc)
  }

  #### Get only scores to write to a file

  PRS_test_PRSet_scores <- colnames(PRS_test_PRSet_PC)[!(colnames(PRS_test_PRSet_PC) %in% colnames(PC_only))]
  PRS_test_PRSet <- PRS_test_PRSet_PC[, PRS_test_PRSet_scores]
  PRS_training_PRSet_scores <- colnames(PRS_training_PRSet_PC)[!(colnames(PRS_training_PRSet_PC) %in% colnames(PC_only))]
  PRS_training_PRSet <- PRS_training_PRSet_PC[, PRS_training_PRSet_scores]

  write.table(PRS_test_PRSet, paste0(out_PRSet, "/", test_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")
  write.table(PRS_training_PRSet, paste0(out_PRSet, "/", training_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")


  ##############################
  #######Linear regression

  Pt_res_sc <- paste(Pt, "_res_sc", sep="")
  PRS_test_PRSet_pheno <- merge(PRS_test_PRSet, pheno, by=c("IID", "FID"))
  PRS_test_PRSet_pheno$PHENO <- as.numeric(PRS_test_PRSet_pheno$PHENO)

  #### Perform the regression

  if (cross_validation=="FALSE"){

    Regression_results_PRSet <- data.frame(matrix(ncol=6, nrow=0))

    for (i in 1:length(Pt_res_sc)) {
      a <- Pt_res_sc[i]
      thresh <- Pt[i]  
      tip <- paste("PHENO", "~", a,sep="")
      alpha <- lm(tip, data = PRS_test_PRSet_pheno)
      R2 <- summary(alpha)$adj.r.squared
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      p <- coef(summary(alpha))[2,4]
      Regression_results_PRSet <- rbind(Regression_results_PRSet, c("PRSet",thresh,beta,SE,p,R2))
      colnames(Regression_results_PRSet) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2")
    }

    Regression_results_PRSet$R2 <- as.numeric(Regression_results_PRSet$R2)
    Regression_results_PRSet$P_value <- as.numeric(Regression_results_PRSet$P_value)
    Regression_results_PRSet$P_value <- ifelse(Regression_results_PRSet$P_value==0, 3.4e-314, Regression_results_PRSet$P_value)

    write.table(Regression_results_PRSet, paste0(out_PRSet, "/Regression_results_PRSet"), row.names=F, quote=F, sep="\t")

  } else if (cross_validation=="TRUE"){

    Regression_results_PRSet <- data.frame(matrix(ncol=7, nrow=0))
    train_control <- trainControl(method = "cv", number = 10)

    for (i in 1:length(Pt_res_sc)) {
      a <- Pt_res_sc[i]
      thresh <- Pt[i]  
      tip <- paste("PHENO", "~", a,sep="")
      alpha <- lm(tip, data = PRS_test_PRSet_pheno)
      R2 <- summary(alpha)$adj.r.squared
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      p <- coef(summary(alpha))[2,4]
      model <- train(as.formula(tip), data=PRS_test_PRSet_pheno, method="lm", trControl=train_control, metric="Rsquared")
      R2_CV <- as.numeric(model$results$Rsquared[1])
      Regression_results_PRSet <- rbind(Regression_results_PRSet, c("PRSet",thresh,beta,SE,p,R2,R2_CV))
      colnames(Regression_results_PRSet) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2", "R2_CV")
    }

    Regression_results_PRSet$R2 <- as.numeric(Regression_results_PRSet$R2)
    Regression_results_PRSet$R2_CV <- as.numeric(Regression_results_PRSet$R2_CV)
    Regression_results_PRSet$P_value <- as.numeric(Regression_results_PRSet$P_value)
    Regression_results_PRSet$P_value <- ifelse(Regression_results_PRSet$P_value==0, 3.4e-314, Regression_results_PRSet$P_value)

    write.table(Regression_results_PRSet, paste0(out_PRSet, "/Regression_results_PRSet"), row.names=F, quote=F, sep="\t")


  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ##############################
  #######Plot

  #### Multiset bar plot
  ## Calculate -log10(P) and filter out unnecessary rows if needed
  if (cross_validation=="FALSE"){
    prs_top <- Regression_results_PRSet %>%
    mutate(Threshold = gsub(".*_([^_]+)$", "\\1", Parameters)) %>%
    mutate(Parameters = gsub("_([^_]+)$", "", Parameters)) %>%
    mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
    arrange(desc(R2)) %>%                   # Order by variance explained
    head(7)                                 # Top 7
  } else if (cross_validation=="TRUE"){    
    prs_top <- Regression_results_PRSet %>%
    mutate(Threshold = gsub(".*_([^_]+)$", "\\1", Parameters)) %>%
    mutate(Parameters = gsub("_([^_]+)$", "", Parameters)) %>%
    mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
    arrange(desc(R2_CV)) %>%                   # Order by variance explained
    head(7)                                 # Top 7
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ## Make multiset bar plot
  if (cross_validation=="FALSE"){  
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2", logP="log", out=paste0(out_PRSet, "/multiset_bar_plot_R2_PRSet.svg"))  
  } else if (cross_validation=="TRUE"){     
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2_CV", logP="log", out=paste0(out_PRSet, "/multiset_bar_plot_R2_PRSet.svg"))
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }
  #### Correlation plot
  # Extract the columns data
  pathway_columns <- grep("_res_sc$", colnames(PRS_test_PRSet_pheno), value = TRUE)
  # Rename pathways: Remove "_res_sc" and replace "_" with " "
  pathway_labels <- gsub("_[^_]+_res_sc$", "_res_sc", pathway_columns)
  pathway_labels <- gsub("_", " ", gsub("_res_sc", "", pathway_labels))
  # Make plot
  correlation_plot(PRS_test_PRSet_pheno, pathway_columns, pathway_labels, paste0(out_PRSet, "/correlation_PRSet.svg"))

  ##############################
  #######Covariates

  if (cov_analysis == TRUE) {
    Regression_covariates_PRSet <- data.frame(matrix(ncol=6, nrow=0))
    for (i in 1:length(Pt_res_sc)) {
      a <- Pt_res_sc[i]
      thresh <- Pt[i]
      for (covs in covariate_list){
        PRS_test_PRSet_cov <- merge(PRS_test_PRSet_pheno, cov_data, by=c("FID", "IID"))
        full_form <- paste("PHENO", "~", a, "+", covs, sep="")
        null_form <- paste("PHENO", "~", covs, sep="")
        full_lm <- lm(full_form, data = PRS_test_PRSet_cov)
        null_lm <- lm(null_form, data = PRS_test_PRSet_cov)
        R2_full <- summary(full_lm)$adj.r.squared
        R2_null <- summary(null_lm)$adj.r.squared
        R2_PRS <- R2_full - R2_null
        Regression_covariates_PRSet <- rbind(Regression_covariates_PRSet, c("PRSet", thresh, R2_full, R2_null, R2_PRS, covs))
        colnames(Regression_covariates_PRSet) <- c("Tool", "Parameters", "R2_PRS_and_cov" , "R2_cov_only", "R2_PRS_only", "included_covariates")
      }
    }
    write.table(Regression_covariates_PRSet, paste0(out_PRSet, "/Regression_results_PRSet_covariates"), row.names=F, quote=F, sep="\t")
  }

  cat("Done getting best set for PRSet.... PRSice....\n")

} else {
  compare_tools = FALSE
  message("Skipping get best set for PRSet.... PRSice....")
}


##################
#####PRSice#######
##################

if (file.exists(paste0(out_PRSice_on_Set, "/", test_prefix, ".all_score"))) {

  cat("Start get best set... PRSice on set...\n")

  #### Scores
  PRS_test_PRSice <- read.table(paste0(out_PRSice_on_Set, "/", test_prefix, ".all_score"), header=T)
  PRS_training_PRSice <- read.table(paste0(out_PRSice_on_Set, "/", training_prefix, ".all_score"), header=T)
  # Keep only the first occurrence of 'FID' and 'IID', and remove subsequent duplicates
  PRS_test_PRSice <- PRS_test_PRSice[, !grepl("^(FID|IID)\\.[0-9]+$", colnames(PRS_test_PRSice))]
  PRS_training_PRSice <- PRS_training_PRSice[, !grepl("^(FID|IID)\\.[0-9]+$", colnames(PRS_training_PRSice))]


  ##############################
  #######PC correction

  #### Get Pt

  get_Pt <- select(PRS_test_PRSice, -starts_with("FID"), -starts_with("IID"))
  Pt <- colnames(get_Pt)
  Pt_res <- paste(Pt, "_res", sep="")

  #### Merge files with PCs

  PRS_test_PRSice_PC <- merge(PRS_test_PRSice, PCs_test, by=c("IID", "FID"))
  PRS_training_PRSice_PC <- merge(PRS_training_PRSice, PCs_training, by=c("IID", "FID"))

  #### Merge test and training for PC correction

  PRS_test_PRSice_PC$Group <- "test"
  PRS_training_PRSice_PC$Group <- "training"
  merged_PRSice <- rbind(PRS_test_PRSice_PC, PRS_training_PRSice_PC)
  ##### Account for NA columns
  na_col_PRSice <- colSums(is.na(merged_PRSice))
  all_na_cols_PRSice <- which(na_col_PRSice == nrow(merged_PRSice))
  if(length(all_na_cols_PRSice) > 0) {
    merged_PRSice <- merged_PRSice[, -all_na_cols_PRSice]
  } else {
    # If all_na_cols is empty, just use the original dataframe
    merged_PRSice <- merged_PRSice
  }
  ##### Account for rows that have all NA
  na_rows_PRSice <- apply(merged_PRSice, 1, function(row) any(is.na(row)))
  merged_PRSice <- merged_PRSice[!na_rows_PRSice,]

  #### PC correction

  for (i in 1:length(Pt)) {
    score <- Pt[i]
    score_res <- Pt_res[i]
    merged_PRSice <- PC_correction(data=merged_PRSice, score=score, score_res=score_res)
  }

  #### Split again to then do standardization
  PRS_test_PRSice_PC <- merged_PRSice[(merged_PRSice$Group=="test"),]
  PRS_test_PRSice_PC <- PRS_test_PRSice_PC[, !colnames(PRS_test_PRSice_PC) %in% "Group"]
  PRS_training_PRSice_PC <- merged_PRSice[(merged_PRSice$Group=="training"),]
  PRS_training_PRSice_PC <- PRS_training_PRSice_PC[, !colnames(PRS_training_PRSice_PC) %in% "Group"]


  ##############################
  #######Standardization

  #### Standardize

  Pt_sc <- c(Pt, Pt_res)
  Pt_sc_2 <- paste(Pt_sc, "_sc", sep="")

  for (i in 1:length(Pt_sc)) {
    score <- Pt_sc[i]
    score_sc <- Pt_sc_2[i]
    PRS_test_PRSice_PC <- standardize(data_training=PRS_training_PRSice_PC, data_test=PRS_test_PRSice_PC, score=score, score_sc=score_sc)
    PRS_training_PRSice_PC <- standardize(data_training=PRS_training_PRSice_PC, data_test=PRS_training_PRSice_PC, score=score, score_sc=score_sc)
  }

  #### Get only scores to write to a file

  PRS_test_PRSice_scores <- colnames(PRS_test_PRSice_PC)[!(colnames(PRS_test_PRSice_PC) %in% colnames(PC_only))]
  PRS_test_PRSice <- PRS_test_PRSice_PC[, PRS_test_PRSice_scores]
  PRS_training_PRSice_scores <- colnames(PRS_training_PRSice_PC)[!(colnames(PRS_training_PRSice_PC) %in% colnames(PC_only))]
  PRS_training_PRSice <- PRS_training_PRSice_PC[, PRS_training_PRSice_scores]

  write.table(PRS_test_PRSice, paste0(out_PRSice_on_Set, "/", test_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")
  write.table(PRS_training_PRSice, paste0(out_PRSice_on_Set, "/", training_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")


  ##############################
  #######Linear regression

  Pt_res_sc <- paste(Pt, "_res_sc", sep="")
  PRS_test_PRSice_pheno <- merge(PRS_test_PRSice, pheno, by=c("IID", "FID"))
  PRS_test_PRSice_pheno$PHENO <- as.numeric(PRS_test_PRSice_pheno$PHENO)

  #### Perform the regression

  if (cross_validation=="FALSE"){

    Regression_results_PRSice <- data.frame(matrix(ncol=6, nrow=0))

    for (i in 1:length(Pt_res_sc)) {
      a <- Pt_res_sc[i]
      thresh <- Pt[i]  
      tip <- paste("PHENO", "~", a,sep="")
      alpha <- lm(tip, data = PRS_test_PRSice_pheno)
      R2 <- summary(alpha)$adj.r.squared
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      p <- coef(summary(alpha))[2,4]
      Regression_results_PRSice <- rbind(Regression_results_PRSice, c("PRSice",thresh,beta,SE,p,R2))
      colnames(Regression_results_PRSice) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2")
    }

    Regression_results_PRSice$R2 <- as.numeric(Regression_results_PRSice$R2)
    Regression_results_PRSice$P_value <- as.numeric(Regression_results_PRSice$P_value)
    Regression_results_PRSice$P_value <- ifelse(Regression_results_PRSice$P_value==0, 3.4e-314, Regression_results_PRSice$P_value)

    write.table(Regression_results_PRSice, paste0(out_PRSice_on_Set, "/Regression_results_PRSice"), row.names=F, quote=F, sep="\t")

  } else if (cross_validation=="TRUE"){

    Regression_results_PRSice <- data.frame(matrix(ncol=7, nrow=0))
    train_control <- trainControl(method = "cv", number = 10)

    for (i in 1:length(Pt_res_sc)) {
      a <- Pt_res_sc[i]
      thresh <- Pt[i]  
      tip <- paste("PHENO", "~", a,sep="")
      alpha <- lm(tip, data = PRS_test_PRSice_pheno)
      R2 <- summary(alpha)$adj.r.squared
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      p <- coef(summary(alpha))[2,4]
      model <- train(as.formula(tip), data=PRS_test_PRSice_pheno, method="lm", trControl=train_control, metric="Rsquared")
      R2_CV <- as.numeric(model$results$Rsquared[1])
      Regression_results_PRSice <- rbind(Regression_results_PRSice, c("PRSice",thresh,beta,SE,p,R2,R2_CV))
      colnames(Regression_results_PRSice) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2", "R2_CV")
    }

    Regression_results_PRSice$R2 <- as.numeric(Regression_results_PRSice$R2)
    Regression_results_PRSice$R2_CV <- as.numeric(Regression_results_PRSice$R2_CV)
    Regression_results_PRSice$P_value <- as.numeric(Regression_results_PRSice$P_value)
    Regression_results_PRSice$P_value <- ifelse(Regression_results_PRSice$P_value==0, 3.4e-314, Regression_results_PRSice$P_value)

    write.table(Regression_results_PRSice, paste0(out_PRSice_on_Set, "/Regression_results_PRSice"), row.names=F, quote=F, sep="\t")

  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ##############################
  #######Plot

  #### Multiset bar plot
  ## Calculate -log10(P) and filter out unnecessary rows if needed
  if (cross_validation=="FALSE"){
    prs_top <- Regression_results_PRSice %>%
      mutate(Threshold = gsub(".*_([^_]+)$", "\\1", Parameters)) %>%
      mutate(Parameters = gsub("_Pt.*$", "", Parameters)) %>%
      mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
      arrange(desc(R2)) %>%                   # Order by variance explained
      head(7)                                 # Top 7
  } else if (cross_validation=="TRUE"){    
      prs_top <- Regression_results_PRSice %>%
      mutate(Threshold = gsub(".*_([^_]+)$", "\\1", Parameters)) %>%
      mutate(Parameters = gsub("_Pt.*$", "", Parameters)) %>%
      mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
      arrange(desc(R2_CV)) %>%                # Order by variance explained
      head(7)
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ## Make multiset bar plot
  if (cross_validation=="FALSE"){
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2", logP="log", out=paste0(out_PRSice_on_Set, "/multiset_bar_plot_R2_PRSice.svg"))
  } else if (cross_validation=="TRUE"){  
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2_CV", logP="log", out=paste0(out_PRSice_on_Set, "/multiset_bar_plot_R2_PRSice.svg"))
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  #### Correlation plot
  # Extract the columns data
  pathway_columns <- grep("_res_sc$", colnames(PRS_test_PRSice_pheno), value = TRUE)
  # Rename pathways: Remove "Pt", "_res_sc" and replace "_" with " "
  pathway_labels <- gsub("_[^_]+_res_sc$", "", pathway_columns)  # Remove the part before _res_sc
  pathway_labels <- gsub("_Pt|_res_sc", "", pathway_labels)      # Remove _Pt and _res_sc
  pathway_labels <- gsub("_", " ", pathway_labels)               # Replace remaining underscores with spaces
  # Make plot
  correlation_plot(PRS_test_PRSice_pheno, pathway_columns, pathway_labels, paste0(out_PRSice_on_Set, "/correlation_PRSice.svg"))

  ##############################
  #######Covariates

  if (cov_analysis == TRUE) {
    Regression_covariates_PRSice <- data.frame(matrix(ncol=6, nrow=0))
    for (i in 1:length(Pt_res_sc)) {
      a <- Pt_res_sc[i]
      thresh <- Pt[i]
      for (covs in covariate_list){
        PRS_test_PRSice_cov <- merge(PRS_test_PRSice_pheno, cov_data, by=c("FID", "IID"))
        full_form <- paste("PHENO", "~", a, "+", covs, sep="")
        null_form <- paste("PHENO", "~", covs, sep="")
        full_lm <- lm(full_form, data = PRS_test_PRSice_cov)
        null_lm <- lm(null_form, data = PRS_test_PRSice_cov)
        R2_full <- summary(full_lm)$adj.r.squared
        R2_null <- summary(null_lm)$adj.r.squared
        R2_PRS <- R2_full - R2_null
        Regression_covariates_PRSice <- rbind(Regression_covariates_PRSice, c("PRSice", thresh, R2_full, R2_null, R2_PRS, covs))
        colnames(Regression_covariates_PRSice) <- c("Tool", "Parameters", "R2_PRS_and_cov" , "R2_cov_only", "R2_PRS_only", "included_covariates")
      }
    }
    write.table(Regression_covariates_PRSice, paste0(out_PRSice_on_Set, "/Regression_results_PRSice_covariates"), row.names=F, quote=F, sep="\t")
  }


  cat("Done getting best set for PRSice.... lassosum....\n")

} else {
  compare_tools = FALSE
  message("Skipping get best set for PRSice.... lassosum....")
}

##################
#####Lassosum#####
##################

if (file.exists(file.path(out_lassosum_on_set, "gene_sets.txt"))) {

  cat("Start get best set... lassosum on set...\n")

  ##############################
  #######Read in the files

  # Read in the gene sets
  gene_set_file <- file.path(out_lassosum_on_set, "gene_sets.txt") # Define the file path
  gene_sets <- read.table(gene_set_file, header = FALSE, stringsAsFactors = FALSE) # Read the file
  gene_sets <- gene_sets$V1 # Convert to a character vector

  # Get best s and L
  lassosum_row <- comparison[grep("^lassosum$", comparison$Tool), ]  # Filter for "lassosum" only
  best_s_lassosum <- strsplit(lassosum_row$Parameters, "_")[[1]][2] # Extract the second part after splitting by "_"
  best_L_lassosum <- strsplit(lassosum_row$Parameters, "_")[[1]][4]

  data_frames_training <- list()
  data_frames_test <- list()
  for (gene_set in gene_sets) {
    file_name <- paste(out_lassosum_on_set, "/", test_prefix, "_", gene_set, "_scores_", "s_", best_s_lassosum, ".txt", sep="")
    data <- read.table(file_name, header=T)
    data <- data[, colSums(data==0) !=(nrow(data))] #Delete columns where all values are 0
    data_frames_test[[gene_set]] <- data
    
    file_name <- paste(out_lassosum_on_set, "/", training_prefix, "_", gene_set, "_scores_", "s_", best_s_lassosum, ".txt", sep="")
    data <- read.table(file_name, header=T)
    data <- data[, colSums(data==0) !=(nrow(data))]
    data_frames_training[[gene_set]] <- data
  }


  ##############################
  #######PC correction

  data_frames_training_PC <- list()
  data_frames_test_PC <- list()
  PC_names <- colnames(PC_only)

  for (j in 1:length(data_frames_training)){
    training <- data_frames_training[[j]]
    training_PC <- merge(training, PCs_training, by=c("FID","IID"))
    training_PC$Group <- "training"
    test <- data_frames_test[[j]]
    test_PC <- merge(test, PCs_test, by=c("FID","IID"))
    test_PC$Group <- "test"
    common_cols <- intersect(names(training_PC), names(test_PC))
    training_PC <- training_PC[, common_cols]
    test_PC <- test_PC[, common_cols]
    merge_training_test <- rbind(training_PC, test_PC)
    na_col_lasso <- colSums(is.na(merge_training_test))
    all_na_cols_lasso <- which(na_col_lasso == nrow(merge_training_test))
    if(length(all_na_cols_lasso) > 0) {
      merge_training_test <- merge_training_test[, -all_na_cols_lasso]
    } else {
      # If all_na_cols is empty, just use the original dataframe
      merge_training_test <- merge_training_test
    }
    na_rows_lasso <- apply(merge_training_test, 1, function(row) any(is.na(row)))
    merge_training_test <- merge_training_test[!na_rows_lasso, ]
    scores <- select(merge_training_test, -IID, -FID, -Group, -all_of(PC_names))
    lambda <- colnames(scores)
    lambda_res <- paste(lambda, "_res", sep="")
  for (i in 1:length(lambda)) {
    score <- lambda[i]
    score_res <- lambda_res[i]
    merge_training_test <- PC_correction(data=merge_training_test, score=score, score_res=score_res)
  }
    training_residuals <- merge_training_test[(merge_training_test$Group=="training"),]
    training_residuals <- training_residuals[, !colnames(training_residuals) %in% "Group"]
    test_residuals <- merge_training_test[(merge_training_test$Group=="test"),]
    test_residuals <- test_residuals[, !colnames(test_residuals) %in% "Group"]

    data_frames_training_PC[[j]] <- training_residuals
    data_frames_test_PC[[j]] <- test_residuals
  }

  ##############################
  #######Standardization

  for (i in 1:length(data_frames_training_PC)) {
    training_PC <- data_frames_training_PC[[i]]
    test_PC <- data_frames_test_PC[[i]]
    training <- data_frames_training[[i]]
    scores <- select(training_PC, -IID, -FID, -all_of(PC_names))
    lambda_sc <- colnames(scores)
    lambda_sc_2 <- paste(lambda_sc, "_sc", sep="")
    
    for (j in 1:length(lambda_sc)){
      score <- lambda_sc[j]
      score_sc <- lambda_sc_2[j]
      training_PC <- standardize(data_training=training_PC, data_test=training_PC, score=score, score_sc=score_sc)
      test_PC <- standardize(data_training=training_PC, data_test=test_PC, score=score, score_sc=score_sc)
    }
    
    training_col <- colnames(training_PC)[!(colnames(training_PC) %in% colnames(PC_only))]
    training_PC <- training_PC[, training_col]
    test_col <- colnames(test_PC)[!(colnames(test_PC) %in% colnames(PC_only))]
    test_PC <- test_PC[, test_col]
    data_frames_training[[i]] <- training_PC
    data_frames_test[[i]] <- test_PC
    
    write.table(training_PC, paste0(out_lassosum_on_set, "/", training_prefix, "_", gene_sets[i], "_scaled_scores"), row.names=F, quote=F, sep="\t")
    write.table(test_PC, paste0(out_lassosum_on_set, "/", test_prefix, "_", gene_sets[i], "_scaled_scores"), row.names=F, quote=F, sep="\t")
  }


  ##############################
  #######Linear regression

  #### Perform the regression

  if (cross_validation=="FALSE"){

    Regression_results_lasso <- data.frame(matrix(ncol=11, nrow=0))

    for (j in 1:length(data_frames_test)){
      test <- data_frames_test[[j]]
      scores <- colnames(test)[grep("_res_sc$", colnames(test))]
      columns <- c("IID", "FID", scores)
      test_pheno <- test[, columns]
      test_pheno <- merge(test_pheno, pheno, by=c("IID", "FID"))
      colnames(test_pheno) <- c("IID", "FID",scores, "PHENO")
      test_pheno$PHENO <- as.numeric(test_pheno$PHENO)
      
      for (i in 1:length(scores)) {
        b <- scores[i]
        thresh <- paste(b)  
        tip <- paste("PHENO", "~", b,sep="")
        alpha <- lm(tip, data = test_pheno)
        R2 <- summary(alpha)$adj.r.squared
        beta <- coef(summary(alpha))[2,1]
        SE <- coef(summary(alpha))[2,2]
        p <- coef(summary(alpha))[2,4]
        Regression_results_lasso <- rbind(Regression_results_lasso, c("lassosum", paste0(gene_sets[j], "_", thresh),gene_sets[j], thresh,beta,SE,p,R2))
        colnames(Regression_results_lasso) <- c("Tool", "Parameters", "Set", "Lambda", "Beta", "SE", "P_value", "R2")
      }
    }

    Regression_results_lasso$R2 <- as.numeric(Regression_results_lasso$R2)
    Regression_results_lasso$P_value <- as.numeric(Regression_results_lasso$P_value)
    Regression_results_lasso$P_value <- ifelse(Regression_results_lasso$P_value==0, 3.4e-314, Regression_results_lasso$P_value)

    write.table(Regression_results_lasso, paste0(out_lassosum_on_set, "/Regression_results_lasso"), row.names=F, quote=F, sep="\t")

  } else if (cross_validation=="TRUE"){

    Regression_results_lasso <- data.frame(matrix(ncol=9, nrow=0))
    train_control <- trainControl(method = "cv", number = 10)

    for (j in 1:length(data_frames_test)){
      test <- data_frames_test[[j]]
      scores <- colnames(test)[grep("_res_sc$", colnames(test))]
      columns <- c("IID", "FID", scores)
      test_pheno <- test[, columns]
      test_pheno <- merge(test_pheno, pheno, by=c("IID", "FID"))
      colnames(test_pheno) <- c("IID", "FID",scores, "PHENO")
      test_pheno$PHENO <- as.numeric(test_pheno$PHENO)
      
      for (i in 1:length(scores)) {
        b <- scores[i]
        thresh <- paste(b)  
        tip <- paste("PHENO", "~", b,sep="")
        alpha <- lm(tip, data = test_pheno)
        R2 <- summary(alpha)$adj.r.squared
        beta <- coef(summary(alpha))[2,1]
        SE <- coef(summary(alpha))[2,2]
        p <- coef(summary(alpha))[2,4]
        model <- train(as.formula(tip), data=test_pheno, method="lm", trControl=train_control, metric="Rsquared")
        R2_CV <- as.numeric(model$results$Rsquared[1])
        Regression_results_lasso <- rbind(Regression_results_lasso, c("lassosum", paste0(gene_sets[j], "_", thresh),gene_sets[j], thresh,beta,SE,p,R2,R2_CV))
        colnames(Regression_results_lasso) <- c("Tool", "Parameters", "Set", "Lambda", "Beta", "SE", "P_value", "R2", "R2_CV")
      }
    }

    Regression_results_lasso$R2 <- as.numeric(Regression_results_lasso$R2)
    Regression_results_lasso$R2_CV <- as.numeric(Regression_results_lasso$R2_CV)
    Regression_results_lasso$P_value <- as.numeric(Regression_results_lasso$P_value)
    Regression_results_lasso$P_value <- ifelse(Regression_results_lasso$P_value==0, 3.4e-314, Regression_results_lasso$P_value)

    write.table(Regression_results_lasso, paste0(out_lassosum_on_set, "/Regression_results_lasso"), row.names=F, quote=F, sep="\t")

  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ##############################
  #######Plot

  #### Multiset bar plot
  ## Calculate -log10(P) and filter out unnecessary rows if needed
  if (cross_validation=="FALSE"){
    prs_top <- Regression_results_lasso %>%
      mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
      arrange(desc(R2)) %>%                   # Order by variance explained
      head(7)                                 # Top 7
  } else if (cross_validation=="TRUE"){    
    prs_top <- Regression_results_lasso %>%
      mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
      arrange(desc(R2_CV)) %>%                # Order by variance explained
      head(7)                                 # Top 7
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }
  
  ## Make multiset bar plot
  if (cross_validation=="FALSE"){    
    multiset_barplot(data=prs_top, Parameters="Set", R2="R2", logP="log", out=paste0(out_lassosum_on_set, "/multiset_bar_plot_R2_lasso.svg"))
  } else if (cross_validation=="TRUE"){    
    multiset_barplot(data=prs_top, Parameters="Set", R2="R2_CV", logP="log", out=paste0(out_lassosum_on_set, "/multiset_bar_plot_R2_lasso.svg"))
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  #### Correlation plot
  ## Create file with all scaled scores and pheno data using data_frames_test
  # Initialize an empty list to store transformed data frames
  transformed_list <- list()
  # Loop through each pathway in the list of data frames
  for (pathway in names(data_frames_test)) {
    df <- data_frames_test[[pathway]]
    # Select relevant columns
    df <- df[, c("FID", "IID", tail(names(df), 1))]
    # Rename the PRS column to the pathway name
    colnames(df)[3] <- pathway
    # Store the transformed dataframe
    transformed_list[[pathway]] <- df
  }
  # Merge all data frames by "FID" and "IID"
  final_data <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID"), all = TRUE), transformed_list)
  # Merge with pheno by "FID" and "IID"
  PRS_test_lasso_pheno <- merge(final_data, pheno, by = c("FID", "IID"), all.x = TRUE)


  # Extract the columns data
  pathway_columns <- setdiff(colnames(PRS_test_lasso_pheno), c("FID", "IID", "PHENO"))
  # Rename pathways: replace "_" with " "
  pathway_labels <- gsub("_", " ", pathway_columns)
  # Make plot
  correlation_plot(PRS_test_lasso_pheno, pathway_columns, pathway_labels, paste0(out_lassosum_on_set, "/correlation_lasso.svg"))

  ##############################
  #######Covariates

  if (cov_analysis == TRUE) {
    
    Regression_covariates_lasso <- data.frame(matrix(ncol=7, nrow=0))

    for (j in 1:length(data_frames_test)){
      test <- data_frames_test[[j]]
      scores <- colnames(test)[grep("_res_sc$", colnames(test))]
      columns <- c("IID", "FID", scores)
      test_pheno <- test[, columns]
      test_pheno <- merge(test_pheno, pheno, by=c("IID", "FID"))
      colnames(test_pheno) <- c("IID", "FID",scores, "PHENO")
      test_pheno$PHENO <- as.numeric(test_pheno$PHENO)

      for (i in 1:length(scores)) {
        b <- scores[i]
        thresh <- paste(b)
        for (covs in covariate_list){
          PRS_test_lasso_cov <- merge(test_pheno, cov_data, by=c("FID", "IID"))
          full_form <- paste("PHENO", "~", b, "+", covs, sep="")
          null_form <- paste("PHENO", "~", covs, sep="")
          full_lm <- lm(full_form, data = PRS_test_lasso_cov)
          null_lm <- lm(null_form, data = PRS_test_lasso_cov)
          R2_full <- summary(full_lm)$adj.r.squared
          R2_null <- summary(null_lm)$adj.r.squared
          R2_PRS <- R2_full - R2_null
          Regression_covariates_lasso <- rbind(Regression_covariates_lasso, c("lassosum",gene_sets[j],thresh, R2_full, R2_null, R2_PRS, covs))
          colnames(Regression_covariates_lasso) <- c("Tool", "Set", "Lambda", "R2_PRS_and_cov" , "R2_cov_only", "R2_PRS_only", "included_covariates")
        }
      }
    }
    write.table(Regression_covariates_lasso, paste0(out_lassosum_on_set, "/Regression_results_lasso_covariates"), row.names=F, quote=F, sep="\t")
  }


  cat("Done getting best set for lassosum.... PRS-CS....\n")

} else {
  compare_tools = FALSE
  message("Skipping get best set for lassosum.... PRS-CS....")
}

##################
#####PRS-CS#######
##################

if (length(list.files(out_PRScs_on_set)) > 0) {

  ##############################
  #######Read in the files

  # Get best phi
  PRSCS_row <- comparison[grep("^PRS-CS$", comparison$Tool), ]  # Filter for "PRS-CS"
  best_phi <- strsplit(PRSCS_row$Parameters, "_")[[1]][2] # Extract the second part after splitting by "_"

  PRS_test_PRS_CS <- PCs_test[,c("IID", "FID")]
  PRS_training_PRS_CS <- PCs_training[,c("IID","FID")]

  for (gene_set in gene_sets) {
    a <- gsub("\\.", "-", best_phi)
    b <- best_phi
    
    #test
    file_name <- paste(out_PRScs_on_set, "/", test_prefix, "_", gene_set, "_phi_", a, ".profile", sep="")
    data <- read.table(file_name, header=T)
    data <- data[,c(1,2,6)]
    colnames(data) <- c("FID", "IID", gene_set)
    PRS_test_PRS_CS <- merge(PRS_test_PRS_CS, data, by=c("IID", "FID"))

    #training
    file_name <- paste(out_PRScs_on_set, "/", training_prefix, "_", gene_set, "_phi_", a, ".profile", sep="")
    data <- read.table(file_name, header=T)
    data <- data[,c(1,2,6)]
    colnames(data) <- c("FID", "IID", gene_set)
    PRS_training_PRS_CS <- merge(PRS_training_PRS_CS, data, by=c("IID", "FID")) 
  }
  #print(head(PRS_test_PRS_CS))

  #### Merge files with PCs

  PRS_test_PRS_CS_PC <- merge(PRS_test_PRS_CS, PCs_test, by=c("IID", "FID"))
  PRS_training_PRS_CS_PC <- merge(PRS_training_PRS_CS, PCs_training, by=c("IID", "FID"))

  #### Merge them for PC correction
  PRS_test_PRS_CS_PC$Group <- "test"
  PRS_training_PRS_CS_PC$Group <- "training"
  merged_PRS_CS <- rbind(PRS_test_PRS_CS_PC, PRS_training_PRS_CS_PC)

  na_col_PRS_CS <- colSums(is.na(merged_PRS_CS))
  all_na_cols_PRS_CS <- which(na_col_PRS_CS == nrow(merged_PRS_CS))
  if(length(all_na_cols_PRS_CS) > 0) {
    merged_PRS_CS <- merged_PRS_CS[, -all_na_cols_PRS_CS]
  } else {
    # If all_na_cols is empty, just use the original dataframe
    merged_PRS_CS <- merged_PRS_CS
  }
  na_rows_PRS_CS <- apply(merged_PRS_CS, 1, function(row) any(is.na(row)))
  merged_PRS_CS <- merged_PRS_CS[!na_rows_PRS_CS, ]
  #print(head(merged_PRS_CS))

  ##############################
  #######PC correction

  for (gene_set in gene_sets) {
    score <- gene_set
    score_res <- paste0(gene_set, "_res")
    merged_PRS_CS <- PC_correction(data=merged_PRS_CS, score=score, score_res=score_res)
  }

  #### Split again to then do standardization
  PRS_test_PRS_CS_PC <- merged_PRS_CS[(merged_PRS_CS$Group=="test"),]
  PRS_test_PRS_CS_PC <- PRS_test_PRS_CS_PC[, !colnames(PRS_test_PRS_CS_PC) %in% "Group"]
  PRS_training_PRS_CS_PC <- merged_PRS_CS[(merged_PRS_CS$Group=="training"),]
  PRS_training_PRS_CS_PC <- PRS_training_PRS_CS_PC[, !colnames(PRS_training_PRS_CS_PC) %in% "Group"]

  ##############################
  #######Standardization

  for (gene_set in gene_sets) {
    score <- paste0(gene_set, "_res")
    score_sc <- paste(paste0(gene_set, "_res"), "_sc", sep="")
    PRS_test_PRS_CS_PC <- standardize(data_training=PRS_training_PRS_CS_PC, data_test=PRS_test_PRS_CS_PC, score=score, score_sc=score_sc)
    PRS_training_PRS_CS_PC <- standardize(data_training=PRS_training_PRS_CS_PC, data_test=PRS_training_PRS_CS_PC, score=score, score_sc=score_sc)
  }

  ### Get only scores to write to a file

  PRS_test_PRS_CS_scores <- colnames(PRS_test_PRS_CS_PC)[!(colnames(PRS_test_PRS_CS_PC) %in% colnames(PC_only))]
  PRS_test_PRS_CS <- PRS_test_PRS_CS_PC[, PRS_test_PRS_CS_scores]
  PRS_training_PRS_CS_scores <- colnames(PRS_training_PRS_CS_PC)[!(colnames(PRS_training_PRS_CS_PC) %in% colnames(PC_only))]
  PRS_training_PRS_CS <- PRS_training_PRS_CS_PC[, PRS_training_PRS_CS_scores]

  write.table(PRS_test_PRS_CS, paste0(out_PRScs_on_set, "/", test_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")
  write.table(PRS_training_PRS_CS, paste0(out_PRScs_on_set, "/", training_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")


  ##############################
  #######Linear regression

  PRS_test_PRS_CS_pheno <- merge(PRS_test_PRS_CS, pheno, by=c("IID", "FID"))
  PRS_test_PRS_CS_pheno$PHENO <- as.numeric(PRS_test_PRS_CS_pheno$PHENO)

  #### Perform the regression

  if (cross_validation=="FALSE"){

    Regression_results_PRS_CS <- data.frame(matrix(ncol=6, nrow=0))

    for (gene_set in gene_sets) {
      a <- paste(gene_set, "_res_sc", sep="")
      thresh <- gene_set
      tip <- paste("PHENO", "~", a, sep="")
      alpha <- lm(tip, data = PRS_test_PRS_CS_pheno)
      R2_full <- summary(alpha)$adj.r.squared
      p <- coef(summary(alpha))[2,4]
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      Regression_results_PRS_CS <- rbind(Regression_results_PRS_CS, c("PRS-CS",thresh,beta,SE,p,R2_full))
      colnames(Regression_results_PRS_CS) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2")
    }

    Regression_results_PRS_CS$R2 <- as.numeric(Regression_results_PRS_CS$R2)
    Regression_results_PRS_CS$P_value <- as.numeric(Regression_results_PRS_CS$P_value)
    Regression_results_PRS_CS$P_value <- ifelse(Regression_results_PRS_CS$P_value==0, 3.4e-314, Regression_results_PRS_CS$P_value)

    write.table(Regression_results_PRS_CS, paste0(out_PRScs_on_set, "/Regression_results_PRScs"), row.names=F, quote=F, sep="\t")

  } else if (cross_validation=="TRUE"){

    Regression_results_PRS_CS <- data.frame(matrix(ncol=7, nrow=0))
    train_control <- trainControl(method = "cv", number = 10)

      for (gene_set in gene_sets) {
        a <- paste(gene_set, "_res_sc", sep="")
        thresh <- gene_set
        tip <- paste("PHENO", "~", a, sep="")
        alpha <- lm(tip, data = PRS_test_PRS_CS_pheno)
        R2 <- summary(alpha)$adj.r.squared
        p <- coef(summary(alpha))[2,4]
        beta <- coef(summary(alpha))[2,1]
        SE <- coef(summary(alpha))[2,2]
        model <- train(as.formula(tip), data=PRS_test_PRS_CS_pheno, method="lm", trControl=train_control, metric="Rsquared")
        R2_CV <- as.numeric(model$results$Rsquared[1])
        Regression_results_PRS_CS <- rbind(Regression_results_PRS_CS, c("PRS-CS",thresh,beta,SE,p,R2,R2_CV))
        colnames(Regression_results_PRS_CS) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2", "R2_CV")
      }

      Regression_results_PRS_CS$R2 <- as.numeric(Regression_results_PRS_CS$R2)
      Regression_results_PRS_CS$R2_CV <- as.numeric(Regression_results_PRS_CS$R2_CV)
      Regression_results_PRS_CS$P_value <- as.numeric(Regression_results_PRS_CS$P_value)
      Regression_results_PRS_CS$P_value <- ifelse(Regression_results_PRS_CS$P_value==0, 3.4e-314, Regression_results_PRS_CS$P_value)

      write.table(Regression_results_PRS_CS, paste0(out_PRScs_on_set, "/Regression_results_PRScs"), row.names=F, quote=F, sep="\t")

  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ##############################
  #######Plot

  #### Multiset bar plot
  ## Calculate -log10(P) and filter out unnecessary rows if needed
  if (cross_validation=="FALSE"){
    prs_top <- Regression_results_PRS_CS %>%
      mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
      arrange(desc(R2)) %>%                   # Order by variance explained
      head(7)                                 # Top 7
  } else if (cross_validation=="TRUE"){  
    prs_top <- Regression_results_PRS_CS %>%
      mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
      arrange(desc(R2_CV)) %>%                   # Order by variance explained
      head(7)                                 # Top 7
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ## Make multiset bar plot
  if (cross_validation=="FALSE"){    
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2", logP="log", out=paste0(out_PRScs_on_set, "/multiset_bar_plot_R2_PRS-CS.svg"))
  } else if (cross_validation=="TRUE"){   
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2_CV", logP="log", out=paste0(out_PRScs_on_set, "/multiset_bar_plot_R2_PRS-CS.svg"))
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  #### Correlation plot
  # Extract the columns data
  pathway_columns <- grep("_res_sc$", colnames(PRS_test_PRS_CS_pheno), value = TRUE)
  # Rename pathways: Remove "_res_sc" and replace "_" with " "
  pathway_labels <- gsub("_", " ", gsub("_res_sc", "", pathway_columns))
  # Make plot
  correlation_plot(PRS_test_PRS_CS_pheno, pathway_columns, pathway_labels, paste0(out_PRScs_on_set, "/correlation_PRScs.svg"))

  ##############################
  #######Covariates

  if (cov_analysis == TRUE) {
    Regression_covariates_PRS_CS <- data.frame(matrix(ncol=6, nrow=0))
    for (gene_set in gene_sets) {
      a <- paste(gene_set, "_res_sc", sep="")
      thresh <- gene_set
      for (covs in covariate_list){
        PRS_test_PRS_CS_cov <- merge(PRS_test_PRS_CS_pheno, cov_data, by=c("FID", "IID"))
        full_form <- paste("PHENO", "~", a, "+", covs, sep="")
        null_form <- paste("PHENO", "~", covs, sep="")
        full_lm <- lm(full_form, data = PRS_test_PRS_CS_cov)
        null_lm <- lm(null_form, data = PRS_test_PRS_CS_cov)
        R2_full <- summary(full_lm)$adj.r.squared
        R2_null <- summary(null_lm)$adj.r.squared
        R2_PRS <- R2_full - R2_null
        Regression_covariates_PRS_CS <- rbind(Regression_covariates_PRS_CS, c("PRS-CS", thresh, R2_full, R2_null, R2_PRS, covs))
        colnames(Regression_covariates_PRS_CS) <- c("Tool", "Parameters", "R2_PRS_and_cov" , "R2_cov_only", "R2_PRS_only", "included_covariates")
      }
    }
    write.table(Regression_covariates_PRS_CS, paste0(out_PRScs_on_set, "/Regression_results_PRScs_covariates"), row.names=F, quote=F, sep="\t")
  }


  cat("Done getting best set for PRS-CS.... LDpred2....\n")

} else {
  compare_tools = FALSE
  message("Skipping get best set for PRS-CS.... LDpred2....")
}

##################
#####LDpred2######
##################

if (length(list.files(out_LDpred2_on_set)) > 0) {

  ##############################
  #######Select best LDpred2 model

  LDpred2_row <- comparison[grep("^LDpred2$", comparison$Tool), ]
  model_type <- ifelse(startsWith(LDpred2_row$Parameters, "inf"), "INF model",
                ifelse(startsWith(LDpred2_row$Parameters, "grid"), "GRID model",
                ifelse(startsWith(LDpred2_row$Parameters, "auto_grid"), "AUTO GRID model",
                ifelse(startsWith(LDpred2_row$Parameters, "auto"), "AUTO model", "Unknown"))))
  print("Best LDpred2 model:")
  print(model_type)


  ##############################
  #######Load scores

  ##### Inf
  if (model_type == "INF model") {
    PRS_inf_test <- PCs_test[,c("FID", "IID")]
    PRS_inf_training <- PCs_training[,c("FID","IID")]

    for (gene_set in gene_sets) {
      test <- read.table(paste0(out_LDpred2_on_set, "/", "pred_inf_", test_prefix_rsID, "_", gene_set), header=T)
      colnames(test) <- c(gene_set, "FID", "IID")
      PRS_inf_test <- merge(PRS_inf_test, test, by=c("IID", "FID"))

      training <- read.table(paste0(out_LDpred2_on_set, "/", "pred_inf_", training_prefix_rsID, "_", gene_set), header=T)
      colnames(training) <- c(gene_set, "FID", "IID")
      PRS_inf_training <- merge(PRS_inf_training, training, by=c("IID", "FID"))
    }
  }

  ##### Grid
  if (model_type == "GRID model") {
    PRS_grid_test <- PCs_test[,c("FID", "IID")]
    PRS_grid_training <- PCs_training[,c("FID","IID")]

    for (gene_set in gene_sets) {
      test <- read.table(paste0(out_LDpred2_on_set, "/", "pred_grid_", test_prefix_rsID, "_", gene_set), header=T)
      colnames(test) <- c(gene_set, "FID", "IID")
      PRS_grid_test <- merge(PRS_grid_test, test, by=c("IID", "FID"))

      training <- read.table(paste0(out_LDpred2_on_set, "/", "pred_grid_", training_prefix_rsID, "_", gene_set), header=T)
      colnames(training) <- c(gene_set, "FID", "IID")
      PRS_grid_training <- merge(PRS_grid_training, training, by=c("IID", "FID"))
    }
  }

  ##### Auto grid
  if (model_type == "AUTO GRID model") {
    PRS_auto_grid_test <- PCs_test[,c("FID", "IID")]
    PRS_auto_grid_training <- PCs_training[,c("FID","IID")]

    for (gene_set in gene_sets) {
      test <- read.table(paste0(out_LDpred2_on_set, "/", "pred_auto_grid_", test_prefix_rsID, "_", gene_set), header=T)
      colnames(test) <- c(gene_set, "FID", "IID")
      PRS_auto_grid_test <- merge(PRS_auto_grid_test, test, by=c("IID", "FID"))

      training <- read.table(paste0(out_LDpred2_on_set, "/", "pred_auto_grid_", training_prefix_rsID, "_", gene_set), header=T)
      colnames(training) <- c(gene_set, "FID", "IID")
      PRS_auto_grid_training <- merge(PRS_auto_grid_training, training, by=c("IID", "FID"))
    }
  }

  ##### Auto
  if (model_type == "AUTO model") {
    PRS_auto_test <- PCs_test[,c("FID", "IID")]
    PRS_auto_training <- PCs_training[,c("FID","IID")]

    for (gene_set in gene_sets) {
      test <- read.table(paste0(out_LDpred2_on_set, "/", "pred_auto_", test_prefix_rsID, "_", gene_set), header=T)
      colnames(test) <- c(gene_set, "FID", "IID")
      PRS_auto_test <- merge(PRS_auto_test, test, by=c("IID", "FID"))

      training <- read.table(paste0(out_LDpred2_on_set, "/", "pred_auto_", training_prefix_rsID, "_", gene_set), header=T)
      colnames(training) <- c(gene_set, "FID", "IID")
      PRS_auto_training <- merge(PRS_auto_training, training, by=c("IID", "FID"))
    }
  }

  #### Merge files with PCs

  ##### Inf
  if (model_type == "INF model") {
    PRS_inf_test_PCs <- merge(PRS_inf_test, PCs_test, by=c("IID", "FID"))
    PRS_inf_training_PCs <- merge(PRS_inf_training, PCs_training, by=c("IID", "FID"))
  }

  ##### Grid
  if (model_type == "GRID model") {
    PRS_grid_test_PCs <- merge(PRS_grid_test, PCs_test, by=c("IID", "FID"))
    PRS_grid_training_PCs <- merge(PRS_grid_training, PCs_training, by=c("IID", "FID"))
  }

  ##### Auto grid
  if (model_type == "AUTO GRID model") {
    PRS_auto_grid_test_PCs <- merge(PRS_auto_grid_test, PCs_test, by=c("IID", "FID"))
    PRS_auto_grid_training_PCs <- merge(PRS_auto_grid_training, PCs_training, by=c("IID", "FID"))
  }

  ##### Auto
  if (model_type == "AUTO model") {
    PRS_auto_test_PCs <- merge(PRS_auto_test, PCs_test, by=c("IID", "FID"))
    PRS_auto_training_PCs <- merge(PRS_auto_training, PCs_training, by=c("IID", "FID"))
  }


  ##############################
  #######PC correction

  #### Merge test and training for PC correction

  ##### Inf
  if (model_type == "INF model") {
    PRS_inf_test_PCs$Group <- "test"
    PRS_inf_training_PCs$Group <- "training"
    merge_inf <- rbind(PRS_inf_test_PCs, PRS_inf_training_PCs)
    #merge_inf <- merge_inf[!is.na(merge_inf$pred_inf), ]
    na_col_inf <- colSums(is.na(merge_inf))
    all_na_cols_inf <- which(na_col_inf == nrow(merge_inf))
    if(length(all_na_cols_inf) > 0) {
      merge_inf <- merge_inf[, -all_na_cols_inf]
    } else {
      # If all_na_cols is empty, just use the original dataframe
      merge_inf <- merge_inf
    }
    na_rows_inf <- apply(merge_inf, 1, function(row) any(is.na(row)))
    merge_inf <- na.omit(merge_inf)
  }

  ##### Grid
  if (model_type == "GRID model") {
    PRS_grid_test_PCs$Group <- "test"
    PRS_grid_training_PCs$Group <- "training"
    merge_grid <- rbind(PRS_grid_test_PCs, PRS_grid_training_PCs)
    na_col_grid <- colSums(is.na(merge_grid))
    all_na_cols_grid <- which(na_col_grid == nrow(merge_grid))
    if(length(all_na_cols_grid) > 0) {
      merge_grid <- merge_grid[, -all_na_cols_grid]
    } else {
      # If all_na_cols is empty, just use the original dataframe
      merge_grid <- merge_grid
    }
    na_rows_grid <- apply(merge_grid, 1, function(row) any(is.na(row)))
    merge_grid <- merge_grid[!na_rows_grid, ]
  }

  ##### Auto grid
  if (model_type == "AUTO GRID model") {
    PRS_auto_grid_test_PCs$Group <- "test"
    PRS_auto_grid_training_PCs$Group <- "training"
    merge_auto_grid <- rbind(PRS_auto_grid_test_PCs, PRS_auto_grid_training_PCs)
    na_col_auto_grid <- colSums(is.na(merge_auto_grid))
    all_na_cols_auto_grid <- which(na_col_auto_grid == nrow(merge_auto_grid))
    if(length(all_na_cols_auto_grid) > 0) {
      merge_auto_grid <- merge_auto_grid[, -all_na_cols_auto_grid]
    } else {
      # If all_na_cols is empty, just use the original dataframe
      merge_auto_grid <- merge_auto_grid
    }
    na_rows_auto_grid <- apply(merge_auto_grid, 1, function(row) any(is.na(row)))
    merge_auto_grid <- merge_auto_grid[!na_rows_auto_grid, ]
  }

  ##### Auto
  if (model_type == "AUTO model") {
    PRS_auto_test_PCs$Group <- "test"
    PRS_auto_training_PCs$Group <- "training"
    merge_auto <- rbind(PRS_auto_test_PCs, PRS_auto_training_PCs)
    #merge_auto <- merge_auto[!is.na(merge_auto$pred_auto),]
    na_col_auto <- colSums(is.na(merge_auto))
    all_na_cols_auto <- which(na_col_auto == nrow(merge_auto))
    if(length(all_na_cols_auto) > 0) {
      merge_auto <- merge_auto[, -all_na_cols_auto]
    } else {
      # If all_na_cols is empty, just use the original dataframe
      merge_auto <- merge_auto
    }
    na_rows_auto <- apply(merge_auto, 1, function(row) any(is.na(row)))
    merge_auto <- na.omit(merge_auto)
  }

  ####PC correction

  ##### Inf
  if (model_type == "INF model") {
    scores_inf <- select(merge_inf, -IID, -FID, -Group, -colnames(PC_only))
    scores_inf_2 <- colnames(scores_inf)
    scores_inf_res <- paste0(scores_inf_2, "_res")
    for (i in 1:length(scores_inf_2)){
      score <- scores_inf_2[i]
      print(score)
      score_res <- scores_inf_res[i]
      print(score_res)
      merge_inf <- PC_correction(data=merge_inf, score=score, score_res=score_res)
    }
  }

  ##### Grid
  if (model_type == "GRID model") {
    scores_grid <- select(merge_grid, -IID, -FID, -Group, -colnames(PC_only))
    scores_grid_2 <- colnames(scores_grid)
    scores_grid_res <- paste0(scores_grid_2, "_res")
    for (i in 1:length(scores_grid_2)){
      score <- scores_grid_2[i]
      score_res <- scores_grid_res[i]
      merge_grid <- PC_correction(data=merge_grid, score=score, score_res=score_res)
    }
  }

  ##### Auto grid
  if (model_type == "AUTO GRID model") {
    scores_auto_grid <- select(merge_auto_grid, -IID, -FID, -Group, -colnames(PC_only))
    scores_auto_grid_2 <- colnames(scores_auto_grid)
    scores_auto_grid_res <- paste0(scores_auto_grid_2, "_res")
    for (i in 1:length(scores_auto_grid_2)){
      score <- scores_auto_grid_2[i]
      score_res <- scores_auto_grid_res[i]
      merge_auto_grid <- PC_correction(data=merge_auto_grid, score=score, score_res=score_res)
    }
  }

  ##### Auto
  if (model_type == "AUTO model") {
    scores_auto <- select(merge_auto, -IID, -FID, -Group, -colnames(PC_only))
    scores_auto_2 <- colnames(scores_auto)
    scores_auto_res <- paste0(scores_auto_2, "_res")
    for (i in 1:length(scores_auto_2)){
      score <- scores_auto_2[i]
      print(score)
      score_res <- scores_auto_res[i]
      print(score_res)
      merge_auto <- PC_correction(data=merge_auto, score=score, score_res=score_res)
    }
  }

  #### Split again to then do standardization

  ##### Inf
  if (model_type == "INF model") {
    PRS_inf_test_PCs <- merge_inf[(merge_inf$Group=="test"),]
    PRS_inf_test_PCs <- PRS_inf_test_PCs[, !colnames(PRS_inf_test_PCs) %in% "Group"]
    PRS_inf_training_PCs <- merge_inf[(merge_inf$Group=="training"),]
    PRS_inf_training_PCs <- PRS_inf_training_PCs[, !colnames(PRS_inf_training_PCs) %in% "Group"]
  }

  ##### Grid
  if (model_type == "GRID model") {
    PRS_grid_test_PCs <- merge_grid[(merge_grid$Group=="test"),]
    PRS_grid_test_PCs <- PRS_grid_test_PCs[, !colnames(PRS_grid_test_PCs) %in% "Group"]
    PRS_grid_training_PCs <- merge_grid[(merge_grid$Group=="training"),]
    PRS_grid_training_PCs <- PRS_grid_training_PCs[, !colnames(PRS_grid_training_PCs) %in% "Group"]
  }

  ##### Auto Grid
  if (model_type == "AUTO GRID model") {
    PRS_auto_grid_test_PCs <- merge_auto_grid[(merge_auto_grid$Group=="test"),]
    PRS_auto_grid_test_PCs <- PRS_auto_grid_test_PCs[, !colnames(PRS_auto_grid_test_PCs) %in% "Group"]
    PRS_auto_grid_training_PCs <- merge_auto_grid[(merge_auto_grid$Group=="training"),]
    PRS_auto_grid_training_PCs <- PRS_auto_grid_training_PCs[, !colnames(PRS_auto_grid_training_PCs) %in% "Group"]
  }

  ##### Auto
  if (model_type == "AUTO model") {
    PRS_auto_test_PCs <- merge_auto[(merge_auto$Group=="test"),]
    PRS_auto_test_PCs <- PRS_auto_test_PCs[, !colnames(PRS_auto_test_PCs) %in% "Group"]
    PRS_auto_training_PCs <- merge_auto[(merge_auto$Group=="training"),]
    PRS_auto_training_PCs <- PRS_auto_training_PCs[, !colnames(PRS_auto_training_PCs) %in% "Group"]
  }


  ##############################
  #######Standardization

  #### Standardize

  ##### Inf
  if (model_type == "INF model") {
    sc_inf <- c(scores_inf_2, scores_inf_res)
    sc2_inf <- paste0(sc_inf, "_sc")
    for (i in 1:length(sc_inf)){
      score <- sc_inf[i]
      score_sc <- sc2_inf[i]
      PRS_inf_training_PCs <- standardize(data_training=PRS_inf_training_PCs, data_test = PRS_inf_training_PCs, score=score, score_sc=score_sc)
      PRS_inf_test_PCs <- standardize(data_training=PRS_inf_training_PCs, data_test = PRS_inf_test_PCs, score=score, score_sc=score_sc)
    }
  }

  ##### Grid
  if (model_type == "GRID model") {
    sc_grid <- c(scores_grid_2, scores_grid_res)
    sc2_grid <- paste0(sc_grid, "_sc")
    for (i in 1:length(sc_grid)){
      score <- sc_grid[i]
      score_sc <- sc2_grid[i]
      PRS_grid_training_PCs <- standardize(data_training=PRS_grid_training_PCs, data_test = PRS_grid_training_PCs, score=score, score_sc=score_sc)
      PRS_grid_test_PCs <- standardize(data_training=PRS_grid_training_PCs, data_test = PRS_grid_test_PCs, score=score, score_sc=score_sc)
    }
  }

  ##### Auto grid
  if (model_type == "AUTO GRID model") {
    sc_auto_grid <- c(scores_auto_grid_2, scores_auto_grid_res)
    sc2_auto_grid <- paste0(sc_auto_grid, "_sc")
    for (i in 1:length(sc_auto_grid)){
      score <- sc_auto_grid[i]
      score_sc <- sc2_auto_grid[i]
      PRS_auto_grid_training_PCs <- standardize(data_training=PRS_auto_grid_training_PCs, data_test = PRS_auto_grid_training_PCs, score=score, score_sc=score_sc)
      PRS_auto_grid_test_PCs <- standardize(data_training=PRS_auto_grid_training_PCs, data_test = PRS_auto_grid_test_PCs, score=score, score_sc=score_sc)
    }
  }

  ##### Auto
  if (model_type == "AUTO model") {
    sc_auto <- c(scores_auto_2, scores_auto_res)
    sc2_auto <- paste0(sc_auto, "_sc")
    for (i in 1:length(sc_auto)){
      score <- sc_auto[i]
      score_sc <- sc2_auto[i]
      PRS_auto_training_PCs <- standardize(data_training=PRS_auto_training_PCs, data_test = PRS_auto_training_PCs, score=score, score_sc=score_sc)
      PRS_auto_test_PCs <- standardize(data_training=PRS_auto_training_PCs, data_test = PRS_auto_test_PCs, score=score, score_sc=score_sc)
    }
  }

  #### Get only scores to write to a file

  ##### Inf
  if (model_type == "INF model") {
    PRS_inf_test_scores <- colnames(PRS_inf_test_PCs)[!(colnames(PRS_inf_test_PCs) %in% colnames(PC_only))]
    PRS_inf_test <- PRS_inf_test_PCs[, PRS_inf_test_scores]
    PRS_inf_training_scores <- colnames(PRS_inf_training_PCs)[!(colnames(PRS_inf_training_PCs) %in% colnames(PC_only))]
    PRS_inf_training <- PRS_inf_training_PCs[, PRS_inf_training_scores]
    write.table(PRS_inf_test, paste0(out_LDpred2_on_set, "/", test_prefix_rsID, "_scaled_scores_inf"), row.names = F, quote = F, sep="\t")
    write.table(PRS_inf_training, paste0(out_LDpred2_on_set, "/", training_prefix_rsID, "_scaled_scores_inf"), row.names = F, quote = F, sep="\t")
  }

  ##### Grid
  if (model_type == "GRID model") {
    PRS_grid_test_scores <- colnames(PRS_grid_test_PCs)[!(colnames(PRS_grid_test_PCs) %in% colnames(PC_only))]
    PRS_grid_test <- PRS_grid_test_PCs[, PRS_grid_test_scores]
    PRS_grid_training_scores <- colnames(PRS_grid_training_PCs)[!(colnames(PRS_grid_training_PCs) %in% colnames(PC_only))]
    PRS_grid_training <- PRS_grid_training_PCs[, PRS_grid_training_scores]
    write.table(PRS_grid_test, paste0(out_LDpred2_on_set, "/", test_prefix_rsID, "_scaled_scores_grid"), row.names = F, quote = F, sep="\t")
    write.table(PRS_grid_training, paste0(out_LDpred2_on_set, "/", training_prefix_rsID, "_scaled_scores_grid"), row.names = F, quote = F, sep="\t")
  }

  ##### Auto grid
  if (model_type == "AUTO GRID model") {
    PRS_auto_grid_test_scores <- colnames(PRS_auto_grid_test_PCs)[!(colnames(PRS_auto_grid_test_PCs) %in% colnames(PC_only))]
    PRS_auto_grid_test <- PRS_auto_grid_test_PCs[, PRS_auto_grid_test_scores]
    PRS_auto_grid_training_scores <- colnames(PRS_auto_grid_training_PCs)[!(colnames(PRS_auto_grid_training_PCs) %in% colnames(PC_only))]
    PRS_auto_grid_training <- PRS_auto_grid_training_PCs[, PRS_auto_grid_training_scores]
    write.table(PRS_auto_grid_test, paste0(out_LDpred2_on_set, "/", test_prefix_rsID, "_scaled_scores_auto_grid"), row.names = F, quote = F, sep="\t")
    write.table(PRS_auto_grid_training, paste0(out_LDpred2_on_set, "/", training_prefix_rsID, "_scaled_scores_auto_grid"), row.names = F, quote = F, sep="\t")
  }

  ##### Auto
    if (model_type == "AUTO model") {
    PRS_auto_test_scores <- colnames(PRS_auto_test_PCs)[!(colnames(PRS_auto_test_PCs) %in% colnames(PC_only))]
    PRS_auto_test <- PRS_auto_test_PCs[, PRS_auto_test_scores]
    PRS_auto_training_scores <- colnames(PRS_auto_training_PCs)[!(colnames(PRS_auto_training_PCs) %in% colnames(PC_only))]
    PRS_auto_training <- PRS_auto_training_PCs[, PRS_auto_training_scores]
    write.table(PRS_auto_test, paste0(out_LDpred2_on_set, "/", test_prefix_rsID, "_scaled_scores_auto"), row.names = F, quote = F, sep="\t")
    write.table(PRS_auto_training, paste0(out_LDpred2_on_set, "/", training_prefix_rsID, "_scaled_scores_auto"), row.names = F, quote = F, sep="\t")
  }


  ##############################
  #######Linear regression

  #### Prepare files for regression

  ##### Inf
  if (model_type == "INF model") {
    scores_regr_inf <- colnames(PRS_inf_test)[grep("_res_sc$", colnames(PRS_inf_test))]
    columns_inf <- c("IID", "FID", scores_regr_inf)
    PRS_inf_for_regression <- PRS_inf_test[, columns_inf]
    PRS_inf_for_regression <- merge(PRS_inf_for_regression, pheno, by=c("IID", "FID"))
    colnames(PRS_inf_for_regression) <- c("IID", "FID",scores_regr_inf, "PHENO")
    PRS_inf_for_regression$PHENO <- as.numeric(PRS_inf_for_regression$PHENO)
  }

  ##### Grid
  if (model_type == "GRID model") {
    scores_regr_grid <- colnames(PRS_grid_test)[grep("_res_sc$", colnames(PRS_grid_test))]
    columns_grid <- c("IID", "FID", scores_regr_grid)
    PRS_grid_for_regression <- PRS_grid_test[, columns_grid]
    PRS_grid_for_regression <- merge(PRS_grid_for_regression, pheno, by=c("IID", "FID"))
    colnames(PRS_grid_for_regression) <- c("IID", "FID",scores_regr_grid, "PHENO")
    PRS_grid_for_regression$PHENO <- as.numeric(PRS_grid_for_regression$PHENO)
  }

  ##### Auto grid
  if (model_type == "AUTO GRID model") {
    scores_regr_auto_grid <- colnames(PRS_auto_grid_test)[grep("_res_sc$", colnames(PRS_auto_grid_test))]
    columns_auto_grid <- c("IID", "FID", scores_regr_auto_grid)
    PRS_auto_grid_for_regression <- PRS_auto_grid_test[, columns_auto_grid]
    PRS_auto_grid_for_regression <- merge(PRS_auto_grid_for_regression, pheno, by=c("IID", "FID"))
    colnames(PRS_auto_grid_for_regression) <- c("IID", "FID",scores_regr_auto_grid, "PHENO")
    PRS_auto_grid_for_regression$PHENO <- as.numeric(PRS_auto_grid_for_regression$PHENO)
  }

  ##### Auto
  if (model_type == "AUTO model") {
    scores_regr_auto <- colnames(PRS_auto_test)[grep("_res_sc$", colnames(PRS_auto_test))]
    columns_auto <- c("IID", "FID", scores_regr_auto)
    PRS_auto_for_regression <- PRS_auto_test[, columns_auto]
    PRS_auto_for_regression <- merge(PRS_auto_for_regression, pheno, by=c("IID", "FID"))
    colnames(PRS_auto_for_regression) <- c("IID", "FID",scores_regr_auto, "PHENO")
    PRS_auto_for_regression$PHENO <- as.numeric(PRS_auto_for_regression$PHENO)
  }

  #### Perform linear regression

  if (cross_validation=="FALSE"){

    ##### Inf
    if (model_type == "INF model") {
      Regression_results_inf <- data.frame(matrix(ncol=10, nrow=0))
      for (i in 1:length(scores_regr_inf)) {
        b <- scores_regr_inf[i]
        thresh <- paste(b)  
        tip <- paste("PHENO", "~", b,sep="")
        alpha <- lm(tip, data = PRS_inf_for_regression)
        R2_full <- summary(alpha)$adj.r.squared
        p <- coef(summary(alpha))[2,4]
        beta <- coef(summary(alpha))[2,1]
        SE <- coef(summary(alpha))[2,2]
        Regression_results_inf <- rbind(Regression_results_inf, c("inf", thresh,beta,SE,p,R2_full, "LDpred2"))
        colnames(Regression_results_inf) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "Tool")
      }
      Regression_results_inf$R2 <- as.numeric(Regression_results_inf$R2)
      Regression_results_inf$P_value <- as.numeric(Regression_results_inf$P_value)
      Regression_results_inf$P_value <- ifelse(Regression_results_inf$P_value==0, 3.4e-314, Regression_results_inf$P_value)

      write.table(Regression_results_inf, paste0(out_LDpred2_on_set, "/Regression_results_LDpred2_inf"), row.names=F, quote=F, sep="\t")
    }

    ##### Grid
    if (model_type == "GRID model") {
      Regression_results_grid <- data.frame(matrix(ncol=10, nrow=0))
      for (i in 1:length(scores_regr_grid)) {
        b <- scores_regr_grid[i]
        thresh <- paste(b)  
        tip <- paste("PHENO", "~", b,sep="")
        alpha <- lm(tip, data = PRS_grid_for_regression)
        R2_full <- summary(alpha)$adj.r.squared
        p <- coef(summary(alpha))[2,4]
        beta <- coef(summary(alpha))[2,1]
        SE <- coef(summary(alpha))[2,2]
        Regression_results_grid <- rbind(Regression_results_grid, c("grid", thresh,beta,SE,p,R2_full, "LDpred2"))
        colnames(Regression_results_grid) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "Tool")
      }

      Regression_results_grid$R2 <- as.numeric(Regression_results_grid$R2)
      Regression_results_grid$P_value <- as.numeric(Regression_results_grid$P_value)
      Regression_results_grid$P_value <- ifelse(Regression_results_grid$P_value==0, 3.4e-314, Regression_results_grid$P_value)

      write.table(Regression_results_grid, paste0(out_LDpred2_on_set, "/Regression_results_LDpred2_grid"), row.names=F, quote=F, sep="\t")
    }

    ##### Auto grid
    if (model_type == "AUTO GRID model") {
      Regression_results_auto_grid <- data.frame(matrix(ncol=10, nrow=0))
      for (i in 1:length(scores_regr_auto_grid)) {
        b <- scores_regr_auto_grid[i]
        thresh <- paste(b)  
        tip <- paste("PHENO", "~", b,sep="")
        alpha <- lm(tip, data = PRS_auto_grid_for_regression)
        R2_full <- summary(alpha)$adj.r.squared
        p <- coef(summary(alpha))[2,4]
        beta <- coef(summary(alpha))[2,1]
        SE <- coef(summary(alpha))[2,2]
        Regression_results_auto_grid <- rbind(Regression_results_auto_grid, c("auto_grid", thresh,beta,SE,p,R2_full, "LDpred2"))
        colnames(Regression_results_auto_grid) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "Tool")
      }

      Regression_results_auto_grid$R2 <- as.numeric(Regression_results_auto_grid$R2)
      Regression_results_auto_grid$P_value <- as.numeric(Regression_results_auto_grid$P_value)
      Regression_results_auto_grid$P_value <- ifelse(Regression_results_auto_grid$P_value==0, 3.4e-314, Regression_results_auto_grid$P_value)

      write.table(Regression_results_auto_grid, paste0(out_LDpred2_on_set, "/Regression_results_LDpred2_auto_grid"), row.names=F, quote=F, sep="\t")
    }

    ##### Auto
    if (model_type == "AUTO model") {
      Regression_results_auto <- data.frame(matrix(ncol=10, nrow=0))
      for (i in 1:length(scores_regr_auto)) {
        b <- scores_regr_auto[i]
        thresh <- paste(b)  
        tip <- paste("PHENO", "~", b,sep="")
        alpha <- lm(tip, data = PRS_auto_for_regression)
        R2_full <- summary(alpha)$adj.r.squared
        p <- coef(summary(alpha))[2,4]
        beta <- coef(summary(alpha))[2,1]
        SE <- coef(summary(alpha))[2,2]
        Regression_results_auto <- rbind(Regression_results_auto, c("auto", thresh,beta,SE,p,R2_full, "LDpred2"))
        colnames(Regression_results_auto) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "Tool")
      }

      Regression_results_auto$R2 <- as.numeric(Regression_results_auto$R2)
      Regression_results_auto$P_value <- as.numeric(Regression_results_auto$P_value)
      Regression_results_auto$P_value <- ifelse(Regression_results_auto$P_value==0, 3.4e-314, Regression_results_auto$P_value)

      write.table(Regression_results_auto, paste0(out_LDpred2_on_set, "/Regression_results_LDpred2_auto"), row.names=F, quote=F, sep="\t")
    }

  } else if (cross_validation=="TRUE"){

    if (model_type == "INF model") {
      PRS_for_regression = PRS_inf_for_regression
      scores_regr = scores_regr_inf
      model_name = "inf"
    }
    if (model_type == "GRID model") {
      PRS_for_regression = PRS_grid_for_regression
      scores_regr = scores_regr_grid
      model_name = "grid"
    }
    if (model_type == "AUTO GRID model") {
      PRS_for_regression = PRS_auto_grid_for_regression
      scores_regr = scores_regr_auto_grid
      model_name = "auto_grid"
    }
    if (model_type == "AUTO model") {
      PRS_for_regression = PRS_auto_for_regression
      scores_regr = scores_regr_auto
      model_name = "auto"
    }

    Regression_results_LDpred2 <- data.frame(matrix(ncol=11, nrow=0))
    train_control <- trainControl(method = "cv", number = 10)

    for (i in 1:length(scores_regr)) {
      b <- scores_regr[i]
      thresh <- paste(b)  
      tip <- paste("PHENO", "~", b,sep="")
      alpha <- lm(tip, data = PRS_for_regression)
      R2 <- summary(alpha)$adj.r.squared
      p <- coef(summary(alpha))[2,4]
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      model <- train(as.formula(tip), data=PRS_for_regression, method="lm", trControl=train_control, metric="Rsquared")
      R2_CV <- as.numeric(model$results$Rsquared[1])
      Regression_results_LDpred2 <- rbind(Regression_results_LDpred2, c("LDpred2",model_name, thresh,beta,SE,p,R2,R2_CV))
      colnames(Regression_results_LDpred2) <- c("Tool", "Model", "Parameters", "Beta", "SE", "P_value", "R2", "R2_CV")
    }

    Regression_results_LDpred2$R2 <- as.numeric(Regression_results_LDpred2$R2)
    Regression_results_LDpred2$R2_CV <- as.numeric(Regression_results_LDpred2$R2_CV)
    Regression_results_LDpred2$P_value <- as.numeric(Regression_results_LDpred2$P_value)
    Regression_results_LDpred2$P_value <- ifelse(Regression_results_LDpred2$P_value==0, 3.4e-314, Regression_results_LDpred2$P_value)

    write.table(Regression_results_LDpred2, paste0(out_LDpred2_on_set, "/Regression_results_LDpred2_", model_name), row.names=F, quote=F, sep="\t")
    assign(paste0("Regression_results_", model_name), Regression_results_LDpred2)

  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ##############################
  #######Plot

  #### Multiset bar plot
  if (model_type == "INF model") {
    Regression_results_LDpred2 = Regression_results_inf
    model_name = "inf"
  }
  if (model_type == "GRID model") {
    Regression_results_LDpred2 = Regression_results_grid
    model_name = "grid"
  }
  if (model_type == "AUTO GRID model") {
    Regression_results_LDpred2 = Regression_results_auto_grid
    model_name = "auto_grid"
  }
  if (model_type == "AUTO model") {
    Regression_results_LDpred2 = Regression_results_auto
    model_name = "auto"
  }

  if (cross_validation=="FALSE"){
    prs_top <- Regression_results_LDpred2 %>%
    mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
    mutate(Parameters = gsub("_res_sc", "", Parameters)) %>%
    arrange(desc(R2)) %>%                   # Order by variance explained
    head(7)   
  } else if (cross_validation=="TRUE"){    
    prs_top <- Regression_results_LDpred2 %>%
    mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
    mutate(Parameters = gsub("_res_sc", "", Parameters)) %>%
    arrange(desc(R2_CV)) %>%                   # Order by variance explained
    head(7)                                 # Top 7
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ## Make multiset bar plot
  if (cross_validation=="FALSE"){
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2", logP="log", out=paste0(out_LDpred2_on_set, "/multiset_bar_plot_R2_LDpred2_", model_name, ".svg"))
  } else if (cross_validation=="TRUE"){      
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2_CV", logP="log", out=paste0(out_LDpred2_on_set, "/multiset_bar_plot_R2_LDpred2_", model_name, ".svg"))
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  #### Correlation plot
  if (model_type == "INF model") {PRS_LDpred2_for_regression = PRS_inf_for_regression}
  if (model_type == "GRID model") {PRS_LDpred2_for_regression = PRS_grid_for_regression}
  if (model_type == "AUTO GRID model") {PRS_LDpred2_for_regression = PRS_auto_grid_for_regression}
  if (model_type == "AUTO model") {PRS_LDpred2_for_regression = PRS_auto_for_regression}
  # Extract the columns data
  pathway_columns <- grep("_res_sc$", colnames(PRS_LDpred2_for_regression), value = TRUE)
  # Rename pathways: Remove "_res_sc" and replace "_" with " "
  pathway_labels <- gsub("_", " ", gsub("_res_sc", "", pathway_columns))
  # Make plot
  correlation_plot(PRS_LDpred2_for_regression, pathway_columns, pathway_labels, paste0(out_LDpred2_on_set, "/correlation_LDpred2.svg"))

  ##############################
  #######Covariates

  if (cov_analysis == TRUE) {
    
    if (model_type == "INF model") {scores_regr = scores_regr_inf}
    if (model_type == "GRID model") {scores_regr = scores_regr_grid}
    if (model_type == "AUTO GRID model") {scores_regr = scores_regr_auto_grid}
    if (model_type == "AUTO model") {scores_regr = scores_regr_auto}

    Regression_covariates_LDpred2 <- data.frame(matrix(ncol=6, nrow=0))
    for (i in 1:length(scores_regr)) {
      b <- scores_regr[i]
      thresh <- paste(b) 
      for (covs in covariate_list){
        PRS_test_LDpred2_cov <- merge(PRS_LDpred2_for_regression, cov_data, by=c("FID", "IID"))
        full_form <- paste("PHENO", "~", b, "+", covs, sep="")
        null_form <- paste("PHENO", "~", covs, sep="")
        full_lm <- lm(full_form, data = PRS_test_LDpred2_cov)
        null_lm <- lm(null_form, data = PRS_test_LDpred2_cov)
        R2_full <- summary(full_lm)$adj.r.squared
        R2_null <- summary(null_lm)$adj.r.squared
        R2_PRS <- R2_full - R2_null
        Regression_covariates_LDpred2 <- rbind(Regression_covariates_LDpred2, c("LDpred2", thresh, R2_full, R2_null, R2_PRS, covs))
        colnames(Regression_covariates_LDpred2) <- c("Tool", "Parameters", "R2_PRS_and_cov" , "R2_cov_only", "R2_PRS_only", "included_covariates")
      }
    }
    write.table(Regression_covariates_LDpred2, paste0(out_LDpred2_on_set, "/Regression_results_LDpred2_covariates"), row.names=F, quote=F, sep="\t")
  }

  cat("Done getting best PRS for LDpred2.... lassosum2....\n")

} else {
  compare_tools = FALSE
  message("Skipping get best set for LDpred2.... lassosum2....")
}

##################
####Lassosum-2####
##################

if (length(list.files(out_lasso2_on_set)) > 0) {

  ##############################
  #######Load scores

  PRS_lasso2_test <- PCs_test[,c("FID", "IID")]
  PRS_lasso2_training <- PCs_training[,c("FID","IID")]

  for (gene_set in gene_sets) {
    test <- read.table(paste0(out_lasso2_on_set, "/", "pred_lasso2_", test_prefix_rsID, "_", gene_set), header=T)
    colnames(test) <- c(gene_set, "FID", "IID")
    PRS_lasso2_test <- merge(PRS_lasso2_test, test, by=c("IID", "FID"))
    
    training <- read.table(paste0(out_lasso2_on_set, "/", "pred_lasso2_", training_prefix_rsID, "_", gene_set), header=T)
    colnames(training) <- c(gene_set, "FID", "IID")
    PRS_lasso2_training <- merge(PRS_lasso2_training, training, by=c("IID", "FID"))
  }

  #### Merge files with PCs
  PRS_lasso2_test_PCs <- merge(PRS_lasso2_test, PCs_test, by=c("IID", "FID"))
  PRS_lasso2_training_PCs <- merge(PRS_lasso2_training, PCs_training, by=c("IID", "FID"))


  ##############################
  #######PC correction

  #### Merge test and training for PC correction

  PRS_lasso2_test_PCs$Group <- "test"
  PRS_lasso2_training_PCs$Group <- "training"
  merge_lasso2 <- rbind(PRS_lasso2_test_PCs, PRS_lasso2_training_PCs)
  na_col_lasso2 <- colSums(is.na(merge_lasso2))
  all_na_cols_lasso2 <- which(na_col_lasso2 == nrow(merge_lasso2))
  if(length(all_na_cols_lasso2) > 0) {
    merge_lasso2 <- merge_lasso2[, -all_na_cols_lasso2]
  } else {
    # If all_na_cols is empty, just use the original dataframe
    merge_lasso2 <- merge_lasso2
  }
  na_rows_lasso2 <- apply(merge_lasso2, 1, function(row) any(is.na(row)))
  merge_lasso2 <- na.omit(merge_lasso2)

  ####PC correction

  scores_lasso2 <- select(merge_lasso2, -IID, -FID, -Group, -colnames(PC_only))
  scores_lasso2_2 <- colnames(scores_lasso2)
  scores_lasso2_res <- paste0(scores_lasso2_2, "_res")
  for (i in 1:length(scores_lasso2_2)){
    score <- scores_lasso2_2[i]
    score_res <- scores_lasso2_res[i]
    merge_lasso2 <- PC_correction(data=merge_lasso2, score=score, score_res=score_res)
  }

  #### Split again to then do standardization

  PRS_lasso2_test_PCs <- merge_lasso2[(merge_lasso2$Group=="test"),]
  PRS_lasso2_test_PCs <- PRS_lasso2_test_PCs[, !colnames(PRS_lasso2_test_PCs) %in% "Group"]
  PRS_lasso2_training_PCs <- merge_lasso2[(merge_lasso2$Group=="training"),]
  PRS_lasso2_training_PCs <- PRS_lasso2_training_PCs[, !colnames(PRS_lasso2_training_PCs) %in% "Group"]

  ##############################
  #######Standardization

  #### Standardize

  sc_lasso2 <- c(scores_lasso2_2, scores_lasso2_res)
  sc2_lasso2 <- paste0(sc_lasso2, "_sc")
  for (i in 1:length(sc_lasso2)){
    score <- sc_lasso2[i]
    score_sc <- sc2_lasso2[i]
    PRS_lasso2_training_PCs <- standardize(data_training=PRS_lasso2_training_PCs, data_test = PRS_lasso2_training_PCs, score=score, score_sc=score_sc)
    PRS_lasso2_test_PCs <- standardize(data_training=PRS_lasso2_training_PCs, data_test = PRS_lasso2_test_PCs, score=score, score_sc=score_sc)
  }

  #### Get only scores to write to a file

  PRS_lasso2_test_scores <- colnames(PRS_lasso2_test_PCs)[!(colnames(PRS_lasso2_test_PCs) %in% colnames(PC_only))]
  PRS_lasso2_test <- PRS_lasso2_test_PCs[, PRS_lasso2_test_scores]
  PRS_lasso2_training_scores <- colnames(PRS_lasso2_training_PCs)[!(colnames(PRS_lasso2_training_PCs) %in% colnames(PC_only))]
  PRS_lasso2_training <- PRS_lasso2_training_PCs[, PRS_lasso2_training_scores]
  write.table(PRS_lasso2_test, paste0(out_lasso2_on_set, "/", test_prefix_rsID, "_scaled_scores_lasso2"), row.names = F, quote = F, sep="\t")
  write.table(PRS_lasso2_training, paste0(out_lasso2_on_set, "/", training_prefix_rsID, "_scaled_scores_lasso2"), row.names = F, quote = F, sep="\t")

  ##############################
  #######Linear regression

  #### Prepare files for regression
  scores_regr_lasso2 <- colnames(PRS_lasso2_test)[grep("_res_sc$", colnames(PRS_lasso2_test))]
  columns_lasso2 <- c("IID", "FID", scores_regr_lasso2)
  PRS_lasso2_for_regression <- PRS_lasso2_test[, columns_lasso2]
  PRS_lasso2_for_regression <- merge(PRS_lasso2_for_regression, pheno, by=c("IID", "FID"))
  colnames(PRS_lasso2_for_regression) <- c("IID", "FID",scores_regr_lasso2, "PHENO")
  PRS_lasso2_for_regression$PHENO <- as.numeric(PRS_lasso2_for_regression$PHENO)

  #### Perform linear regression

  if (cross_validation=="FALSE"){

    Regression_results_lasso2 <- data.frame(matrix(ncol=6, nrow=0))
    
    for (i in 1:length(scores_regr_lasso2)) {
      b <- scores_regr_lasso2[i]
      thresh <- paste(b)  
      tip <- paste("PHENO", "~", b,sep="")
      alpha <- lm(tip, data = PRS_lasso2_for_regression)
      R2_full <- summary(alpha)$adj.r.squared
      p <- coef(summary(alpha))[2,4]
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      Regression_results_lasso2 <- rbind(Regression_results_lasso2, c("lassosum2", thresh,beta,SE,p,R2_full))
      colnames(Regression_results_lasso2) <- c("Tool", "Parameters", "Beta", "SE", "P_value", "R2")
    }

    Regression_results_lasso2$R2 <- as.numeric(Regression_results_lasso2$R2)
    Regression_results_lasso2$P_value <- as.numeric(Regression_results_lasso2$P_value)
    Regression_results_lasso2$P_value <- ifelse(Regression_results_lasso2$P_value==0, 3.4e-314, Regression_results_lasso2$P_value)

    write.table(Regression_results_lasso2, paste0(out_lasso2_on_set, "/Regression_results_lassosum2"), row.names=F, quote=F, sep="\t")


  } else if (cross_validation=="TRUE"){

    Regression_results_lasso2 <- data.frame(matrix(ncol=7, nrow=0))
    train_control <- trainControl(method = "cv", number = 10)
    
    for (i in 1:length(scores_regr_lasso2)) {
      b <- scores_regr_lasso2[i]
      thresh <- paste(b)  
      tip <- paste("PHENO", "~", b,sep="")
      alpha <- lm(tip, data = PRS_lasso2_for_regression)
      R2 <- summary(alpha)$adj.r.squared
      p <- coef(summary(alpha))[2,4]
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      model <- train(as.formula(tip), data=PRS_lasso2_for_regression, method="lm", trControl=train_control, metric="Rsquared")
      R2_CV <- as.numeric(model$results$Rsquared[1])
      Regression_results_lasso2 <- rbind(Regression_results_lasso2, c("lassosum2", thresh,beta,SE,p,R2,R2_CV))
      colnames(Regression_results_lasso2) <- c("Tool", "Parameters", "Beta", "SE", "P_value", "R2", "R2_CV")
    }

    Regression_results_lasso2$R2 <- as.numeric(Regression_results_lasso2$R2)
    Regression_results_lasso2$R2_CV <- as.numeric(Regression_results_lasso2$R2_CV)
    Regression_results_lasso2$P_value <- as.numeric(Regression_results_lasso2$P_value)
    Regression_results_lasso2$P_value <- ifelse(Regression_results_lasso2$P_value==0, 3.4e-314, Regression_results_lasso2$P_value)

    write.table(Regression_results_lasso2, paste0(out_lasso2_on_set, "/Regression_results_lassosum2"), row.names=F, quote=F, sep="\t")

  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ##############################
  #######Plot

  #### Multiset bar plot
  ## Calculate -log10(P) and filter out unnecessary rows if needed
  if (cross_validation=="FALSE"){
    prs_top <- Regression_results_lasso2 %>%
      mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
      mutate(Parameters = gsub("_res_sc", "", Parameters)) %>%
      arrange(desc(R2)) %>%                   # Order by variance explained
      head(7)                                 # Top 7
  } else if (cross_validation=="TRUE"){    
    prs_top <- Regression_results_lasso2 %>%
      mutate(logP = -log10(P_value)) %>%      # Create -log10(P) column
      mutate(Parameters = gsub("_res_sc", "", Parameters)) %>%
      arrange(desc(R2_CV)) %>%                # Order by variance explained
      head(7)                                 # Top 7
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  ## Make multiset bar plot
  if (cross_validation=="FALSE"){
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2", logP="log", out=paste0(out_lasso2_on_set, "/multiset_bar_plot_R2_lasso2.svg"))
  } else if (cross_validation=="TRUE"){    
    multiset_barplot(data=prs_top, Parameters="Parameters", R2="R2_CV", logP="log", out=paste0(out_lasso2_on_set, "/multiset_bar_plot_R2_lasso2.svg"))
  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  #### Correlation plot
  # Extract the columns data
  pathway_columns <- grep("_res_sc$", colnames(PRS_lasso2_for_regression), value = TRUE)
  # Rename pathways: Remove "_res_sc" and replace "_" with " "
  pathway_labels <- gsub("_", " ", gsub("_res_sc", "", pathway_columns))
  # Make plot
  correlation_plot(PRS_lasso2_for_regression, pathway_columns, pathway_labels, paste0(out_lasso2_on_set, "/correlation_lasso2.svg"))

  ##############################
  #######Covariates

  if (cov_analysis == TRUE) {
    Regression_covariates_lasso2 <- data.frame(matrix(ncol=6, nrow=0))
    for (i in 1:length(scores_regr_lasso2)) {
      b <- scores_regr_lasso2[i]
      thresh <- paste(b) 
      for (covs in covariate_list){
        PRS_test_lasso2_cov <- merge(PRS_lasso2_for_regression, cov_data, by=c("FID", "IID"))
        full_form <- paste("PHENO", "~", b, "+", covs, sep="")
        null_form <- paste("PHENO", "~", covs, sep="")
        full_lm <- lm(full_form, data = PRS_test_lasso2_cov)
        null_lm <- lm(null_form, data = PRS_test_lasso2_cov)
        R2_full <- summary(full_lm)$adj.r.squared
        R2_null <- summary(null_lm)$adj.r.squared
        R2_PRS <- R2_full - R2_null
        Regression_covariates_lasso2 <- rbind(Regression_covariates_lasso2, c("lassosum2", thresh, R2_full, R2_null, R2_PRS, covs))
        colnames(Regression_covariates_lasso2) <- c("Tool", "Parameters", "R2_PRS_and_cov" , "R2_cov_only", "R2_PRS_only", "included_covariates")
      }
    }
    write.table(Regression_covariates_lasso2, paste0(out_lasso2_on_set, "/Regression_results_lassosum2_covariates"), row.names=F, quote=F, sep="\t")
  }


  cat("Done getting best set for lassosum2...\n")
  cat("Done getting best set per tool!\nGet best PRS per set...\n")

} else {
  compare_tools = FALSE
  message("Skipping get best set for lassosum2....")
  cat("Done getting best set per tool!\nGet best PRS per set...\n")
}


##################
####Comparison####
##################

if (compare_tools == TRUE) {

  #### Compare R2 per set

  ## Creating score files lassosum 
  data_frames_training <- list()
  data_frames_test <- list()
  for (gene_set in gene_sets) {
    file_name <- paste(out_lassosum_on_set, "/", test_prefix, "_", gene_set, "_scaled_scores", sep="")
    data <- read.table(file_name, header=T)
    data <- data[, colSums(data==0) !=(nrow(data))] #Delete columns where all values are 0
    data_frames_test[[gene_set]] <- data
    
    file_name <- paste(out_lassosum_on_set, "/", training_prefix, "_", gene_set, "_scaled_scores", sep="")
    data <- read.table(file_name, header=T)
    data <- data[, colSums(data==0) !=(nrow(data))]
    data_frames_training[[gene_set]] <- data
  }

  test_scaled_scores <- NULL  # To collect merged data
  training_scaled_scores <- NULL
  for (gene_set in gene_sets) {
    test <- data_frames_test[[gene_set]]
    training <- data_frames_training[[gene_set]]
    
    # Select FID, IID, and the _res_sc column
    res_sc_col <- grep("_res_sc$", colnames(test), value = TRUE)
    selected <- test[, c("FID", "IID", res_sc_col)]
    # Rename the _res_sc column to [gene_set]_res_sc
    colnames(selected)[3] <- paste0(gene_set, "_res_sc")
    # Merge with the main dataframe
    if (is.null(test_scaled_scores)) {
      test_scaled_scores <- selected
    } else {
      test_scaled_scores <- merge(test_scaled_scores, selected, by = c("FID", "IID"), all = TRUE)
    }

    # Select FID, IID, and the _res_sc column
    res_sc_col <- grep("_res_sc$", colnames(training), value = TRUE)
    selected <- training[, c("FID", "IID", res_sc_col)]
    # Rename the _res_sc column to [gene_set]_res_sc
    colnames(selected)[3] <- paste0(gene_set, "_res_sc")
    # Merge with the main dataframe
    if (is.null(training_scaled_scores)) {
      training_scaled_scores <- selected
    } else {
      training_scaled_scores <- merge(training_scaled_scores, selected, by = c("FID", "IID"), all = TRUE)
    }
  }
  write.table(test_scaled_scores, paste0(out_lassosum_on_set, "/", test_prefix, "_scaled_scores"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(training_scaled_scores, paste0(out_lassosum_on_set, "/", training_prefix, "_scaled_scores"), sep = "\t", row.names = FALSE, quote = FALSE)

  ## Select LDpred2 model
  if (model_type == "INF model") {
    Regression_results_LDpred2 = Regression_results_inf
    model_name = "inf"
  }
  if (model_type == "GRID model") {
    Regression_results_LDpred2 = Regression_results_grid
    model_name = "grid"
  }
  if (model_type == "AUTO GRID model") {
    Regression_results_LDpred2 = Regression_results_auto_grid
    model_name = "auto_grid"
  }
  if (model_type == "AUTO model") {
    Regression_results_LDpred2 = Regression_results_auto
    model_name = "auto"
  }

  for (gene_set in gene_sets) {

      # Select regression result for set
      Regression_best_per_tool <- data.frame()
      Regression_result_PRSet <- Regression_results_PRSet[grepl(gene_set, Regression_results_PRSet$Parameters), ]
      Regression_result_PRSice <- Regression_results_PRSice[grepl(gene_set, Regression_results_PRSice$Parameters), ]
      Regression_result_lasso <- Regression_results_lasso[grepl(gene_set, Regression_results_lasso$Parameters), ]
      Regression_result_PRScs <- Regression_results_PRS_CS[grepl(gene_set, Regression_results_PRS_CS$Parameters), ]
      Regression_result_LDpred2 <- Regression_results_LDpred2[grepl(gene_set, Regression_results_LDpred2$Parameters), ]
      Regression_result_lasso2 <- Regression_results_lasso2[grepl(gene_set, Regression_results_lasso2$Parameters), ]

      Regression_result_LDpred2 <- Regression_result_LDpred2 %>%
        select(-Model) %>%       # remove Model
        relocate(Tool, .before = everything())
      
      Regression_result_lasso <- Regression_result_lasso %>%
        select(-Set, -Lambda) %>%
        relocate(Tool, Parameters, .before = everything())

      Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_result_PRSet)
      Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_result_PRSice)
      Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_result_lasso)
      Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_result_PRScs)

      Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_result_LDpred2)
      Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_result_lasso2)

      Regression_best_per_tool$Tool <- factor(Regression_best_per_tool$Tool, levels = c("PRSet", "PRSice", "PRS-CS", "LDpred2","lassosum", "lassosum2"))
      
      if (cross_validation=="FALSE"){
        Comparison_bar_plot <- ggplot(data=Regression_best_per_tool, aes(x=Tool, y=R2, fill=Tool)) +
          geom_bar(stat="identity", position="dodge") +
            scale_fill_manual(values=c("purple4", "#CC79A7", "#56B4E9","cyan3", "#E69F00","coral"), name="Tool") +
            theme_minimal()+
            theme(legend.position = "none") +
            theme(axis.text.x = element_text(size=12, face = "bold"),
                    axis.text.y = element_text(size=12, face = "bold"),
                    axis.title.y = element_text(size=15,face="bold"),
                    axis.title.x = element_text(size=15,face="bold")) +
            xlab("") + ylab("R²\n") + ggtitle("Comparison R² per tool")
        ggsave(paste0(out_comparison, "/Bar_plot_comparison_R2_per_tool_", gene_set, ".svg"), width=7, height=5)
      
      } else if (cross_validation=="TRUE"){
        Comparison_bar_plot <- ggplot(data=Regression_best_per_tool, aes(x=Tool, y=R2_CV, fill=Tool)) +
            geom_bar(stat="identity", position="dodge") +
            scale_fill_manual(values=c("purple4", "#CC79A7", "#56B4E9","cyan3", "#E69F00","coral"), name="Tool") +
            theme_minimal()+
            theme(legend.position = "none") +
            theme(axis.text.x = element_text(size=12, face = "bold"),
                    axis.text.y = element_text(size=12, face = "bold"),
                    axis.title.y = element_text(size=15,face="bold"),
                    axis.title.x = element_text(size=15,face="bold")) +
            xlab("") + ylab("R²\n") + ggtitle("Comparison R² per tool")
        ggsave(paste0(out_comparison, "/Bar_plot_comparison_R2_per_tool_", gene_set, ".svg"), width=7, height=5)

      } else {
        "Fill in TRUE or FALSE for cross-validation"
      }
  }


  ### Compare R2 over all sets

  # Collect all results across gene sets
  All_results <- data.frame()

  for (gene_set in gene_sets) {
      
      # Select regression result for set
      Regression_result_PRSet   <- Regression_results_PRSet[grepl(gene_set, Regression_results_PRSet$Parameters), ]
      Regression_result_PRSice  <- Regression_results_PRSice[grepl(gene_set, Regression_results_PRSice$Parameters), ]
      Regression_result_lasso   <- Regression_results_lasso[grepl(gene_set, Regression_results_lasso$Parameters), ]
      Regression_result_PRScs   <- Regression_results_PRS_CS[grepl(gene_set, Regression_results_PRS_CS$Parameters), ]
      Regression_result_LDpred2 <- Regression_results_LDpred2[grepl(gene_set, Regression_results_LDpred2$Parameters), ]
      Regression_result_lasso2  <- Regression_results_lasso2[grepl(gene_set, Regression_results_lasso2$Parameters), ]

      Regression_result_LDpred2 <- Regression_result_LDpred2 %>%
        select(-Model) %>%       # remove Model
        relocate(Tool, .before = everything())
      
      Regression_result_lasso <- Regression_result_lasso %>%
        select(-Set, -Lambda) %>%
        relocate(Tool, Parameters, .before = everything())
      
      # Combine results and add gene set label
      combined <- rbind(
          Regression_result_PRSet,
          Regression_result_PRSice,
          Regression_result_lasso,
          Regression_result_PRScs,
          Regression_result_LDpred2,
          Regression_result_lasso2
      )
      combined$GeneSet <- gene_set

      All_results <- rbind(All_results, combined)
  }

  # Set factor levels
  All_results$Tool <- factor(All_results$Tool, levels = c("PRSet", "PRSice", "PRS-CS", "LDpred2", "lassosum", "lassosum2"))
  All_results$GeneSet <- factor(All_results$GeneSet, levels = gene_sets)  # preserve input order

  # Plot: Gene sets on x-axis, tool as fill

  if (cross_validation=="FALSE"){
    Comparison_combined_plot <- ggplot(All_results, aes(x = GeneSet, y = R2, fill = Tool)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
        scale_fill_manual(values = c("purple4", "#CC79A7", "#56B4E9", "cyan3", "#E69F00", "coral"), name = "Tool") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 15, face = "bold"),
            axis.title.x = element_text(size = 15, face = "bold"),
            legend.position = "right"
        ) +
        xlab("Gene Set") + 
        ylab("Nagelkerke R²")

  } else if (cross_validation=="TRUE"){
    Comparison_combined_plot <- ggplot(All_results, aes(x = GeneSet, y = R2_CV, fill = Tool)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
        scale_fill_manual(values = c("purple4", "#CC79A7", "#56B4E9", "cyan3", "#E69F00", "coral"), name = "Tool") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 15, face = "bold"),
            axis.title.x = element_text(size = 15, face = "bold"),
            legend.position = "right"
        ) +
        xlab("Gene Set") + 
        ylab("Nagelkerke R²")

  } else {
    "Fill in TRUE or FALSE for cross-validation"
  }

  # Save
  ggsave(paste0(out_comparison, "/Bar_plot_comparison_R2_sets.svg"), plot = Comparison_combined_plot, width = 10, height = 6)

  cat("Done comparing tools!\nGene-set extension completed!\n")

} else {
  message("Skipping comparing tools, please run all tools if you want a comparison.\nGene-set extension completed!")
}
