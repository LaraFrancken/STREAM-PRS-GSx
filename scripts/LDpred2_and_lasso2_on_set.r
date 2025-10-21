#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(remotes)
list.of.packages <- c("bigsnpr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(bigsnpr)
library(bigreadr)
library(data.table)
library(magrittr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(dplyr)
library(ggplot2)

set.seed(564)

ss <- args[1]
training_prefix <- args[2]
training_file <- args[3]
test_prefix <- args[4]
test_file <- args[5]
SNP_col <- args[6]
CHR_col <- args[7]
POS_col <- args[8]
effect_A <- args[9]
other_A <- args[10]
N_col <- args[11]
beta_se_col <- args[12]
allele_freq_col <- args[13]
out_LDpred2 <- args[14]
out_lasso2 <- args[15]
ldref_hm3_plus <- args[16]
NCORES <- args[17]
NCORES <- as.numeric(NCORES)
out_SNPs_per_set <- args[18]

#Select best model
out_comparison <- args[19]
file_path_out_comparison <- file.path(out_comparison, "Regression_results_best_per_tool") # Read the comparison file
comparison <- read.table(file_path_out_comparison, header=TRUE, sep="\t", stringsAsFactors=FALSE)
LDpred2_row <- comparison[grep("^LDpred2$", comparison$Tool), ]
model_type <- ifelse(startsWith(LDpred2_row$Parameters, "inf"), "INF model",
              ifelse(startsWith(LDpred2_row$Parameters, "grid"), "GRID model",
              ifelse(startsWith(LDpred2_row$Parameters, "auto_grid"), "AUTO GRID model",
              ifelse(startsWith(LDpred2_row$Parameters, "auto"), "AUTO model", "Unknown"))))
print(model_type)

#Grid model parameters
if (model_type == "GRID model") {
  N_h2_best <- as.numeric(strsplit(LDpred2_row$Parameters, "_")[[1]][5]) # Extract the fifth part after splitting by "_"
  N_p_best <- as.numeric(strsplit(LDpred2_row$Parameters, "_")[[1]][3])
  sparse_best <- strsplit(LDpred2_row$Parameters, "_")[[1]][7]
}

#Auto grid model parameters
shrink <- 0.95
shrink <- as.numeric(shrink)
if (model_type == "AUTO GRID model") {
  initial_p_best <- as.numeric(strsplit(LDpred2_row$Parameters, "_")[[1]][5])
  #cat("Extracted initial_p_best:", initial_p_best, "\n")
}

#Lassosum2 parameters
lassosum2_row <- comparison[grep("^lassosum2$", comparison$Tool), ]
best_delta <- as.numeric(strsplit(lassosum2_row$Parameters, "_")[[1]][2])
best_lambda <- as.numeric(strsplit(lassosum2_row$Parameters, "_")[[1]][4])
min_ratio <- args[20]
min_ratio <- as.numeric(min_ratio)

map_path <- args[21]

cat("Starting LDpred2...\n")

map <- readRDS(paste0(map_path))
  
sumstats <- bigreadr::fread2(ss)
sumstats$chr <- as.numeric(pull(sumstats, CHR_col))
sumstats$pos <- pull(sumstats, POS_col)
sumstats$a1 <- pull(sumstats, effect_A)
sumstats$a0 <- pull(sumstats, other_A)
sumstats$freq <- pull(sumstats, allele_freq_col)
sumstats$n_eff <- pull(sumstats, N_col)
sumstats$beta <- sumstats$beta_num
sumstats$beta_se <- pull(sumstats, beta_se_col)
sumstats$rsid <- pull(sumstats, SNP_col)

# Matching and QC

info_snp <- as_tibble(snp_match(sumstats, map, return_flip_and_rev = TRUE, join_by_pos=FALSE))
df_beta <- info_snp
df_beta2 <- data.frame(df_beta)

cat("SNPs matched! Building correlation matrix...\n")

# Build correlation matrix

tmp <- tempfile(tmpdir = paste0(out_LDpred2, "/tmp-data"))
for (chr in 1:22) {

    cat(chr, ".. ", sep = "")

    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in 'corr_chr'
    ind.chr3 <- match(ind.chr2, which(map$chr == chr))

    corr_chr <- readRDS(paste0(ldref_hm3_plus, "/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]

    if (chr == 1) {
      corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
    } else {
      corr$add_columns(corr_chr, nrow(corr))
    }
  }
  corr

cat("Correlation matrix built! Estimating heritability...\n")

# Estimate heritability

(ldsc <- snp_ldsc2(corr, df_beta, blocks = 200, intercept = NULL, ncores = NCORES))

ldsc_h2_est <- ldsc[["h2"]]

if (ldsc_h2_est < 0) {
  ldsc_h2_est <- 0.001
  cat("Heritability estimate is negative, changed value to 0.001\n")
}


cat("Heritability estimated! Loading in training and test data...\n")

# Training data

# List all *_snps.snplist files in the out_SNPs_per_set directory
snplist_files <- list.files(out_SNPs_per_set, pattern = "_rs.snplist$", full.names = TRUE)

# Define the output file for gene set names
gene_set_file <- file.path(out_LDpred2, "gene_sets.txt")
# Create or clear the file
file.create(gene_set_file)  

# Create a temporary directory for intermediate files
temp_dir <- file.path(out_LDpred2, "tmp-data")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

# Loop over each snplist file
for (snp_file in snplist_files) {
  set_name <- tools::file_path_sans_ext(basename(snp_file))
  
  cat(set_name, "\n", file = gene_set_file, append = TRUE)
  cat("Extracting SNPs for training data: ", set_name, "\n")

  training_bfile <- file.path(temp_dir, paste0(training_prefix, "_", set_name))
    
  training_plink_command <- paste("plink --bfile", training_file, 
  "--extract", snp_file, 
  "--make-bed", 
  "--out", training_bfile)
  system(training_plink_command)

  training_rds_file <- paste0(training_bfile, ".rds")

  if (!file.exists(training_rds_file)) {
    cat("training .rds file for set", set_name, "doesn't exist. Creating it now...\n")
    snp_readBed(paste0(training_bfile, ".bed"))
  } else {
    cat("training .rds file for set", set_name, "already exists. No new file created.\n")
  }

  #obj.bigSNP_training <- snp_attach(training_rds_file)
  assign(paste0("obj.bigSNP_training_", set_name), snp_attach(training_rds_file))

  # Access the newly created object for this SNP set
  obj.bigSNP_training <- get(paste0("obj.bigSNP_training_", set_name))

  # Impute genotypes
  assign(paste0("G_training_", set_name), obj.bigSNP_training$genotypes)
  G_training_imp <- snp_fastImputeSimple(get(paste0("G_training_", set_name)), method = "mean2", ncores = NCORES)
  assign(paste0("G_training_imp_", set_name), G_training_imp)

  # Extract sample IDs
  assign(paste0("sample_ids_training_", set_name), obj.bigSNP_training$fam$sample.ID)

  # Create map data for training
  assign(paste0("map_bigSNP_training_", set_name), dplyr::transmute(obj.bigSNP_training$map,
                          chr = as.integer(chromosome), pos = physical.pos,
                          a0 = allele1, a1 = allele2, rsid=marker.ID))

  # Map PGS and match
  map_pgs_training <- df_beta[1:4]; map_pgs_training$beta <- 1
  assign(paste0("map_pgs2_training_", set_name), snp_match(map_pgs_training, get(paste0("map_bigSNP_training_", set_name)), join_by_pos=FALSE))
}

# Test data

for (snp_file in snplist_files) {
  set_name <- tools::file_path_sans_ext(basename(snp_file))

  cat("Processing test data for SNP set: ", set_name, "\n")

  test_bfile <- file.path(temp_dir, paste0(test_prefix, "_", set_name))
    
  test_plink_command <- paste("plink --bfile", test_file, 
                                  "--extract", snp_file, 
                                  "--make-bed", 
                                  "--out", test_bfile)
  system(test_plink_command)

  test_rds_file <- paste0(test_bfile, ".rds")

  if (!file.exists(test_rds_file)) {
      cat("Test .rds file for set", set_name, "doesn't exist. Creating it now...\n")
      snp_readBed(paste0(test_bfile, ".bed"))
  } else {
      cat("Test .rds file for set", set_name, "already exists. No new file created.\n")
  }

  # Create dynamic variable names for each SNP set
  assign(paste0("obj.bigSNP_test_", set_name), snp_attach(test_rds_file))

  # Access the dynamically created object for this SNP set
  obj.bigSNP_test <- get(paste0("obj.bigSNP_test_", set_name))

  # Impute genotypes for the test data
  assign(paste0("G_test_", set_name), obj.bigSNP_test$genotypes)
  G_test_imp <- snp_fastImputeSimple(get(paste0("G_test_", set_name)), method = "mean2", ncores = NCORES)
  assign(paste0("G_test_imp_", set_name), G_test_imp)

  # Extract sample IDs for test data
  assign(paste0("sample_ids_test_", set_name), obj.bigSNP_test$fam$sample.ID)

  # Create map data for test
  assign(paste0("map_bigSNP_test_", set_name), dplyr::transmute(obj.bigSNP_test$map,
                                          chr = as.integer(chromosome), pos = physical.pos,
                                          a0 = allele1, a1 = allele2, rsid=marker.ID))

  # Map PGS and match test data with training data
  assign(paste0("map_pgs_test_", set_name), snp_match(
      get(paste0("map_pgs2_training_", set_name)), 
      get(paste0("map_bigSNP_test_", set_name)), 
      join_by_pos = FALSE
  ))
}

cat("Data loaded and SNPs matched!\n")

#### INF model
if (model_type == "INF model") {
  for (snp_file in snplist_files) {
    set_name <- tools::file_path_sans_ext(basename(snp_file))

    # Retrieve objects dynamically
    G_training <- get(paste0("G_training_", set_name))
    G_test <- get(paste0("G_test_", set_name))
    map_pgs2_training <- get(paste0("map_pgs2_training_", set_name))
    map_pgs_test <- get(paste0("map_pgs_test_", set_name))
    sample_ids_training <- get(paste0("sample_ids_training_", set_name))
    sample_ids_test <- get(paste0("sample_ids_test_", set_name))

    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = ldsc_h2_est)
    pred_inf_training <- big_prodVec(G_training, beta_inf[map_pgs2_training[["_NUM_ID_.ss"]]], ind.col = map_pgs2_training[["_NUM_ID_"]])
    pred_inf_test <- big_prodVec(G_test, beta_inf[map_pgs_test[["_NUM_ID_.ss"]]], ind.col = map_pgs_test[["_NUM_ID_"]])

    beta_inf2 <- as.data.frame(beta_inf)
    beta_inf2$rsid <- df_beta2$rsid

    pred_inf_training2 <- as.data.frame(pred_inf_training)
    colnames(pred_inf_training2) <- "pred_inf"
    pred_inf_training2$FID <- sample_ids_training
    pred_inf_training2$IID <- sample_ids_training

    pred_inf_test2 <- as.data.frame(pred_inf_test)
    colnames(pred_inf_test2) <- "pred_inf"
    pred_inf_test2$FID <- sample_ids_test
    pred_inf_test2$IID <- sample_ids_test

    set_name <- sub("_rs$", "", tools::file_path_sans_ext(basename(snp_file)))
    write.table(beta_inf2, paste0(out_LDpred2, "/beta_inf", "_", set_name), row.names = F, quote = F, sep="\t")
    write.table(pred_inf_training2, paste0(out_LDpred2, "/pred_inf_", training_prefix, "_", set_name), row.names = F, quote = F, sep="\t")
    write.table(pred_inf_test2, paste0(out_LDpred2, "/pred_inf_", test_prefix, "_", set_name), row.names = F, quote = F, sep="\t")
  }
  cat("Infinitesimal model completed!\n")
}


#### GRID model
if (model_type == "GRID model") {
  for (snp_file in snplist_files) {
    set_name <- tools::file_path_sans_ext(basename(snp_file))

    # Retrieve objects dynamically
    G_training <- get(paste0("G_training_", set_name))
    G_test <- get(paste0("G_test_", set_name))
    map_pgs2_training <- get(paste0("map_pgs2_training_", set_name))
    map_pgs_test <- get(paste0("map_pgs_test_", set_name))
    sample_ids_training <- get(paste0("sample_ids_training_", set_name))
    sample_ids_test <- get(paste0("sample_ids_test_", set_name))

    h2_seq <- round(ldsc_h2_est * N_h2_best, 4)
    p_seq <- signif(N_p_best, 2)
    params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = sparse_best)

    beta_grid <- snp_ldpred2_grid(corr, df_beta, params)
    pred_grid_training <- big_prodMat(G_training, as.matrix(beta_grid[map_pgs2_training[["_NUM_ID_.ss"]], ]), ind.col = map_pgs2_training[["_NUM_ID_"]])
    pred_grid_test <- big_prodMat(G_test, as.matrix(beta_grid[map_pgs_test[["_NUM_ID_.ss"]], ]), ind.col = map_pgs_test[["_NUM_ID_"]])

    params$colname <- paste0("p_",params$p,"_h2_",params$h2,"_sparse_",params$sparse)

    beta_grid_withcolnames <- as.data.frame(beta_grid)
    colnames_betagrid <- as.character(params[,4])
    colnames(beta_grid_withcolnames) <- colnames_betagrid
    beta_grid_withcolnames$rsid <- df_beta2$rsid

    pred_grid_withcolnames_training <- as.data.frame(pred_grid_training)
    colnames(pred_grid_withcolnames_training) <- colnames_betagrid
    pred_grid_withcolnames_training$FID <- sample_ids_training
    pred_grid_withcolnames_training$IID <- sample_ids_training

    pred_grid_withcolnames_test <- as.data.frame(pred_grid_test)
    colnames(pred_grid_withcolnames_test) <- colnames_betagrid
    pred_grid_withcolnames_test$FID <- sample_ids_test
    pred_grid_withcolnames_test$IID <- sample_ids_test

    set_name <- sub("_rs$", "", tools::file_path_sans_ext(basename(snp_file)))
    write.table(beta_grid_withcolnames, paste0(out_LDpred2, "/beta_grid", "_", set_name), row.names = F, quote = F, sep="\t")
    write.table(pred_grid_withcolnames_training, paste0(out_LDpred2, "/pred_grid_", training_prefix, "_", set_name), row.names = F, quote = F, sep="\t")
    write.table(pred_grid_withcolnames_test, paste0(out_LDpred2, "/pred_grid_", test_prefix, "_", set_name), row.names = F, quote = F, sep="\t")
  }
  cat("Grid model completed!\n")
}


#### AUTO GRID model
if (model_type == "AUTO GRID model") {
  for (snp_file in snplist_files) {
    set_name <- tools::file_path_sans_ext(basename(snp_file))

    # Retrieve objects dynamically
    G_training <- get(paste0("G_training_", set_name))
    G_test <- get(paste0("G_test_", set_name))
    map_pgs2_training <- get(paste0("map_pgs2_training_", set_name))
    map_pgs_test <- get(paste0("map_pgs_test_", set_name))
    sample_ids_training <- get(paste0("sample_ids_training_", set_name))
    sample_ids_test <- get(paste0("sample_ids_test_", set_name))

    coef_shrink <- shrink
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc_h2_est,
                        vec_p_init = initial_p_best,
                        allow_jump_sign = FALSE, shrink_corr = coef_shrink,
                        ncores = NCORES)
    auto <- multi_auto[[1]]

    range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
    keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

    beta_auto_grid <- sapply(multi_auto[keep], function(auto) auto$beta_est)
    pred_auto_grid_training <- big_prodMat(G_training, as.matrix(beta_auto_grid[map_pgs2_training[["_NUM_ID_.ss"]], ]), ind.col = map_pgs2_training[["_NUM_ID_"]])
    pred_auto_grid_test <- big_prodMat(G_test, as.matrix(beta_auto_grid[map_pgs_test[["_NUM_ID_.ss"]], ]), ind.col = map_pgs_test[["_NUM_ID_"]])

    colnames_auto <- initial_p_best[keep]
    colnames_auto <- as.data.frame(colnames_auto)
    colnames_auto$colname <- paste0("p_init_",colnames_auto$colnames_auto)

    beta_auto_grid_withcolnames <- as.data.frame(beta_auto_grid)
    colnames_beta_auto_grid <- as.character(colnames_auto[,2])
    colnames(beta_auto_grid_withcolnames) <- colnames_beta_auto_grid
    beta_auto_grid_withcolnames$rsid <- df_beta2$rsid

    pred_auto_grid_withcolnames_training <- as.data.frame(pred_auto_grid_training)
    colnames(pred_auto_grid_withcolnames_training) <- colnames_beta_auto_grid
    pred_auto_grid_withcolnames_training$FID <- sample_ids_training
    pred_auto_grid_withcolnames_training$IID <- sample_ids_training

    pred_auto_grid_withcolnames_test <- as.data.frame(pred_auto_grid_test)
    colnames(pred_auto_grid_withcolnames_test) <- colnames_beta_auto_grid
    pred_auto_grid_withcolnames_test$FID <- sample_ids_test
    pred_auto_grid_withcolnames_test$IID <- sample_ids_test

    set_name <- sub("_rs$", "", tools::file_path_sans_ext(basename(snp_file)))
    write.table(beta_auto_grid_withcolnames, paste0(out_LDpred2, "/beta_auto_grid", "_", set_name), row.names = F, quote = F, sep="\t")
    write.table(pred_auto_grid_withcolnames_training, paste0(out_LDpred2,"/pred_auto_grid_", training_prefix, "_", set_name), row.names = F, quote = F, sep="\t")
    write.table(pred_auto_grid_withcolnames_test, paste0(out_LDpred2, "/pred_auto_grid_", test_prefix, "_", set_name), row.names = F, quote = F, sep="\t")
  }
  cat("Auto grid model completed!\n")
}


#### AUTO model
if (model_type == "AUTO model") {
  for (snp_file in snplist_files) {
    set_name <- tools::file_path_sans_ext(basename(snp_file))

    # Retrieve objects dynamically
    G_training <- get(paste0("G_training_", set_name))
    G_test <- get(paste0("G_test_", set_name))
    map_pgs2_training <- get(paste0("map_pgs2_training_", set_name))
    map_pgs_test <- get(paste0("map_pgs_test_", set_name))
    sample_ids_training <- get(paste0("sample_ids_training_", set_name))
    sample_ids_test <- get(paste0("sample_ids_test_", set_name))

    beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    pred_auto_training <- big_prodVec(G_training, beta_auto[map_pgs2_training[["_NUM_ID_.ss"]]], ind.col = map_pgs2_training[["_NUM_ID_"]])
    pred_auto_test <- big_prodVec(G_test, beta_auto[map_pgs_test[["_NUM_ID_.ss"]]], ind.col = map_pgs_test[["_NUM_ID_"]])

    beta_auto2 <- as.data.frame(beta_auto)
    beta_auto2$rsid <- df_beta2$rsid

    pred_auto2_training <- as.data.frame(pred_auto_training)
    colnames(pred_auto2_training) <- "pred_auto"
    pred_auto2_training$FID <- sample_ids_training
    pred_auto2_training$IID <- sample_ids_training

    pred_auto2_test <- as.data.frame(pred_auto_test)
    colnames(pred_auto2_test) <- "pred_auto"
    pred_auto2_test$FID <- sample_ids_test
    pred_auto2_test$IID <- sample_ids_test

    set_name <- sub("_rs$", "", tools::file_path_sans_ext(basename(snp_file)))

    write.table(beta_auto2, paste0(out_LDpred2, "/beta_auto", "_", set_name), row.names = F, quote = F, sep="\t")
    write.table(pred_auto2_training, paste0(out_LDpred2, "/pred_auto_", training_prefix, "_", set_name), row.names = F, quote = F, sep="\t")
    write.table(pred_auto2_test, paste0(out_LDpred2, "/pred_auto_", test_prefix, "_", set_name), row.names = F, quote = F, sep="\t")
  }
  cat("Auto model completed!\n")
}


#### LASSOSUM 2

cat("Running lassosum2...\n")
for (snp_file in snplist_files) {
  set_name <- tools::file_path_sans_ext(basename(snp_file))
  cat("Running snp_lassosum2 for SNP set:", set_name, "\n")

  # Save set_name into gene_sets.txt, one per line
  write(set_name, 
        file = file.path(out_lasso2, "gene_sets.txt"), 
        append = TRUE)

  # Retrieve objects dynamically
  G_training <- get(paste0("G_training_", set_name))
  G_test <- get(paste0("G_test_", set_name))
  map_pgs2_training <- get(paste0("map_pgs2_training_", set_name))
  map_pgs_test <- get(paste0("map_pgs_test_", set_name))
  sample_ids_training <- get(paste0("sample_ids_training_", set_name))
  sample_ids_test <- get(paste0("sample_ids_test_", set_name))

  # Compute correlation matrix
  #corr <- snp_cor(G_training, ind.col = map_pgs2_training[["_NUM_ID_"]], ncores = NCORES)

  beta_lassosum2 <- snp_lassosum2(corr, df_beta, delta = best_delta, nlambda = best_lambda, lambda.min.ratio = min_ratio)
  pred_lassosum2_training <- big_prodMat(G_training, as.matrix(beta_lassosum2[map_pgs2_training[["_NUM_ID_.ss"]], ]), ind.col = map_pgs2_training[["_NUM_ID_"]])
  pred_lassosum2_test <- big_prodMat(G_test, as.matrix(beta_lassosum2[map_pgs_test[["_NUM_ID_.ss"]], ]), ind.col = map_pgs_test[["_NUM_ID_"]])

  params_lasso <- attr(beta_lassosum2, "grid_param")
  params_lasso$colname <- paste0("delta_",params_lasso$delta,"_lambda_",params_lasso$lambda)

  beta_lasso_withcolnames <- as.data.frame(beta_lassosum2)
  colnames_betalasso <- as.character(params_lasso[,6])
  colnames(beta_lasso_withcolnames) <- colnames_betalasso
  beta_lasso_withcolnames$rsid <- df_beta2$rsid

  pred_lasso_withcolnames_training <- as.data.frame(pred_lassosum2_training)
  colnames(pred_lasso_withcolnames_training) <- colnames_betalasso
  pred_lasso_withcolnames_training$FID <- sample_ids_training
  pred_lasso_withcolnames_training$IID <- sample_ids_training

  pred_lasso_withcolnames_test <- as.data.frame(pred_lassosum2_test)
  colnames(pred_lasso_withcolnames_test) <- colnames_betalasso
  pred_lasso_withcolnames_test$FID <- sample_ids_test
  pred_lasso_withcolnames_test$IID <- sample_ids_test

  set_name <- sub("_rs$", "", tools::file_path_sans_ext(basename(snp_file)))
  write.table(beta_lasso_withcolnames, file.path(out_lasso2, paste0("beta_lasso2_", set_name)), row.names = F, quote = F, sep = "\t")
  write.table(pred_lasso_withcolnames_training, file.path(out_lasso2, paste0("pred_lasso2_", training_prefix, "_", set_name)), row.names = F, quote = F, sep = "\t")
  write.table(pred_lasso_withcolnames_test, file.path(out_lasso2, paste0("pred_lasso2_", test_prefix, "_", set_name)), row.names = F, quote = F, sep = "\t")
}

cat("Lassosum2 completed!\n")
