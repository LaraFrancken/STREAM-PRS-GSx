#!/usr/bin/bash

#Script that combines all the other scripts: calculates PRS for PRSice, lassosum, PRS-cs, LDpred-2 and lassosum2.

#Parameters - to adapt

#Tools
PRSice_path=""
PRScs_path=""

#Model 
binary_trait="TRUE" #(TRUE=binary, FALSE=continuous)
cross_validation="FALSE" #(TRUE=use cross validation, FALSE=don't use cross validation)
covariates="" #Comma separated list of covariates to include in regression model e.g. "Age+Sex,Sex" will do regression analysis for PHENO ~ PRS + Age+Sex and for PHENO ~ PRS + Sex; if e.g. "Age,Sex" it will do regression analysis for PHENO ~ PRS + Age and for PHENO ~ PRS + Sex It will NOT do PHENO ~ PRS + Age + Sex; leave empty if you have no covariates.
gene_set="FALSE" #(TRUE=yes, FALSE=no)

#Output path
out=""

#GWAS summary statistics
#If chr_pos_ID and rsID are in the same file, fill in the same GWAS twice
GWAS_file=""
GWAS_file_rsID=""

#GWAS column names
snp_col=""
rsID_col=""
chromosome_col=""
pos_col=""
effect_allele=""
noneffect_allele=""
beta_col=""
beta_se_col=""
p_value_col=""
N_col=""
allele_freq_col=""

#Size of GWAS (cases+controls)
GWAS_size=""

#Training data
training_file=""
training_file_prefix=""
training_file_rsID=""
training_file_rsID_prefix=""

#Test data
test_file=""
test_file_prefix=""
test_file_rsID=""
test_file_rsID_prefix=""

#Phenotype file (FID, IID, PHENO)
pheno=""

#PC files (FID, IID, PC1, PC2, ..., PCN)
PC_training=""
PC_test=""

#Covariate file (FID, IID, cov_1, cov_2, ..., cov_N): this file is optional
cov_file=""

#Gene set file (Name set, URL, Genes): this file is optional
#Make sure there is no missing newline at end of the file!
MSigDB=""
#GTF file with the genome boundary of the genetic elements within the human genome: this file is optional
GTF=""
#Gene boundary padding: padding to add upstream (5') and downstream (3') to the gene region. 
#Values are in kilobases (kb).
pad5="5"
pad3="10"
#"Base" is a file containing the SNPs of all gene sets together
base="TRUE" #(TRUE=make "Base" gene set, FALSE=no need for "Base")

#Cores
cores="1"

#PRSice specific parameters
pval_thresholds="5e-08,1e-05,0.0001,0.001,0.05,0.1,0.5,1"
clumpkb="250"
clumpp="1.000000"
clumpr2="0.100000"

#Lassosum specific parameters
ld_blocks_lasso="EUR.hg38"
lasso_thresholding_values="0.1,0.2,0.4,0.5,0.7,0.9,1"
lambda_values="exp(seq(log(0.001), log(0.1), length.out = 20))"

#PRS-CS specific parameters
reference_files_PRScs=""
phi_values="1e+00,1e-02,1e-04,1e-06"

#LDpred2 and lassosum2 specific parameters
ldref_hm3_plus=""
values_h2_grid="c(0.3, 0.7, 1, 1.4)"
values_p_grid="seq_log(1e-5, 1, length.out = 21)"
initial_p_auto_grid="seq_log(1e-4, 0.2, length.out = 30)"
delta="c(0.001, 0.01, 0.1, 1)"
nlambda="30"
lambda_min_ratio="0.01"

#Parameters - not to adapt

#Output paths
out_PRSice=${out}/PRSice
out_lasso=${out}/lasso
out_PRScs=${out}/PRScs
out_LDpred2=${out}/LDpred2
out_lasso2=${out}/lasso2
out_comparison=${out}/comparison

mkdir -p ${out_PRSice}
mkdir -p ${out_lasso}
mkdir -p ${out_PRScs}
mkdir -p ${out_LDpred2}
mkdir -p ${out_lasso2}
mkdir -p ${out_comparison}

if [ "$gene_set" == "TRUE" ]; then
    out_SNPs_per_set=${out}/SNPS_per_set
    out_PRSet=${out}/PRSet
    out_PRSice_on_Set=${out}/PRSice_on_set
    out_lassosum_on_set=${out}/lassosum_on_set
    out_PRScs_on_set=${out}/PRScs_on_set
    out_LDpred2_on_set=${out}/LDpred2_on_set
    out_lasso2_on_set=${out}/lasso2_on_set
fi

if [ "$gene_set" == "TRUE" ]; then
    mkdir -p ${out_SNPs_per_set}
    mkdir -p ${out_PRSet}
    mkdir -p ${out_PRSice_on_Set}
    mkdir -p ${out_lassosum_on_set}
    mkdir -p ${out_PRScs_on_set}
    mkdir -p ${out_LDpred2_on_set}
    mkdir -p ${out_lasso2_on_set}
fi

#GWAS summary statistics
edited_GWAS_file=${GWAS_file}_edited_by_pipeline
edited_GWAS_for_PRS_CS=${GWAS_file}_edited_by_pipeline_for_PRS_CS
edited_GWAS_for_LDpred2=${GWAS_file}_edited_by_pipeline_for_LDpred2

#Run scripts

#Edit GWAS
Rscript edit_GWAS.r "$GWAS_file" "$GWAS_file_rsID" "$p_value_col" "$beta_col" "$snp_col" "$rsID_col" "$effect_allele" "$noneffect_allele" "$pos_col"

#Run PRSice
./prsice.bash "$PRSice_path" "$edited_GWAS_file" "$training_file" "$training_file_prefix" "$snp_col" "$chromosome_col" "$pos_col" "$effect_allele" "$noneffect_allele" "$clumpkb" "$clumpp" "$clumpr2" "$pval_thresholds" "$out_PRSice" "$test_file" "$test_file_prefix"

#Run lassosum
Rscript lassosum.r "$test_file_prefix" "$test_file" "$edited_GWAS_file" "$training_file" "$ld_blocks_lasso" "$chromosome_col" "$pos_col" "$effect_allele" "$noneffect_allele" "$N_col" "$lasso_thresholding_values" "$out_lasso" "$training_file_prefix" "$lambda_values"

#Run PRS-CS
./PRS_CS.bash "$PRScs_path" "$edited_GWAS_for_PRS_CS" "$test_file_rsID" "$test_file_prefix" "$training_file_rsID" "$training_file_prefix" "$out_PRScs" "$reference_files_PRScs" "$GWAS_size" "$phi_values"

#Run LDPred-2 and lassosum2
Rscript LDpred2_and_lassosum2.r "$edited_GWAS_for_LDpred2" "$training_file_rsID_prefix" "$training_file_rsID" "$test_file_rsID_prefix" "$test_file_rsID" "$rsID_col" "$chromosome_col" "$pos_col" "$effect_allele" "$noneffect_allele" "$N_col" "$beta_se_col" "$allele_freq_col" "$out_LDpred2" "$out_lasso2" "$ldref_hm3_plus" "$cores" "$values_h2_grid" "$values_p_grid" "$initial_p_auto_grid" "$delta" "$nlambda" "$lambda_min_ratio"

#Get best PRS
if [ "$binary_trait" == "TRUE" ]; then
	Rscript get_best_PRS.r "$out_PRSice" "$out_lasso" "$out_PRScs" "$out_LDpred2" "$out_lasso2" "$out_comparison" "$test_file_prefix" "$test_file_rsID_prefix" "$training_file_prefix" "$training_file_rsID_prefix" "$PC_test" "$PC_training" "$pheno" "$lasso_thresholding_values" "$phi_values" "$cross_validation" "$covariates" "$cov_file"
elif [ "$binary_trait" == "FALSE" ]; then
		Rscript get_best_PRS_linear.r "$out_PRSice" "$out_lasso" "$out_PRScs" "$out_LDpred2" "$out_lasso2" "$out_comparison" "$test_file_prefix" "$test_file_rsID_prefix" "$training_file_prefix" "$training_file_rsID_prefix" "$PC_test" "$PC_training" "$pheno" "$lasso_thresholding_values" "$phi_values" "$cross_validation" "$covariates" "$cov_file"
else 
	echo "Indicate if trait is binary (TRUE) or not (FALSE)"
fi

#Run gene-set extension scripts

#SNP extraction
if [ "$gene_set" == "TRUE" ]; then
    ./get_SNPs_per_set.bash "$MSigDB" "$GTF" "$out_SNPs_per_set" "$training_file" "$training_file_rsID" "$pad5" "$pad3" "$base"
fi

#Run PRSet
if [ "$gene_set" == "TRUE" ]; then
    ./PRSet.bash "$PRSice_path" "$edited_GWAS_file" "$MSigDB" "$GTF" "$training_file" "$training_file_prefix" "$snp_col" "$chromosome_col" "$pos_col" "$effect_allele" "$noneffect_allele" "$clumpkb" "$clumpp" "$clumpr2" "$pval_thresholds" "$out_PRSet" "$test_file" "$test_file_prefix" "$out_comparison"
fi

#Run PRSice on every gene set
if [ "$gene_set" == "TRUE" ]; then
   ./PRSice_on_set.bash "$PRSice_path" "$edited_GWAS_file" "$MSigDB" "$GTF" "$training_file" "$training_file_prefix" "$snp_col" "$chromosome_col" "$pos_col" "$effect_allele" "$noneffect_allele" "$clumpkb" "$clumpp" "$clumpr2" "$pval_thresholds" "$out_PRSice_on_Set" "$test_file" "$test_file_prefix" "$out_comparison" "$out_SNPs_per_set"
fi

#Run lassosum on every gene set
if [ "$gene_set" == "TRUE" ]; then
   ./lassosum_on_set.bash "$test_file_prefix" "$test_file" "$edited_GWAS_file" "$training_file" "$ld_blocks_lasso" "$chromosome_col" "$pos_col" "$effect_allele" "$noneffect_allele" "$N_col" "$out_lasso" "$training_file_prefix" "$lambda_values" "$out_comparison" "$out_lassosum_on_set" "$out_SNPs_per_set"
fi

#Run PRS-CS on every gene set
if [ "$gene_set" == "TRUE" ]; then
    ./PRScs_on_set.bash "$PRScs_path" "$edited_GWAS_for_PRS_CS" "$test_file_rsID" "$test_file_prefix" "$training_file_rsID" "$training_file_prefix" "$out_PRScs_on_set" "$reference_files_PRScs" "$GWAS_size" "$out_comparison" "$out_SNPs_per_set"
fi

#Run LDPred-2 and lassosum2 on every gene set
if [ "$gene_set" == "TRUE" ]; then
    Rscript LDpred2_and_lasso2_on_set.r "$edited_GWAS_for_LDpred2" "$training_file_prefix" "$training_file_rsID" "$test_file_prefix" "$test_file_rsID" "$rsID_col" "$chromosome_col" "$pos_col" "$effect_allele" "$noneffect_allele" "$N_col" "$beta_se_col" "$allele_freq_col" "$out_LDpred2_on_set" "$out_lasso2_on_set" "$ldref_hm3_plus" "$cores" "$out_SNPs_per_set" "$out_comparison" "$lambda_min_ratio"   
fi

#Get best gene set
if [ "$gene_set" == "TRUE" ]; then
    if [ "$binary_trait" == "TRUE" ]; then
        Rscript get_best_set.r "$out_PRSet" "$out_PRSice_on_Set" "$out_lassosum_on_set" "$out_PRScs_on_set" "$out_LDpred2_on_set" "$out_lasso2_on_set" "$out_comparison" "$test_file_prefix" "$test_file_rsID_prefix" "$training_file_prefix" "$training_file_rsID_prefix" "$PC_test" "$PC_training" "$pheno" "$covariates" "$cov_file" "$cross_validation"
    elif [ "$binary_trait" == "FALSE" ]; then
	    Rscript get_best_set_linear.r "$out_PRSet" "$out_PRSice_on_Set" "$out_lassosum_on_set" "$out_PRScs_on_set" "$out_LDpred2_on_set" "$out_lasso2_on_set" "$out_comparison" "$test_file_prefix" "$test_file_rsID_prefix" "$training_file_prefix" "$training_file_rsID_prefix" "$PC_test" "$PC_training" "$pheno" "$covariates" "$cov_file" "$cross_validation"
    else 
	    echo "Indicate if trait is binary (TRUE) or not (FALSE)"
    fi  
fi
