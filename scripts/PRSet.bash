#!/usr/bin/bash

PRSice_path=$1
GWAS_file=$2
MSigDB_file=$3
GTF_file=$4
training_data=$5
training_prefix=$6
snp_col=$7
chromosome_col=$8
pos_col=$9
effect_allele=${10}
noneffect_allele=${11}
clumpkb=${12}
clumpp=${13}
clumpr2=${14}
pval_thresholds=${15}
output_PRSet=${16}
test_data=${17}
test_prefix=${18}
out_comparison=${19}

# Get the best threshold from PRSice
best_threshold=$(grep "PRSice" ${out_comparison}/Regression_results_best_per_tool | awk '{split($2, a, "_"); print a[length(a)]}')

echo "Best threshold PRSice is ${best_threshold}"
echo "Start PRSet..."


#PRSet for training data
Rscript ${PRSice_path}/PRSice.R --prsice ${PRSice_path}/PRSice_linux --base ${GWAS_file} --target ${training_data} --msigdb ${MSigDB_file} --gtf ${GTF_file} --snp ${snp_col} --chr ${chromosome_col} --bp ${pos_col} --A1 ${effect_allele} --A2 ${noneffect_allele} --stat beta_num --beta --pvalue P_num --clump-kb ${clumpkb} --clump-p ${clumpp} --clump-r2 ${clumpr2} --bar-levels ${best_threshold} --fastscore --no-regress --no-full --print-snp --out ${output_PRSet}/${training_prefix}

echo "PRSice for training data completed...."

#PRSet for test data
Rscript ${PRSice_path}/PRSice.R --prsice ${PRSice_path}/PRSice_linux --base ${GWAS_file} --target ${test_data} --msigdb ${MSigDB_file} --gtf ${GTF_file} --snp ${snp_col} --chr ${chromosome_col} --bp ${pos_col} --A1 ${effect_allele} --A2 ${noneffect_allele} --stat beta_num --beta --pvalue P_num --extract ${output_PRSet}/${training_prefix}.snp --no-clump --bar-levels ${best_threshold} --fastscore --no-regress --no-full --print-snp --out ${output_PRSet}/${test_prefix}


echo "PRSet for test data completed...."
echo "PRSet completed!"

