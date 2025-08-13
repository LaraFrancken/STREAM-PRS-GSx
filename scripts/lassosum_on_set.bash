#!/usr/bin/bash

test_prefix=$1
test_data_chr_pos_ID=$2
edited_GWAS_file=$3
ld_file=$4
ld_blocks_lasso=$5
chromosome_col=$6
pos_col=$7
effect_allele=$8
noneffect_allele=$9
N_col=${10}
out_lasso=${11}
training_prefix=${12}
lambda_values=${13}
out_comparison=${14}
out_lassosum_on_set=${15}
out_SNPs_per_set=${16}

# Get best s and L
best_s_lassosum=$(grep -w "lassosum" ${out_comparison}/Regression_results_best_per_tool | awk '{split($2, a, "_"); print a[2]}')
best_L_lassosum=$(grep -w "lassosum" ${out_comparison}/Regression_results_best_per_tool | awk '{split($2, a, "_"); print a[4]}')
echo "Best thresholds lassosum: s = ${best_s_lassosum} and L = ${best_L_lassosum}"

# Create a temporary directory
temp_dir=$(mktemp -d)

# Define the output file for gene set names
gene_set_file="${out_lassosum_on_set}/gene_sets.txt"
> "$gene_set_file"  # Create or clear the file

for snp_file in ${out_SNPs_per_set}/*_snps.snplist; do
    set_name=$(basename "$snp_file" _snps.snplist)

    # Save the gene set name to the file
    echo "$set_name" >> "$gene_set_file"
    
    echo "Extracting SNPs for test data: $set_name"
    plink --bfile ${test_data_chr_pos_ID} --extract "$snp_file" --make-bed --out ${temp_dir}/${test_prefix}_${set_name}
    
    echo "Extracting SNPs for training data: $set_name"
    plink --bfile ${ld_file} --extract "$snp_file" --make-bed --out ${temp_dir}/${training_prefix}_${set_name}
    
    echo "Running lassosum for set: $set_name"
    Rscript lassosum.r "${test_prefix}_${set_name}" "${temp_dir}/${test_prefix}_${set_name}" "$edited_GWAS_file" "${temp_dir}/${training_prefix}_${set_name}" "$ld_blocks_lasso" "$chromosome_col" "$pos_col" "$effect_allele" "$noneffect_allele" "$N_col" "$best_s_lassosum" "$out_lassosum_on_set" "${training_prefix}_${set_name}" "$best_L_lassosum"

    # Remove temporary files after lassosum runs
    echo "Cleaning up temporary files for $set_name"
    rm -rf ${temp_dir}/${test_prefix}_${set_name}*
    rm -rf ${temp_dir}/${training_prefix}_${set_name}*

done

# Remove the temporary directory
rmdir "$temp_dir"
echo "All lassosum analyses completed!"



