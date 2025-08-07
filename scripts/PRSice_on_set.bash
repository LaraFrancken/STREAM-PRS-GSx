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
output_PRSice_on_Set=${16}
test_data=${17}
test_prefix=${18}
out_comparison=${19}
out_SNPs_per_set=${20}

### Step 1: Get the best threshold from PRSice
best_threshold=$(grep "PRSice" ${out_comparison}/Regression_results_best_per_tool | awk '{split($2, a, "_"); print a[length(a)]}')

echo "Best threshold PRSice is ${best_threshold}"

### Step 2: Loop through the generated SNP files and run PRSice
for snp_file in ${out_SNPs_per_set}/*_snps.snplist; do
    set_name=$(basename "$snp_file" _snps.snplist)

    echo "Running PRSice for training data on set: $set_name"

    Rscript "${PRSice_path}/PRSice.R" \
        --prsice "${PRSice_path}/PRSice_linux" \
        --base "${GWAS_file}" \
        --target "${training_data}" \
        --extract "$snp_file" \
        --snp "${snp_col}" \
        --chr "${chromosome_col}" \
        --bp "${pos_col}" \
        --A1 "${effect_allele}" \
        --A2 "${noneffect_allele}" \
        --stat "beta_num" \
        --beta \
        --pvalue "P_num" \
        --clump-kb "${clumpkb}" \
        --clump-p "${clumpp}" \
        --clump-r2 "${clumpr2}" \
        --bar-levels "${best_threshold}" \
        --fastscore \
        --no-regress \
        --no-full \
        --print-snp \
        --out "${output_PRSice_on_Set}/${set_name}_training"

    echo "PRSice for training data on set $set_name completed!"

    echo "Running PRSice for test data on set: $set_name"

    Rscript "${PRSice_path}/PRSice.R" \
        --prsice "${PRSice_path}/PRSice_linux" \
        --base "${GWAS_file}" \
        --target "${test_data}" \
        --extract "${output_PRSice_on_Set}/${set_name}_training.snp" \
        --snp "${snp_col}" \
        --chr "${chromosome_col}" \
        --bp "${pos_col}" \
        --A1 "${effect_allele}" \
        --A2 "${noneffect_allele}" \
        --stat "beta_num" \
        --beta \
        --pvalue "P_num" \
        --no-clump \
        --bar-levels "${best_threshold}" \
        --fastscore \
        --no-regress \
        --no-full \
        --print-snp \
        --out "${output_PRSice_on_Set}/${set_name}_test"

    echo "PRSice for test data on set $set_name completed!"

done

echo "All PRSice analyses completed!"

### Step 3: Merge files
merged_training="${output_PRSice_on_Set}/${training_prefix}.all_score"
merged_test="${output_PRSice_on_Set}/${test_prefix}.all_score"

# Initialize the training file with only FID and IID columns
first_training_file=$(ls ${output_PRSice_on_Set}/*_training.all_score | head -n 1)
awk '{print $1, $2}' "$first_training_file" > "$merged_training"
# Initialize the test file with only FID and IID columns
first_test_file=$(ls ${output_PRSice_on_Set}/*_test.all_score | head -n 1)
awk '{print $1, $2}' "$first_test_file" > "$merged_test"

#echo "First 5 lines of the initialized training file:"
#head -n 5 "$merged_training"
#echo "First 5 lines of the initialized test file:"
#head -n 5 "$merged_test"

# Function to process and merge files
merge_files() {
    local filetype=$1  # "training" or "test"
    local merged_file=$2

    for file in ${output_PRSice_on_Set}/*_${filetype}.all_score; do
        # Extract base set name, removing _training or _test suffix
        set_name=$(basename "$file" | sed "s/_${filetype}.all_score//")
        echo "Merging ${filetype} set: $set_name"

        # Add prefix to columns and merge on FID, IID
        awk -v prefix="${set_name}_" 'NR==1 {
            # Prefix all columns starting from the third
            for (i=3; i<=NF; i++) $i=prefix $i
        } 1' "$file" | paste "$merged_file" - > "${merged_file}.tmp"

        mv "${merged_file}.tmp" "$merged_file"
    done
}

# Merge training and test files
merge_files "training" "$merged_training"
merge_files "test" "$merged_test"

echo "Merging completed."
echo "PRSice on set completed!"
