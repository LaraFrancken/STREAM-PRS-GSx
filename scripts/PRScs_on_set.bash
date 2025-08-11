#!/usr/bin/bash

PRScs_path=$1
GWAS_file=$2
test_data=$3
test_prefix=$4
training_data=$5
training_prefix=$6
out_PRScs_on_set=$7
PRS_CS_ref_files=$8
GWAS_size=$9
out_comparison=${10}
out_SNPs_per_set=${11}

# Get number of cores for parallel execution
export MKL_NUM_THREADS=$SLURM_CPUS_ON_NODE
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_ON_NODE
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

### Step 1: Get the best phi from Regression_results_best_per_tool
best_phi=$(grep "PRS-CS" "${out_comparison}/Regression_results_best_per_tool" | awk -F'phi_' '{print $2}' | awk '{print $1}' | sed 's/\./-/g')
echo "Best phi PRS-CS is ${best_phi}"

# Create temp directory if it doesn't exist
mkdir -p "${out_PRScs_on_set}/temp"

# Define the output file for gene set names
gene_set_file="${out_PRScs_on_set}/gene_sets.txt"
> "$gene_set_file"  # Create or clear the file

### Step 2: Extracting SNPs for test and training data
temp_dir=$(mktemp -d)
# Loop through the gene sets in out_SNPs_per_set directory
for snp_file in ${out_SNPs_per_set}/*_rs.snplist; do
    set_name=$(basename "$snp_file" _rs.snplist)

    # Save the gene set name to the file
    echo "$set_name" >> "$gene_set_file"

    echo "Extracting SNPs for test data: $set_name"
    plink --bfile ${test_data} --extract "$snp_file" --make-bed --out ${temp_dir}/${test_prefix}_${set_name}
    
    echo "Extracting SNPs for training data: $set_name"
    plink --bfile ${training_data} --extract "$snp_file" --make-bed --out ${temp_dir}/${training_prefix}_${set_name}

    echo "Extracting SNPs for sst_file: $set_name"
    head -n 1 ${GWAS_file} > ${temp_dir}/sst_${set_name} # Print the header first
    # Then extract the SNPs that match, excluding any that are not in the snp_file
    awk 'NR==FNR {a[$1]; next} $1 in a {print}' ${snp_file} ${GWAS_file} >> ${temp_dir}/sst_${set_name}
    #echo "First 5 lines of ${temp_dir}/sst_${set_name}:"
    #head -n 5 ${temp_dir}/sst_${set_name}

    ### Step 3: Run PRS-CS for each gene set
    echo "Running PRS-CS for gene set: $set_name"

    # Step 3.1: Run PRS-CS for each gene set using the training files
    seq 1 22 | parallel -j 22 "python ${PRScs_path}/PRScs.py \
        --ref_dir=${PRS_CS_ref_files} \
        --bim_prefix=${temp_dir}/${training_prefix}_${set_name} \
        --sst_file=${temp_dir}/sst_${set_name} \
        --n_gwas=${GWAS_size} \
        --chrom={} \
        --phi=${best_phi} \
        --out_dir=${out_PRScs_on_set}/temp/${set_name}"

    # Function to check missing files
    check_missing_files() {
        missing_files=()
        for i in {1..22}; do
            file="${out_PRScs_on_set}/temp/${set_name}_pst_eff_a1_b0.5_phi${best_phi}_chr${i}.txt"
            if [ ! -f "$file" ]; then
                missing_files+=("$i")
            fi
        done
    }

    # Step 3.2: Retry PRS-CS for missing chromosomes (if any)
    j=1
    while true; do
        check_missing_files
        
        if [ ${#missing_files[@]} -eq 0 ]; then
            echo "PRS-CS completed for gene set ${set_name} with phi ${best_phi}"
            break
        
        elif [ $j -eq 10 ]; then
            echo "Tried PRS-CS ten times for ${set_name} and still no estimates for all chromosomes... exiting..."
            break
            
        else
            # Output the missing files
            echo "The following chromosomes are missing for ${set_name} and phi ${best_phi}: ${missing_files[@]}"

            # Run PRS-CS for each missing chromosome
            for number in "${missing_files[@]}"; do
                echo "Running PRScs.py for gene set ${set_name} and chromosome ${number}..."
                python "${PRScs_path}/PRScs.py" \
                    --ref_dir="${PRS_CS_ref_files}" \
                    --bim_prefix="${temp_dir}/${training_prefix}_${set_name}" \
                    --sst_file="${temp_dir}/sst_${set_name}" \
                    --n_gwas="${GWAS_size}" \
                    --chrom="${number}" \
                    --phi="${best_phi}" \
                    --out_dir="${out_PRScs_on_set}/temp/${set_name}"
            done
            ((j++))
        fi
    done

    # Step 3.3: Concatenate the results from all chromosomes
    for i in {1..22}; do
        cat ${out_PRScs_on_set}/temp/${set_name}_pst_eff_a1_b0.5_phi${best_phi}_chr${i}.txt >> ${out_PRScs_on_set}/${set_name}_pst_eff_a1_b0.5_phi${best_phi}_all.txt
    done

    # Step 3.4: Run PLINK on the concatenated results
    plink --bfile ${temp_dir}/${training_prefix}_${set_name} --score ${out_PRScs_on_set}/${set_name}_pst_eff_a1_b0.5_phi${best_phi}_all.txt 2 4 6 --out ${out_PRScs_on_set}/${training_prefix}_${set_name}_phi_${best_phi}
    plink --bfile ${temp_dir}/${test_prefix}_${set_name} --score ${out_PRScs_on_set}/${set_name}_pst_eff_a1_b0.5_phi${best_phi}_all.txt 2 4 6 --out ${out_PRScs_on_set}/${test_prefix}_${set_name}_phi_${best_phi}
    
done

# Clean up temporary files
rm -r ${out_PRScs_on_set}/temp

echo "Finished processing all gene sets for PRS-CS."