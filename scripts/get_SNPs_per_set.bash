#!/usr/bin/bash

MSigDB_file=$1
GTF_file=$2
out_SNPs_per_set=$3
training_file=$4
training_file_rsID=$5
pad5=$6
pad3=$7
base=$8

# Convert padding from kb to base pairs
pad5_bp=$((pad5 * 1000))
pad3_bp=$((pad3 * 1000))

### Extract chromosome position per gene set

echo "Creating region files and extract SNPs ..."

while IFS=$'\t' read -r set_name url genes; do # Make sure the file ends in a new line!
    # 'IFS=$'\t'' ensures splitting by tab, which is typical in MSigDB files
    # 'set_name' will hold the first column (gene set name)
    # 'url' will hold the second column
    # 'genes' will hold the rest of the columns (the gene names)

    #echo "Processing gene set: $set_name"

    # Create a temporary file to store the regions for this gene set
    temp_regions_file="${out_SNPs_per_set}/${set_name}_regions.txt"
    > "$temp_regions_file"  # Clear any existing content

    # Split the genes into an array to process each gene individually
    IFS=' ' read -ra gene_array <<< "$genes"
    gene_array=($genes)
    #echo "Genes: ${gene_array[@]}"

    # Convert the list of genes into a regex pattern
    gene_pattern=$(echo "${gene_array[@]}" | tr ' ' '|')
    # Use `awk` to extract all regions for the genes in the list at once
    awk -v pattern="$gene_pattern" -v pad5="$pad5" -v pad3="$pad3" -F'\t' '
    $3 == "gene" && $9 ~ pattern {
        match($9, /gene_name "([^"]+)"/, arr)   # Extract gene name from the 9th column
        if (arr[1] != "") {
            # Get the original start and end coordinates
            start = $4
            end = $5
            strand = $7

            #print "Start before padding:", start > "/dev/stderr"
            #print "Strand: ", strand > "/dev/stderr"

            # Apply padding based on the strand
            if (strand == "+") {
                start = start - pad5
                end = end + pad3
            } else if (strand == "-") {
                start = start - pad3
                end = end + pad5
            }

            # Ensure coordinates are not negative
            if (start < 0) {
                start = 0
            }

            #print "Start after padding:", start > "/dev/stderr"

            # Output the new, padded coordinates
            print $1, start, end, arr[1]    # Output: Chromosome, NEW Start, NEW End, Set ID (gene name)
        }
    }' "$GTF_file" > "$temp_regions_file"
    #echo "Regions for gene set $set_name saved in $temp_regions_file"

    # Extract SNPs in regions using PLINK
    plink --bfile "${training_file}" \
          --extract range "${temp_regions_file}" \
          --write-snplist \
          --out "${out_SNPs_per_set}/${set_name}_snps"

    #echo "SNP extraction for $set_name completed!"

    # Delete the temporary region file after SNP extraction
    rm -f "$temp_regions_file"

done < "$MSigDB_file"

### Extract rs ID's per gene set
echo "Extracting SNPs per gene set..."
python get_SNPs_per_set_rs.py "$out_SNPs_per_set" "$training_file_rsID"

### Make "Base" files
if [ "$base" == "TRUE" ]; then
    cat "$out_SNPs_per_set"/*_snps.snplist | sort -u > "$out_SNPs_per_set"/Base_snps.snplist
    cat "$out_SNPs_per_set"/*_rs.snplist | sort -u > "$out_SNPs_per_set"/Base_rs.snplist
fi

echo "All SNP extraction completed!"