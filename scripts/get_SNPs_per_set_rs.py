import os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from multiprocessing import Pool
from performance import get_time, data
from time import time

# --------------------- 1. Load BIM File ---------------------
def load_bim_file(training_data):
    """Loads the BIM file and returns it as a Pandas DataFrame."""
    print(f"Loading BIM file: {training_data}.bim", flush=True)
    
    t0 = time()
    bim_file = pd.read_csv(
        f"{training_data}.bim", sep='\s+', header=None, usecols=[0, 1, 3], 
        names=['chrom', 'snp_id', 'pos'], dtype={'chrom': np.int32, 'snp_id': str, 'pos': np.int32}
    )
    print("reading bim file: ", time() - t0)
    to = time()
    bim_dict = {(row.chrom, row.pos): row.snp_id for row in bim_file.itertuples(index=False)}
    print(f"Loaded {len(bim_dict)} SNPs from BIM file.", flush=True)
    print("second thingy: ", time() - t0)

    return bim_dict  # Store SNPs in a dictionary for O(1) lookup

# --------------------- 2. Extract SNP IDs ---------------------
@get_time
def extract_snp_ids(chrom, pos, bim_dict):
    """Extracts SNP ID from preloaded BIM dictionary."""
    return bim_dict.get((chrom, pos), None)  # O(1) lookup time

# --------------------- 4. Process Each Gene Set ---------------------
@get_time
def process_gene_set(snp_file, bim_dict):
    t0 = time()
    """Processes a single gene set."""
    set_name = os.path.basename(snp_file).replace("_snps.snplist", "")
    print(f"\nProcessing gene set: {set_name}", flush=True)

    snp_data = {chrom: [] for chrom in range(1, 23)}  # A dictionary to hold SNPs for each chromosome
    snp_ids = []  # List to store all the SNP rsIDs for the pathway file

    with open(snp_file, 'r') as infile:
        snp_count = 0
        for line in infile:
            snp_pos = line.strip()
            try:
                chrom, pos = map(int, snp_pos.split(":"))  # Extract chromosome and position
                snp_id = extract_snp_ids(chrom, pos, bim_dict)  # Extract SNP ID
                if snp_id:
                    snp_data[chrom].append((snp_id, pos))  # Collect SNPs for each chromosome
                    snp_ids.append(snp_id)  # Add SNP rsID to the list for this gene set
                    snp_count += 1
            except Exception as e:
                print(f"Error processing SNP {snp_pos}: {e}", flush=True)

    print(f"Extracted {snp_count} SNPs for gene set {set_name}.", flush=True)
    data.show()
    print("Done reading SNP file:", time() - t0, flush=True)

    # Write the list of SNP rsIDs to a file (pathwayname_rs.snplist)
    snp_ids_file = os.path.join(out_snps_per_set, f"{set_name}_rs.snplist")
    with open(snp_ids_file, 'w') as snp_outfile:
        for snp_id in snp_ids:
            snp_outfile.write(f"{snp_id}\n")

# --------------------- 6. Main Function ---------------------
def main(out_snps_per_set, training_data):
    print("Starting SNP extraction and LD computation...", flush=True)
    
    # Load BIM data once and pass as dictionary for fast lookup
    bim_dict = load_bim_file(training_data)

    snp_files = [os.path.join(out_snps_per_set, f) for f in os.listdir(out_snps_per_set) if f.endswith('_snps.snplist')]
    print(f"Found {len(snp_files)} gene sets to process.", flush=True)

    # Parallel execution of gene set processing
    Parallel(n_jobs= len(snp_files))(delayed(process_gene_set)(snp_file, bim_dict) for snp_file in snp_files)

# --------------------- 7. Entry Point ---------------------
if __name__ == "__main__":
    import sys
    out_snps_per_set = sys.argv[1]
    training_data = sys.argv[2]
    
    main(out_snps_per_set, training_data)