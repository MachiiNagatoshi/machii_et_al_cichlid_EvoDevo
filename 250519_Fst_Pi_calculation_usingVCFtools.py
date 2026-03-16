# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 250519 vcftoolsによるFstとπの計算

# %% [markdown]
# 入力するvcfの関係でpixyが使えない可能性が出てきたので再度計算を行う。

# %% [markdown]
# ## Fst、piの計算

# %%
import pandas as pd
import subprocess
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
import glob
import allel

def bgzip(filename):
    cmd = f"bgzip -@ 8 {filename}"
    p=subprocess.check_output(cmd, shell=True)

def cp(filename, copiedfile):
    cmd = f"cp -n {filename} {copiedfile}"
    p=subprocess.check_output(cmd, shell=True)

def rm(filename):
    cmd = f"rm {filename}"
    p=subprocess.check_output(cmd, shell=True)

def mkdir(dirname):
    cmd = f"mkdir -p {dirname}"
    p=subprocess.check_output(cmd, shell=True)

def sed(file_path, before, after):
    cmd = f"sed -i \"s/{before}/{after}/g\" {file_path}"
    p=subprocess.check_output(cmd, shell=True)

def calculate_Fst(vcf_path, pop1_path, pop2_path, window_size, window_step, outpath):
    cmd = f"vcftools --gzvcf {vcf_path} --weir-fst-pop {pop1_path} --weir-fst-pop {pop2_path} --fst-window-size {window_size} --fst-window-step {window_step} --out {outpath}"
    p=subprocess.check_output(cmd, shell=True)

def calculate_pi(vcf_path, pop_path, window_size, window_step, outpath):
    cmd = f"vcftools --gzvcf {vcf_path} --keep {pop_path} --window-pi {window_size} --window-pi-step {window_step} --out {outpath}"
    p = subprocess.check_output(cmd, shell=True)

def calculate_tajimaD(vcf_path, pop_path, window_size, outpath):
    cmd = f"vcftools --gzvcf {vcf_path} --keep {pop_path} --TajimaD {window_size} --out {outpath}"
    p = subprocess.check_output(cmd, shell=True)

def calculate_tajimaD_withScikit_allel(vcf_path, pop_path, window_size, window_step, outpath, min_maf=0.0001):
    """Calculate Tajima's D using scikit-allel, processing by chromosome"""
    
    # Read population samples
    with open(pop_path, 'r') as f:
        pop_samples = [line.strip() for line in f if line.strip()]
    
    print(f"  Reading VCF for {len(pop_samples)} samples...")
    
    # Read VCF file
    callset = allel.read_vcf(vcf_path, samples=pop_samples)
    
    # Extract data
    gt = allel.GenotypeArray(callset['calldata/GT'])
    pos = callset['variants/POS']
    chrom = callset['variants/CHROM']
    
    print(f"  Found {len(pos)} variants across chromosomes")
    
    # Calculate allele counts for filtering
    ac_all = gt.count_alleles()
    
    # Filter out invariant sites and apply MAF filter
    # Remove sites where all samples have the same genotype (invariant)
    is_segregating = ac_all.is_segregating()
    
    # Calculate MAF for segregating sites
    af = ac_all.to_frequencies()
    # MAF is the minimum of reference allele freq and alt allele freq
    maf = np.minimum(af[:, 0], 1 - af[:, 0])
    
    # Combined filter: segregating sites with MAF >= threshold
    variant_filter = is_segregating & (maf >= min_maf)
    
    print(f"  After filtering (MAF >= {min_maf}): {np.sum(variant_filter)} variants ({np.sum(variant_filter)/len(pos)*100:.1f}%)")
    
    # Apply filter to all data
    gt = gt[variant_filter]
    pos = pos[variant_filter]
    chrom = chrom[variant_filter]
    
    # Get unique chromosomes
    unique_chroms = np.unique(chrom)
    print(f"  Processing {len(unique_chroms)} chromosomes: {unique_chroms[:5]}...")
    
    all_results = []
    
    # Process each chromosome separately
    for chr_name in unique_chroms:
        chr_mask = chrom == chr_name
        chr_pos = pos[chr_mask]
        chr_gt = gt[chr_mask]
        
        if len(chr_pos) < 3:  # Skip chromosomes with too few variants
            continue
        
        # Calculate allele counts for this chromosome
        chr_ac = chr_gt.count_alleles()
        
        try:
            # Calculate Tajima's D for this chromosome
            D, windows, counts = allel.windowed_tajima_d(
                chr_pos, chr_ac, 
                size=window_size, 
                step=window_step, 
                min_sites=3
            )
            
            # Create results dataframe for this chromosome
            chr_results = pd.DataFrame({
                'CHROM': [chr_name] * len(D),
                'BIN_START': windows[:, 0],
                'BIN_END': windows[:, 1],
                'N_VARIANTS': counts,
                'TajimaD': D
            })
            
            all_results.append(chr_results)
            
        except Exception as e:
            print(f"    Error processing chromosome {chr_name}: {e}")
            continue
    
    # Combine results from all chromosomes
    if all_results:
        results = pd.concat(all_results, ignore_index=True)
        print(f"  Calculated {len(results)} total windows")
    else:
        # Create empty dataframe if no results
        results = pd.DataFrame(columns=['CHROM', 'BIN_START', 'BIN_END', 'N_VARIANTS', 'TajimaD'])
        print("  No windows calculated")
    
    # Save results
    results.to_csv(outpath, sep='\t', index=False)
    return results



# %%
### main

# set variable
WINDOW_SIZE = 10000
STEP_SIZE = 1000

# path
MAIN_DIR = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis"
VCF_PATH = "/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_missing0.9_maf0.05.vcf.recode.vcf.gz"
Fst_LIST_PATH = "/Volumes/research_backup/genome_analysis/240501_ReAnalysis/fst/Fst_list"
OUT_PATH = f"{MAIN_DIR}/vcftools_result_window{WINDOW_SIZE}"

os.chdir(MAIN_DIR)
mkdir(f"{OUT_PATH}")

## calculate Fst
# Hchi vs Hmic
calculate_Fst(VCF_PATH, f"{Fst_LIST_PATH}/Hchi.txt", f"{Fst_LIST_PATH}/Hmic.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/fst_{WINDOW_SIZE}_Hchi_Hmic.txt")
    
# Hchi vs Hsau
calculate_Fst(VCF_PATH, f"{Fst_LIST_PATH}/Hchi.txt", f"{Fst_LIST_PATH}/Hsau.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/fst_{WINDOW_SIZE}_Hchi_Hsau.txt")

# Hchi vs Lruf
calculate_Fst(VCF_PATH, f"{Fst_LIST_PATH}/Hchi.txt", f"{Fst_LIST_PATH}/Lruf.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/fst_{WINDOW_SIZE}_Hchi_Lruf.txt")
    
# Hchi vs Pnye
calculate_Fst(VCF_PATH, f"{Fst_LIST_PATH}/Hchi.txt", f"{Fst_LIST_PATH}/Pnye.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/fst_{WINDOW_SIZE}_Hchi_Pnye.txt")
    
# Hchi vs Ppun
calculate_Fst(VCF_PATH, f"{Fst_LIST_PATH}/Hchi.txt", f"{Fst_LIST_PATH}/Ppun.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/fst_{WINDOW_SIZE}_Hchi_Ppun.txt")

## calculate Pi
# Hchi
calculate_pi(VCF_PATH, f"{Fst_LIST_PATH}/Hchi.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/pi_{WINDOW_SIZE}_Hchi.txt")

# Hmic
calculate_pi(VCF_PATH, f"{Fst_LIST_PATH}/Hmic.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/pi_{WINDOW_SIZE}_Hmic.txt")

# Hsau
calculate_pi(VCF_PATH, f"{Fst_LIST_PATH}/Hsau.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/pi_{WINDOW_SIZE}_Hsau.txt")

# Lruf
calculate_pi(VCF_PATH, f"{Fst_LIST_PATH}/Lruf.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/pi_{WINDOW_SIZE}_Lruf.txt")

# Pnye
calculate_pi(VCF_PATH, f"{Fst_LIST_PATH}/Pnye.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/pi_{WINDOW_SIZE}_Pnye.txt")

# Ppun
calculate_pi(VCF_PATH, f"{Fst_LIST_PATH}/Ppun.txt", WINDOW_SIZE, STEP_SIZE, f"{OUT_PATH}/pi_{WINDOW_SIZE}_Ppun.txt")

# %%
# パス設定
base_path = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/vcftools_result_window10000"

# .fstファイルを圧縮
for file in glob.glob(f"{base_path}/*.fst"):
    subprocess.check_output(f"bgzip -f {file} -@ 128", shell=True)

# .piファイルを圧縮  
for file in glob.glob(f"{base_path}/*.pi"):
    subprocess.check_output(f"bgzip -f {file} -@ 128", shell=True)

# %%
# Set variables
WINDOW_SIZE = 10000
STEP_SIZE = 1000

# Paths
MAIN_DIR = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis"
VCF_PATH = "/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_missing0.9_maf0.05.vcf.recode.vcf.gz"
Fst_LIST_PATH = "/Volumes/research_backup/genome_analysis/240501_ReAnalysis/fst/Fst_list"
OUT_PATH = f"{MAIN_DIR}/tajimad_result_window{WINDOW_SIZE}"

os.chdir(MAIN_DIR)
mkdir(OUT_PATH)

# Calculate Tajima's D for each population
populations = ['Hchi', 'Hmic', 'Hsau', 'Lruf', 'Pnye', 'Ppun']

for pop in populations:
    print(f"Processing {pop}...")
    pop_file = f"{Fst_LIST_PATH}/{pop}.txt"
    output_file = f"{OUT_PATH}/tajimad_{WINDOW_SIZE}_{pop}.txt"
    calculate_tajimaD(VCF_PATH, pop_file, WINDOW_SIZE, output_file)

# %%
# Set variables
WINDOW_SIZE = 10000
STEP_SIZE = 1000

# Paths
MAIN_DIR = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis"
VCF_PATH = "/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_missing0.9_maf0.05.vcf.recode.vcf.gz"
Fst_LIST_PATH = "/Volumes/research_backup/genome_analysis/240501_ReAnalysis/fst/Fst_list"
OUT_PATH = f"{MAIN_DIR}/tajimad_result_window{WINDOW_SIZE}"

os.chdir(MAIN_DIR)
mkdir(OUT_PATH)

# Calculate Tajima's D for each population
populations = ['Hchi', 'Hmic', 'Hsau', 'Lruf', 'Pnye', 'Ppun']

for pop in populations:
    print(f"Processing {pop}...")
    pop_file = f"{Fst_LIST_PATH}/{pop}.txt"
    output_file = f"{OUT_PATH}/tajimad_{WINDOW_SIZE}_{pop}.txt"
    
    try:
        result = calculate_tajimaD_withScikit_allel(VCF_PATH, pop_file, WINDOW_SIZE, STEP_SIZE, output_file, min_maf=0.0001)
        if len(result) > 0:
            valid_D = result['TajimaD'][~np.isnan(result['TajimaD'])]
            print(f"  Windows: {len(result)}, Valid D: {len(valid_D)}, Mean D: {valid_D.mean():.4f}")
        else:
            print(f"  No valid results for {pop}")
    except Exception as e:
        print(f"  Error processing {pop}: {e}")
        continue

print("Completed!")

# %%
# パス設定
base_path = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/tajimad_result_window10000"

# tajimaファイルを圧縮
for file in glob.glob(f"{base_path}/*.D"):
    subprocess.check_output(f"bgzip -f {file} -@ 128", shell=True)


# %%
