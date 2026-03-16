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
# # Fst and π calculation using vcftools

# %%
import subprocess
import os
import glob

def mkdir(dirname):
    cmd = f"mkdir -p {dirname}"
    subprocess.check_output(cmd, shell=True)

def calculate_Fst(vcf_path, pop1_path, pop2_path, window_size, window_step, outpath):
    cmd = f"vcftools --gzvcf {vcf_path} --weir-fst-pop {pop1_path} --weir-fst-pop {pop2_path} --fst-window-size {window_size} --fst-window-step {window_step} --out {outpath}"
    subprocess.check_output(cmd, shell=True)

def calculate_pi(vcf_path, pop_path, window_size, window_step, outpath):
    cmd = f"vcftools --gzvcf {vcf_path} --keep {pop_path} --window-pi {window_size} --window-pi-step {window_step} --out {outpath}"
    subprocess.check_output(cmd, shell=True)

def calculate_tajimaD(vcf_path, pop_path, window_size, outpath):
    cmd = f"vcftools --gzvcf {vcf_path} --keep {pop_path} --TajimaD {window_size} --out {outpath}"
    subprocess.check_output(cmd, shell=True)


# %%
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


## Calculate Tajima's D
# Hchi
calculate_tajimaD(VCF_PATH, f"{Fst_LIST_PATH}/Hchi.txt", WINDOW_SIZE, f"{OUT_PATH}/tajimad_{WINDOW_SIZE}_Hchi.txt")

# Hmic
calculate_tajimaD(VCF_PATH, f"{Fst_LIST_PATH}/Hmic.txt", WINDOW_SIZE, f"{OUT_PATH}/tajimad_{WINDOW_SIZE}_Hmic.txt")

# Hsau
calculate_tajimaD(VCF_PATH, f"{Fst_LIST_PATH}/Hsau.txt", WINDOW_SIZE, f"{OUT_PATH}/tajimad_{WINDOW_SIZE}_Hsau.txt")

# Lruf
calculate_tajimaD(VCF_PATH, f"{Fst_LIST_PATH}/Lruf.txt", WINDOW_SIZE, f"{OUT_PATH}/tajimad_{WINDOW_SIZE}_Lruf.txt")

# Pnye
calculate_tajimaD(VCF_PATH, f"{Fst_LIST_PATH}/Pnye.txt", WINDOW_SIZE, f"{OUT_PATH}/tajimad_{WINDOW_SIZE}_Pnye.txt")

# Ppun
calculate_tajimaD(VCF_PATH, f"{Fst_LIST_PATH}/Ppun.txt", WINDOW_SIZE, f"{OUT_PATH}/tajimad_{WINDOW_SIZE}_Ppun.txt")

# %%
# Set path

# bgzip fst
for file in glob.glob(f"{OUT_PATH}/*.fst"):
    subprocess.check_output(f"bgzip -f {file} -@ 128", shell=True)

# bgzip pi
for file in glob.glob(f"{OUT_PATH}/*.pi"):
    subprocess.check_output(f"bgzip -f {file} -@ 128", shell=True)

# bgzip tajima
for file in glob.glob(f"{OUT_PATH}/*.D"):
    subprocess.check_output(f"bgzip -f {file} -@ 128", shell=True)

