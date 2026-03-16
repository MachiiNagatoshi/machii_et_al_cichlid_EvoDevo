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
# # 250407 Make phyloogenetic tree of candidate regions

# %%
import subprocess
import os

def extract_region(vcf_file, output_prefix, region):
    # 領域の情報を分解
    chrom, start, end = region.split('_')
    
    # 出力ファイル名を生成（領域情報を含む）
    output_name = f"{output_prefix}{chrom}_{start}_{end}"
        
    # VCFToolsコマンドを構築
    cmd_vcftools = f"vcftools --gzvcf {vcf_file} --chr {chrom} --from-bp {start} --to-bp {end} --min-alleles 2 --recode --recode-INFO-all --out {output_name}"
    
    print(f"Executing VCFTools command: {cmd_vcftools}")
    
    # VCFToolsコマンドを実行
    subprocess.run(cmd_vcftools, shell=True, check=True)
    print(f"Extraction completed for region: {region}")
    
    # bgzipコマンドを構築
    cmd_bgzip = f"bgzip -@ 16 {output_name}.recode.vcf"
    
    print(f"Executing bgzip command: {cmd_bgzip}")
    
    # bgzipコマンドを実行
    subprocess.run(cmd_bgzip, shell=True, check=True)
    print(f"Compression completed for region: {region}")
    

# 入力VCFファイルのパス（適宜変更してください）
input_vcf = "/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset2_missing0.8_maf0.01.vcf.recode.vcf.gz"

# 出力ディレクトリ（適宜変更してください）
output_dir = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/make_phylogenitic_tree/vcf_phased/"

# 領域のリスト
regions = ['LG8_16845565_16880833']

# 各領域に対して関数を呼び出す
for region in regions:
    try:
        extract_region(input_vcf, output_dir, region)
    except subprocess.CalledProcessError as e:
        print(f"Error processing region {region}: {e}")

print("All regions processed.")

# %%
import gzip
from collections import defaultdict

def read_vcf(file_path):
    sequences = defaultdict(str)
    sample_names = []
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    with open_func(file_path, mode) as vcf_file:
        for line in vcf_file:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                sample_index = headers.index('FORMAT') + 1
                sample_names = headers[sample_index:]
                continue
            fields = line.strip().split('\t')
            chrom, pos, id, ref, alt, qual, filter, info, format, *samples = fields
            alleles = [ref] + alt.split(',')
            for sample, name in zip(samples, sample_names):
                genotype = sample.split(':')[0]
                allele1, allele2 = genotype.split('|')
                sequences[f"{name}_A"] += alleles[int(allele1)]
                sequences[f"{name}_B"] += alleles[int(allele2)]
    return sequences

def write_phylip(sequences, output_file):
    num_sequences = len(sequences)
    seq_length = len(next(iter(sequences.values())))
    
    with open(output_file, 'w') as phylip_file:
        phylip_file.write(f"{num_sequences} {seq_length}\n")
        for name, seq in sequences.items():
            phylip_file.write(f"{name} {seq}\n")



# %%

# 入力ファイルのリスト
vcf_files = ['LG8_16845565_16880833.recode.vcf.gz']

# 入力ファイルと出力ファイルのベースディレクトリ
base_dir = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/make_phylogenitic_tree/vcf_phased"

# 各ファイルに対して処理を実行
for vcf_file in vcf_files:
    input_path = os.path.join(base_dir, vcf_file)
    output_file = os.path.join(base_dir, vcf_file.replace('.vcf.gz', '.phy'))
    
    sequences = read_vcf(input_path)
    write_phylip(sequences, output_file)
    print(f"変換が完了しました。PHYLIPファイルが {output_file} に書き込まれました。")

print("すべてのファイルの処理が完了しました。")

# %% [markdown]
# Run followinf python file

# %%
# #!/usr/bin/env python3
import subprocess
import os

def convert_vcf_to_phylip(input_vcf):
    output_name = os.path.basename(input_vcf).replace('.recode.vcf.gz', '')
    output_dir = os.path.dirname(input_vcf)
    output_prefix = os.path.join(output_dir, output_name)
    cmd = f"plink2 --vcf {input_vcf} --snps-only --export phylip-phased --out {output_prefix}"
    print(f"Executing command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print(f"Conversion completed for file: {input_vcf}")
    return f"{output_prefix}.phy"

def run_modeltest(input_phylip, modeltest_output_dir):
    output_prefix = os.path.join(modeltest_output_dir, os.path.basename(input_phylip).replace('.phy', ''))
    cmd = f"modeltest-ng -i {input_phylip} -o {output_prefix}_modeltest -d nt -p 16 -T raxml"
    print(f"Executing command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print(f"Modeltest completed for file: {input_phylip}")
    return f"{output_prefix}_modeltest.log"

def parse_modeltest_output(modeltest_output):
    best_model = None
    with open(modeltest_output, 'r') as f:
        for line in f:
            if line.strip().startswith("> raxml-ng"):
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == "--model":
                        if i + 1 < len(parts):
                            best_model = parts[i + 1]
                            print(f"Best model found: {best_model}")
                            return best_model
    
    if best_model is None:
        print(f"Failed to find best model in {modeltest_output}")
    
    return best_model

def run_raxml(input_phylip, model, raxml_output_dir):
    output_prefix = os.path.join(raxml_output_dir, os.path.basename(input_phylip).replace('.phy', ''))
    cmd = f"raxml-ng --msa {input_phylip} --model {model} --prefix {output_prefix}_raxml --bs-trees 1000 --all --threads 16 --seed 12345"
    print(f"Executing command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print(f"RAxML-NG analysis completed for file: {input_phylip}")

def process_region(args):
    region, input_dir, modeltest_output_dir, raxml_output_dir = args
    phylip_file = os.path.join(input_dir, f"{region}.recode.phy")
    if not os.path.exists(phylip_file):
        print(f"VCF file not found for region: {region}")
        return

    try:
        
        # Run ModelTest-NG
        modeltest_output = run_modeltest(phylip_file, modeltest_output_dir)
        
        # Parse ModelTest-NG output
        best_model = parse_modeltest_output(modeltest_output)
        if best_model is None:
            print(f"Failed to determine best model for region: {region}")
            print(f"Please check the ModelTest-NG output file: {modeltest_output}")
            return
        
        print(f"Best model for region {region}: {best_model}")
        
        # Run RAxML-NG
        run_raxml(phylip_file, best_model, raxml_output_dir)
    except subprocess.CalledProcessError as e:
        print(f"Error processing region {region}: {e}")
        print(f"Command that failed: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.output}")
    except Exception as e:
        print(f"Unexpected error processing region {region}: {e}")
        import traceback
        traceback.print_exc()

def main():
    input_dir = "/home/nakaharu/HDD5/machii/250927_make_phylogenetic_tree/vcf_phased"
    modeltest_output_dir = "/home/nakaharu/HDD5/machii/250927_make_phylogenetic_tree/modeltest"
    raxml_output_dir = "/home/nakaharu/HDD5/machii/250927_make_phylogenetic_tree/raxml"

    os.makedirs(modeltest_output_dir, exist_ok=True)
    os.makedirs(raxml_output_dir, exist_ok=True)

    regions = ['LG8_16845565_16880833']

    for region in regions:
        args = (region, input_dir, modeltest_output_dir, raxml_output_dir)
        process_region(args)
    print("All regions processed.")

if __name__ == "__main__":
    main()
