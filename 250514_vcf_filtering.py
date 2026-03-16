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
# # 240604 VCF filtering and something else

# %% [markdown]
# ## 1. Compare coverage and mapping rate
#
# We reanalyzed mosdepth and bamtools output processed by multiqc

# %%
# #! /bin/bash

sample_list=(
    HC_2020_1
    HC_2020_3
    HC_2020_4
    HC_2021_1
    HC_2021_10
    HC_2021_11
    HC_2021_12
    HC_2021_2
    HC_2021_3
    HC_2021_4
    HC_2021_5
    HC_2021_6
    HC_2021_7
    HC_2021_8
    HC_2021_9
    SAMD00269291
    SAMD00269292
    SAMD00269293
    SAMD00269294
    SAMD00269295
    SAMD00269296
    SRR12394623
    SRR12700909
    HMic_2021_1
    HMic_2022_2_1
    HMic_2021_2
    HMic_2022_2_2
    HMic_2021_3
    HMic_2022_2_3
    HMic_2022_2_4
    SAMD00269297
    SAMD00269298
    SAMD00269299
    SAMD00269300
    SAMD00269301
    SAMD00269302
    SRR12394583
    SRR12700891
    HMh_2020_6
    HMh_2020_1
    HMh_2020_2
    HMh_2022_2_1
    HMh_2022_2_2
    SAMD00269303
    SAMD00269304
    SAMD00269305
    SAMD00269306
    SAMD00269307
    SAMD00269308
    PN_2022_2_1
    PN_2022_2_2
    PN_2022_2_3
    SRR4169621
    SRR4169622
    SRR4169623
    SRR4169625
    SRR12394594
    SRR12700883
    SRR4169617
    SRR4169618
    SRR4169619
    SRR4169633
    SRR12700840
    PBB_2022_2_1
    PBB_2022_2_2
    PBB_2022_2_3
    PBB_2022_2_4
    PBB_2022_2_5
    PBB_2022_2_6
    PBB_2022_2_7
    ERR3469453
    ERR3469455
    ERR3469621
    ERR3469499
    ERR3469504
    ERR3469466
    ERR3469509
    ERR3469517
    ERR3469524
    ERR3479817
    ERR3479819
    ERR3479820
    ERR3479830
    ERR3479895
    ERR3479850
    ERR3479785
    ERR3469473
    ERR3469640
    ERR3469644
    ERR3469485
    Llab01
    Llab02
    Llab03
    Llab04
    Llab05
    Llab06
    Llab07
    Llab08
    Llab09
    Llab10
    SRR9665643
    SRR9665640
    SRR7662386
    SRR7662400
    SRR7662397
    SRR7662507
    SRR7662509
    SRR7662530
    SRR7662422
    SRR7662504
    SRR7662503
    Ldar01
    Ldar02
    Ldar03
    Ldar04
    Ldar05
    SRR9657565
    SRR9657566
    SAMN00139653
    SRR9673979
    SRR9673986
    SRR9657493
    SRR9657510
    SRR9657541
    SRR9657544
    SRR9674047
    SRR9674048
    SRR9674000
    SRR9657526
    SRR9665680
    SRR9665704
    SRR9675360
    SRR9675363
    SRR9657561
    SRR9657564
    SRR9673877
    SRR9673941
    SRR9673898
    SRR9673901
    SRR9673878
    SRR9673879
    SRR9665702
    SRR9665705
    SRR9657496
    SRR9657517
    SRR9675322
    SRR9675321
    SRR9675371
    SRR9675370
    SRR9673846
    SRR9673847
    Tmoo01
    Tmoo02
    Tmoo03
    Tmoo04
    Tmoo05
    SRR9673968
    SRR9674025
    SRR9673910
    SRR9673913
    SRR9673959
    SRR9673964
    SRR9673833
    SRR9673896
    SRR9673884
    SRR9673885
    SRR9673922
    SRR9673923
    SRR9673918
    SRR9673919
    SRR9673813
    SRR9673955
    SRR9673820
    SRR9673821
    SRR9665655
    SRR9665658
    SRR17262985
    ERR1749533
    SAMEA1877480
    ERR1749477
    ERR1749479
    ERR1749481
    ERR1749483
    ERR1749485
    ERR1749487
    ERR1749476
    ERR1749478
    ERR1749480
    ERR1749482
    ERR1749484
    ERR1749486
    ERR715538
    ERR715539
    ERR715540
    ERR1749539
    ERR1749541
    ERR1749540
    ERR1749542
    SRR17068917
    SRR17068915
    SRR17068914
    ERR1749503
    ERR1749511
    ERR1749513
    ERR1749512
    ERR1749514
    SAMEA1877409
    SAMEA1877451
    ERR1749421
    ERR1749423
    ERR1749422
    ERR715514
    ERR1081372
    SAMEA1877469
    ERR1749491
    ERR1749497
    ERR1749501
    ERR1749490
    ERR1749496
    ERR1749498
    ERR1749502
)

# path
BeforeMV_PATH=/Volumes/research_backup/230220_VariantCall_results
AfterMV_PATH=/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/quality_check

# mkdir '${AfterMV_PATH}'/mosdepth
# mkdir '${AfterMV_PATH}'/bamtools_stats

# cd '${BeforeMV_PATH}'

for fname in '${sample_list[@]}'
do
    if [ -d ${fname} ]; then
        cp -n '${fname}'/mosdepth/*.txt '${AfterMV_PATH}'/mosdepth
        cp -n '${fname}'/bamtools_stats/*_bamstat '${AfterMV_PATH}'/bamtools_stats
    fi
done

# cd '${AfterMV_PATH}'/mosdepth
multiqc .

# cd '${AfterMV_PATH}'/bamtools_stats
multiqc .

# %% [markdown]
# ## Visualization

# %%
#　python
# import library
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %%
### read file

# set path
bamtools_stats_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/quality_check/multiqc_bamtools_stats.txt'
mosdepth_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/quality_check/multiqc_mosdepth_general_stats.txt'
label_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/quality_check/sample_label.txt'


# read tsv
df_bamtools_stats = pd.read_csv(bamtools_stats_path, sep='\t')
df_mosdepth = pd.read_csv(mosdepth_path, sep='\t')
df_label = pd.read_csv(label_path, sep='\t')


print(df_label)
print(df_bamtools_stats)

# %%
### merge datas

# merge
df_mapping_quality = pd.merge(df_label, df_bamtools_stats, how='left', left_on='Run', right_on='Sample')
df_mapping_quality = pd.merge(df_mapping_quality, df_mosdepth, how='left', left_on='Run', right_on='Sample')

df_mapping_quality = df_mapping_quality.dropna(subset='total_reads')
df_mapping_quality = df_mapping_quality.drop_duplicates()

df_mapping_quality.describe()

# %%
# plot
fig = plt.subplots(figsize=(40,5), dpi=100)
sns.swarmplot(x='species_id', y='mapped_reads_pct', data=df_mapping_quality[['species_id', 'mapped_reads_pct']], size=4, edgecolor='black', linewidth=0.5)
plt.xticks(rotation=-45, fontsize=18)
plt.show()


# plot
fig = plt.subplots(figsize=(40,5), dpi=100)
sns.swarmplot(x='species_id', y='total_reads', data=df_mapping_quality[['species_id', 'total_reads']], size=4, edgecolor='black', linewidth=0.5)
plt.xticks(rotation=-45, fontsize=18)
plt.show()

# plot
fig = plt.subplots(figsize=(40,5), dpi=100)
sns.swarmplot(x='species_id', y='mosdepth-mean_coverage', data=df_mapping_quality[['species_id', 'mosdepth-mean_coverage']], size=4, edgecolor='black', linewidth=0.5)
plt.xticks(rotation=-45, fontsize=18)
plt.show()

# %% [markdown]
# ##  Filetering each VCF with depth

# %%
### import library and set functions

import pandas as pd
import subprocess
import os
import seaborn as sns
import matplotlib.pyplot as plt
import glob

def bgzip(filename):
    cmd = f'bgzip -@ 8 {filename}'
    p=subprocess.check_output(cmd, shell=True)

def cp(filename, copiedfile):
    cmd = f'cp -n {filename} {copiedfile}'
    p=subprocess.check_output(cmd, shell=True)

def rm(filename):
    cmd = f'rm {filename}'
    p=subprocess.check_output(cmd, shell=True)

def mkdir(dirname):
    cmd = f'mkdir -p {dirname}'
    p=subprocess.check_output(cmd, shell=True)

def tabix(vcfpath):
    cmd = f'tabix -p vcf {vcfpath}'
    p=subprocess.check_output(cmd, shell=True)

def bcftools_index(vcfpath):
    cmd = f'bcftools index {vcfpath}'
    p=subprocess.check_output(cmd, shell=True)

def bcftools_concat(snppath, indelpath, outpath):
    cmd = f'bcftools concat -a -d all -o {outpath} -O v {snppath} {indelpath}'
    p=subprocess.check_output(cmd, shell=True)

def bcftools_consensus(vcfpath, refpath, outpath):
    cmd = f'bcftools consensus --include \'FILTER=\'PASS\'\' --fasta-ref {refpath} -H I --output {outpath} {vcfpath}'
    # -H I (マージしたbamがない際にはこのオプションをつける)
    p=subprocess.check_output(cmd, shell=True)

def SelectVariants(vcfpath, outpath):
    cmd = f'gatk SelectVariants --exclude-filtered -V {vcfpath} -O {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def VariantFiltration(vcfpath, outpath, refgenome):
    # calculate cut off threshold
    min_depth = 4.0
    max_depth = 100.0

    # filterがPASSじゃない行を除く
    cmd = f'gatk SelectVariants --exclude-filtered -V {vcfpath} -O {vcfpath}_onlyPASS.vcf'
    p=subprocess.check_output(cmd, shell=True)
    
    # depthの低いサイトをミッシングにする
    cmd = f'gatk VariantFiltration -R {refgenome} -V {vcfpath}_onlyPASS.vcf -O {vcfpath}_cutdepth.vcf ' + \
          f'--genotype-filter-expression 'DP < {min_depth}' --genotype-filter-name 'minDP' ' + \
          f'--genotype-filter-expression 'DP > {max_depth}' --genotype-filter-name 'maxDP''
    p=subprocess.check_output(cmd, shell=True)

    cmd = f'gatk SelectVariants --set-filtered-gt-to-nocall -V {vcfpath}_cutdepth.vcf -O {outpath}'
    p=subprocess.check_output(cmd, shell=True)

    cmd = f'rm {vcfpath}_onlyPASS.vcf'
    p=subprocess.check_output(cmd, shell=True)

    cmd = f'rm {vcfpath}_cutdepth.vcf'
    p=subprocess.check_output(cmd, shell=True)
    
    cmd = f'rm {vcfpath}_cutdepth.vcf.idx'
    p=subprocess.check_output(cmd, shell=True)
    
def samtools_faidx(genome_path, pos, outpath):
    cmd = f'samtools faidx {genome_path} \'{pos}\' > {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def cat_fasta(fasta_dir, outpath):
    cmd = f'cat {fasta_dir}/*fasta > {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def mafft(fasta_dir, outpath):
    cmd = f'mafft --auto --thread 16 {fasta_dir} > {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def bcftools_view(vcf_path, pos, outpath):
    cmd = f'bcftools view -r {pos} {vcf_path} > {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def sed(file_path, before, after):
    cmd = f'sed -i \'s/{before}/{after}/g\' {file_path}'
    p=subprocess.check_output(cmd, shell=True)

def calculate_Fst(vcf_path, pop1_path, pop2_path, window_size, window_step, outpath):
    cmd = f'vcftools --gzvcf {vcf_path} --weir-fst-pop {pop1_path} --weir-fst-pop {pop2_path} --fst-window-size {window_size} --fst-window-step {window_step} --out {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def bcftools_merge(vcf_dir, outpath):
    cmd = f'bcftools merge --merge all --force-samples --threads 8 --no-index {vcf_dir}/*.vcf.gz > {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def vcftools_filter_missing(vcf_path, max_missing, outpath):
    cmd = f'vcftools --gzvcf {vcf_path} --remove-indels --recode --max-missing {max_missing} --out {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def calculate_target_range(target_range, target_diff):
    # split target range
    chromosome = target_range.split(':')[0]
    start_pos = int(target_range.split(':')[1].split('-')[0])
    end_pos = int(target_range.split(':')[1].split('-')[1])
    # calculate target range
    start_pos -= target_diff
    end_pos += target_diff
    # make new target_range
    target_range_new = f'{chromosome}:{start_pos}-{end_pos}'
    return target_range_new

def vcftools_remove_indels_missing(vcf_path, outpath, max_missing):
    cmd = f'vcftools --gzvcf {vcf_path} --remove-indels --recode --max-missing {max_missing} --out {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def vcftools_filter_maf(vcf_path, outpath, maf):
    cmd = f'vcftools --gzvcf {vcf_path} --recode --maf {maf} --out {outpath}'
    p=subprocess.check_output(cmd, shell=True)

def pca(vcf,out):
    cmd = f'plink --vcf {vcf} --allow-extra-chr --const-fid --out {out} --pca'
# --const-fid: 全サンプルに同一のダミーFID (0) を割り当て、VCFのサンプルIDをIIDとして使用。家系情報がない場合や区別不要時に用いる。
    p=subprocess.check_output(cmd, shell=True)


# %%
### set list and path

# path
RAW_PATH = '/Volumes/research_backup/230220_VariantCall_results'
MAIN_DIR = '/Volumes/Transcend/250514_vcf_filtered'
FILTER_1_DIR = '/Volumes/Transcend/250514_vcf_filtered/vcf_filtered1'
FILTER_2_DIR = '/Volumes/Transcend/250514_vcf_filtered/vcf_filtered2'
REF_GENOME = '/Volumes/research_backup/genome/Genome4VariantCall/M_zebra/Maylandia_zebra.M_zebra_UMD2a.dna_sm.toplevel.masked.fa'
sample_list = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/quality_check/dataset2.txt'

mkdir(MAIN_DIR)
mkdir(FILTER_1_DIR)
mkdir(FILTER_2_DIR)
os.chdir(MAIN_DIR)

# read average depth info
df_sample_list = pd.read_csv(sample_list, sep='\t')
df_sample_list = df_sample_list.set_index('sample_id')

sample_list = df_sample_list['Run'].to_list()

df_sample_list

# %%
### filter vcf with depth

for sample in sample_list:
    # check directory and file existence
    if os.path.isdir(f'{RAW_PATH}/{sample}/'):
 
        # SNPsファイルが存在し、対応するFILTER_1のファイルが存在しないか確認
        if len(glob.glob(f'{RAW_PATH}/{sample}/vcf/{sample}_*_snps_filtered.vcf.gz'))>0 and len(glob.glob(f'{FILTER_1_DIR}/{sample}_*_snps_filtered.vcf.gz'))==0:
 
            # copy file
            cp(f'{RAW_PATH}/{sample}/vcf/{sample}_*_snps_filtered.vcf.gz', f'{FILTER_1_DIR}')
            cp(f'{RAW_PATH}/{sample}/vcf/{sample}_*_indels_filtered.vcf.gz', f'{FILTER_2_DIR}')
    
            # make index
            tabix(f'{FILTER_1_DIR}/{sample}_*_snps_filtered.vcf.gz')
            
        if len(glob.glob(f'{FILTER_1_DIR}/{sample}_*_snps_filtered.vcf.gz'))>0 and len(glob.glob(f'{FILTER_2_DIR}/{sample}_snps_filtered_2.vcf.gz'))==0:
            # filter vcf with depth
            VariantFiltration(f'{FILTER_1_DIR}/{sample}_*_snps_filtered.vcf.gz', f'{FILTER_2_DIR}/{sample}_snps_filtered_2.vcf', REF_GENOME)

            # compress vcf file
            bgzip(f'{FILTER_2_DIR}/{sample}_snps_filtered_2.vcf')

            # make bcftools index
            bcftools_index(f'{FILTER_2_DIR}/{sample}_snps_filtered_2.vcf.gz')

            # remove needless index
            rm(f'{FILTER_2_DIR}/{sample}_snps_filtered_2.vcf.idx')


# %% [markdown]
# ## merge vcf

# %%
## merge vcf
def vcftools_merge(vcf_dir, outpath, sample_list=None):

    if sample_list:
        # サンプルリストが提供された場合、指定されたサンプルのVCFファイルパスを構築
        vcf_files = [f'{vcf_dir}/{sample}_snps_filtered_2.vcf.gz' for sample in sample_list]
        
        # 存在チェック（オプション）
        for vcf_file in vcf_files:
            if not os.path.exists(vcf_file):
                print(f'警告: ファイル {vcf_file} が見つかりません')
        
        # 存在するファイルのみをスペース区切りで連結
        vcf_files_str = ' '.join([f for f in vcf_files if os.path.exists(f)])
        
        # ファイルが1つもない場合はエラー
        if not vcf_files_str:
            raise FileNotFoundError('指定されたサンプルのVCFファイルが見つかりません')
    else:
        # サンプルリストが提供されない場合は、すべてのVCFファイルを使用
        vcf_files_str = f'{vcf_dir}/*_snps_filtered_2.vcf.gz'
    
    # マージコマンドを構築して実行
    cmd = f'vcf-merge -s -R 0/0 {vcf_files_str} | bgzip -c > {outpath}'
    print(cmd)
    
    # サブプロセスを実行し、エラー処理を追加
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        print(f'マージが完了しました。出力ファイル: {outpath}')
        return result
    except subprocess.CalledProcessError as e:
        print(f'エラー: マージコマンドの実行に失敗しました: {e}')
        print(f'エラーメッセージ: {e.output.decode() if hasattr(e, 'output') else '不明'}')
        raise

# make directory
VCF_MERGE_DIR = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged'
mkdir(VCF_MERGE_DIR)

# read sample list
dataset1_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/quality_check/dataset1.txt'
df_dataset1 = pd.read_csv(dataset1_path, sep='\t')
dataset1_sample_list = df_dataset1['Run'].to_list()

# merge vcf using vcf-merge in vcftools
vcftools_merge(FILTER_2_DIR, VCF_MERGE_DIR + '/vcf_dataset1_merged.vcf.gz', dataset1_sample_list)


# %%
# read sample list
dataset2_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/quality_check/dataset2.txt'
df_dataset2 = pd.read_csv(dataset2_path, sep='\t')
dataset2_sample_list = df_dataset2['Run'].to_list()

# merge vcf using vcf-merge in vcftools
vcftools_merge(FILTER_2_DIR, VCF_MERGE_DIR + '/vcf_dataset2_merged.vcf.gz', dataset2_sample_list)

# %% [markdown]
# ## filtering missing sites

# %%
### merge vcf
## make directory
VCF_MERGE_dataset1_PATH = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_merged.vcf.gz'
VCF_MERGE_dataset2_PATH = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset2_merged.vcf.gz'

VCF_MERGE_DIR = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged'

## remove indels and missing sites

vcftools_remove_indels_missing(VCF_MERGE_dataset1_PATH, VCF_MERGE_DIR + '/vcf_merged_dataset1_missing' + str(1) + '.vcf', 1)
vcftools_remove_indels_missing(VCF_MERGE_dataset2_PATH, VCF_MERGE_DIR + '/vcf_merged_dataset2_missing' + str(0.8) + '.vcf', 0.8)

# compress vcf file
bgzip(VCF_MERGE_DIR + '/vcf_merged_dataset1_missing' + str(1) + '.vcf.recode.vcf')
bgzip(VCF_MERGE_DIR + '/vcf_merged_dataset2_missing' + str(0.8) + '.vcf.recode.vcf')

# make index
tabix(VCF_MERGE_DIR + '/vcf_merged_dataset1_missing' + str(1) + '.vcf.recode.vcf.gz')
tabix(VCF_MERGE_DIR + '/vcf_merged_dataset2_missing' + str(0.8) + '.vcf.recode.vcf.gz')


# %%
### merge vcf
## make directory
VCF_MERGE_dataset1_PATH = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_merged.vcf.gz'
VCF_MERGE_dataset2_PATH = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset2_merged.vcf.gz'

VCF_MERGE_DIR = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged'

## remove indels and missing sites

vcftools_remove_indels_missing(VCF_MERGE_dataset1_PATH, VCF_MERGE_DIR + '/vcf_merged_dataset1_missing' + str(0.9) + '.vcf', 0.9)
vcftools_remove_indels_missing(VCF_MERGE_dataset2_PATH, VCF_MERGE_DIR + '/vcf_merged_dataset2_missing' + str(0.8) + '.vcf', 0.8)

# compress vcf file
bgzip(VCF_MERGE_DIR + '/vcf_merged_dataset1_missing' + str(0.9) + '.vcf.recode.vcf')
bgzip(VCF_MERGE_DIR + '/vcf_merged_dataset2_missing' + str(0.8) + '.vcf.recode.vcf')

# make index
tabix(VCF_MERGE_DIR + '/vcf_merged_dataset1_missing' + str(0.9) + '.vcf.recode.vcf.gz')
tabix(VCF_MERGE_DIR + '/vcf_merged_dataset2_missing' + str(0.8) + '.vcf.recode.vcf.gz')


# %% [markdown]
# ## filtring maf

# %%
### merge vcf
## make directory
VCF_MERGE_dataset1_filtermissing_PATH = VCF_MERGE_DIR + '/vcf_merged_dataset1_missing' + str(0.9) + '.vcf.recode.vcf.gz'
VCF_MERGE_dataset2_filtermissing_PATH = VCF_MERGE_DIR + '/vcf_merged_dataset2_missing' + str(0.8) + '.vcf.recode.vcf.gz'

maf_threshold_dataset1 = 0.05
maf_threshold_dataset2 = 0.01

vcftools_filter_maf(VCF_MERGE_dataset1_filtermissing_PATH, VCF_MERGE_DIR + '/vcf_dataset1_missing' + str(0.9) + '_maf' + str(maf_threshold_dataset1) + '.vcf', maf_threshold_dataset1)
vcftools_filter_maf(VCF_MERGE_dataset2_filtermissing_PATH, VCF_MERGE_DIR + '/vcf_dataset2_missing' + str(0.8) + '_maf' + str(maf_threshold_dataset2) + '.vcf', maf_threshold_dataset2)

# compress vcf file
bgzip(VCF_MERGE_DIR + '/vcf_dataset1_missing' + str(0.9) + '_maf' + str(maf_threshold_dataset1) + '.vcf.recode.vcf')
bgzip(VCF_MERGE_DIR + '/vcf_dataset2_missing' + str(0.8) + '_maf' + str(maf_threshold_dataset2) + '.vcf.recode.vcf')

# make index
tabix(VCF_MERGE_DIR + '/vcf_dataset1_missing' + str(0.9) + '_maf' + str(maf_threshold_dataset1) + '.vcf.recode.vcf.gz')
tabix(VCF_MERGE_DIR + '/vcf_dataset2_missing' + str(0.8) + '_maf' + str(maf_threshold_dataset2) + '.vcf.recode.vcf.gz')


# %% [markdown]
# ## phasing

# %% [markdown]
# ### install beagle

# %%
(base) ~ ❯❯❯ conda install -c bioconda beagle
Retrieving notices: done
Channels:
 - bioconda
 - conda-forge
Platform: osx-64
Collecting package metadata (repodata.json): done
Solving environment: done


==> WARNING: A newer version of conda exists. <==
    current version: 24.11.3
    latest version: 25.3.1

Please update conda by running

    $ conda update -n base -c conda-forge conda



## Package Plan ##

  environment location: /Users/machiinagatoshi/mambaforge

  added / updated specs:
    - beagle


The following packages will be downloaded:

    package                    |            build
    ---------------------------|-----------------
    _python_abi3_support-1.0   |       hd8ed1ab_2           8 KB  conda-forge
    beagle-5.4_29Oct24.c8e     |       hdfd78af_0         317 KB  bioconda
    certifi-2025.4.26          |     pyhd8ed1ab_0         154 KB  conda-forge
    cpython-3.10.17            |  py310hd8ed1ab_0          49 KB  conda-forge
    dask-2025.5.0              |     pyhd8ed1ab_0           8 KB  conda-forge
    dask-core-2025.5.0         |     pyhd8ed1ab_0         970 KB  conda-forge
    distributed-2025.5.0       |     pyhd8ed1ab_0         782 KB  conda-forge
    dna_features_viewer-3.1.5  |     pyh7e72e81_0          31 KB  bioconda
    openjdk-22.0.1             |       h2d185b6_0       168.9 MB  conda-forge
    polars-1.29.0              |default_ha02c4c5_1           7 KB  conda-forge
    polars-default-1.29.0      |   py39he27f111_1        27.4 MB  conda-forge
    python-gil-3.10.17         |       hd8ed1ab_0          49 KB  conda-forge
    raxml-ng-1.2.2             |       hb53313b_2         785 KB  bioconda
    ------------------------------------------------------------
                                           Total:       199.4 MB

The following NEW packages will be INSTALLED:

  _python_abi3_supp~ conda-forge/noarch::_python_abi3_support-1.0-hd8ed1ab_2 
  beagle             bioconda/noarch::beagle-5.4_29Oct24.c8e-hdfd78af_0 
  cpython            conda-forge/noarch::cpython-3.10.17-py310hd8ed1ab_0 
  openjdk            conda-forge/osx-64::openjdk-22.0.1-h2d185b6_0 
  polars-default     conda-forge/osx-64::polars-default-1.29.0-py39he27f111_1 
  python-gil         conda-forge/noarch::python-gil-3.10.17-hd8ed1ab_0 

The following packages will be UPDATED:

  ca-certificates    conda-forge/osx-64::ca-certificates-2~ --> conda-forge/noarch::ca-certificates-2025.4.26-hbd8a1cb_0 
  certifi                            2025.1.31-pyhd8ed1ab_0 --> 2025.4.26-pyhd8ed1ab_0 
  dask                                2025.3.0-pyhd8ed1ab_0 --> 2025.5.0-pyhd8ed1ab_0 
  dask-core                           2025.3.0-pyhd8ed1ab_0 --> 2025.5.0-pyhd8ed1ab_0 
  distributed                         2025.3.0-pyhd8ed1ab_0 --> 2025.5.0-pyhd8ed1ab_0 
  dna_features_view~                     3.1.3-pyh7cba7a3_0 --> 3.1.5-pyh7e72e81_0 
  openssl                                  3.4.1-hc426f3f_0 --> 3.5.0-hc426f3f_1 
  polars                             1.26.0-py310hc09f4ba_0 --> 1.29.0-default_ha02c4c5_1 
  raxml-ng                                 1.2.2-hb53313b_1 --> 1.2.2-hb53313b_2 


Proceed ([y]/n)? y


Downloading and Extracting Packages:
                                                                                                                                                                                     
Preparing transaction: done                                                                                                                                                          
Verifying transaction: done                                                                                                                                                          
Executing transaction: done                                                                                                                                                          
(base) ~ ❯❯❯ beagle
beagle.29Oct24.c8e.jar (version 5.4)
Copyright (C) 2014-2022 Brian L. Browning
Usage: java -jar beagle.29Oct24.c8e.jar [arguments]

data parameters ...
  gt=<VCF file with GT FORMAT field>                 (required)
  ref=<bref3 or VCF file with phased genotypes>      (optional)
  out=<output file prefix>                           (required)
  map=<PLINK map file with cM units>                 (optional)
  chrom=<[chrom] or [chrom]:[start]-[end]>           (optional)
  excludesamples=<file with 1 sample ID per line>    (optional)
  excludemarkers=<file with 1 marker ID per line>    (optional)

phasing parameters ...
  burnin=<max burnin iterations>                     (default=3)
  iterations=<phasing iterations>                    (default=12)
  phase-states=<model states for phasing>            (default=280)

imputation parameters ...
  impute=<impute ungenotyped markers (true/false)>   (default=true)
  imp-states=<model states for imputation>           (default=1600)
  cluster=<max cM in a marker cluster>               (default=0.005)
  ap=<print posterior allele probabilities>          (default=false)
  gp=<print posterior genotype probabilities>        (default=false)

general parameters ...
  ne=<effective population size>                     (default=100000)
  err=<allele mismatch probability>                  (default: data dependent)
  em=<estimate ne and err parameters (true/false)>   (default=true)
  window=<window length in cM>                       (default=40.0)
  overlap=<window overlap in cM>                     (default=2.0)
  seed=<random seed>                                 (default=-99999)
  nthreads=<number of threads>                       (default: machine dependent)

# %% [markdown]
# ### Execution

# %%
# #!/bin/bash

# Set paht
INPUT_VCF='/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_dataset1_missing0.9_maf0.05.vcf.recode.vcf.gz'
OUTPUT_DIR='/Volumes/Transcend/250514_vcf_filtered/vcf_phased/'
FINAL_OUTPUT='/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05.vcf.recode.vcf.gz'

CHROMOSOMES=(
  'LG1' 'LG2' 'LG3' 'LG4' 'LG5' 'LG6' 'LG7' 'LG8' 'LG9' 'LG10' 
  'LG11' 'LG12' 'LG13' 'LG14' 'LG15' 'LG16' 'LG17' 'LG18' 'LG19' 'LG20' 
  'LG22' 'LG23' 'AGTA05000023.1' 'AGTA05000024.1' 'AGTA05000025.1' 'AGTA05000026.1' 
  'AGTA05000027.1' 'AGTA05000028.1' 'AGTA05000029.1' 'AGTA05000030.1' 'AGTA05000031.1' 
  'AGTA05000032.1' 'AGTA05000033.1' 'AGTA05000034.1' 'AGTA05000035.1' 'AGTA05000036.1' 
  'AGTA05000037.1' 'AGTA05000038.1' 'AGTA05000039.1' 'AGTA05000040.1' 'AGTA05000041.1' 
  'AGTA05000042.1' 'AGTA05000043.1' 'AGTA05000044.1' 'AGTA05000045.1' 'AGTA05000046.1' 
  'AGTA05000047.1' 'AGTA05000048.1' 'AGTA05000049.1' 'AGTA05000050.1' 'AGTA05000051.1' 
  'AGTA05000052.1' 'AGTA05000053.1' 'AGTA05000054.1' 'AGTA05000055.1' 'AGTA05000056.1' 
  'AGTA05000057.1' 'AGTA05000058.1' 'AGTA05000059.1' 'AGTA05000060.1' 'AGTA05000061.1' 
  'AGTA05000062.1' 'AGTA05000063.1' 'AGTA05000064.1' 'AGTA05000065.1' 'AGTA05000066.1' 
  'AGTA05000067.1' 'AGTA05000068.1' 'AGTA05000069.1' 'AGTA05000070.1' 'AGTA05000071.1' 
  'AGTA05000072.1' 'AGTA05000073.1' 'AGTA05000074.1' 'AGTA05000075.1' 'AGTA05000076.1' 
  'AGTA05000077.1' 'AGTA05000078.1' 'AGTA05000079.1' 'AGTA05000080.1' 'AGTA05000081.1' 
  'AGTA05000082.1' 'AGTA05000083.1' 'AGTA05000085.1' 'AGTA05000084.1' 'AGTA05000086.1' 
  'AGTA05000087.1' 'AGTA05000088.1' 'AGTA05000089.1' 'AGTA05000090.1' 'AGTA05000091.1' 
  'AGTA05000092.1' 'AGTA05000093.1' 'AGTA05000094.1' 'AGTA05000095.1' 'AGTA05000096.1' 
  'AGTA05000097.1' 'AGTA05000098.1' 'AGTA05000099.1' 'AGTA05000100.1' 'AGTA05000101.1' 
  'AGTA05000102.1' 'AGTA05000103.1' 'AGTA05000104.1' 'AGTA05000105.1' 'AGTA05000106.1' 
  'AGTA05000107.1' 'AGTA05000108.1' 'AGTA05000109.1' 'AGTA05000111.1' 'AGTA05000110.1' 
  'AGTA05000112.1' 'AGTA05000113.1' 'AGTA05000114.1' 'AGTA05000115.1' 'AGTA05000116.1' 
  'AGTA05000117.1' 'AGTA05000118.1' 'AGTA05000119.1' 'AGTA05000120.1' 'AGTA05000121.1' 
  'AGTA05000122.1' 'AGTA05000123.1' 'AGTA05000124.1' 'AGTA05000125.1' 'AGTA05000126.1' 
  'AGTA05000127.1' 'AGTA05000128.1' 'AGTA05000129.1' 'AGTA05000130.1' 'AGTA05000131.1' 
  'AGTA05000133.1' 'AGTA05000132.1' 'AGTA05000134.1' 'AGTA05000135.1' 'AGTA05000137.1' 
  'AGTA05000136.1' 'AGTA05000138.1' 'AGTA05000139.1' 'AGTA05000140.1' 'AGTA05000141.1' 
  'AGTA05000142.1' 'AGTA05000143.1' 'AGTA05000144.1' 'AGTA05000145.1' 'AGTA05000146.1' 
  'AGTA05000147.1' 'AGTA05000148.1' 'AGTA05000149.1' 'AGTA05000150.1' 'AGTA05000151.1' 
  'AGTA05000152.1' 'AGTA05000153.1' 'AGTA05000154.1' 'AGTA05000155.1' 'AGTA05000156.1' 
  'AGTA05000157.1' 'AGTA05000158.1' 'AGTA05000159.1' 'AGTA05000160.1' 'AGTA05000161.1' 
  'AGTA05000162.1' 'AGTA05000163.1' 'AGTA05000164.1' 'AGTA05000165.1' 'AGTA05000166.1' 
  'AGTA05000167.1' 'AGTA05000168.1' 'AGTA05000169.1' 'AGTA05000170.1' 'AGTA05000171.1' 
  'AGTA05000172.1' 'AGTA05000173.1' 'AGTA05000174.1' 'AGTA05000175.1' 'AGTA05000176.1' 
  'AGTA05000177.1' 'AGTA05000178.1' 'AGTA05000179.1' 'AGTA05000180.1' 'AGTA05000181.1' 
  'AGTA05000182.1' 'AGTA05000183.1' 'AGTA05000184.1' 'AGTA05000185.1' 'AGTA05000186.1' 
  'AGTA05000187.1' 'AGTA05000188.1' 'AGTA05000189.1' 'AGTA05000190.1' 'AGTA05000191.1' 
  'AGTA05000192.1' 'AGTA05000195.1' 'AGTA05000193.1' 'AGTA05000194.1' 'AGTA05000196.1' 
  'AGTA05000197.1' 'AGTA05000198.1' 'AGTA05000199.1' 'AGTA05000200.1' 'AGTA05000201.1' 
  'AGTA05000202.1' 'AGTA05000203.1' 'AGTA05000204.1' 'AGTA05000205.1' 'AGTA05000206.1' 
  'AGTA05000207.1' 'AGTA05000208.1' 'AGTA05000209.1' 'AGTA05000210.1' 'AGTA05000211.1' 
  'AGTA05000212.1' 'AGTA05000213.1' 'AGTA05000214.1' 'AGTA05000215.1' 'AGTA05000216.1' 
  'AGTA05000217.1' 'AGTA05000218.1' 'AGTA05000219.1' 'AGTA05000220.1' 'AGTA05000221.1' 
  'AGTA05000222.1' 'AGTA05000223.1' 'AGTA05000224.1' 'AGTA05000225.1' 'AGTA05000226.1' 
  'AGTA05000227.1' 'AGTA05000228.1' 'AGTA05000229.1' 'AGTA05000230.1' 'AGTA05000231.1' 
  'AGTA05000232.1' 'AGTA05000233.1' 'AGTA05000234.1' 'AGTA05000236.1' 'AGTA05000237.1' 
  'AGTA05000239.1' 'AGTA05000238.1' 'AGTA05000235.1' 'AGTA05000240.1' 'AGTA05000241.1' 
  'AGTA05000242.1' 'AGTA05000243.1' 'AGTA05000244.1' 'AGTA05000245.1' 'AGTA05000246.1' 
  'AGTA05000247.1' 'AGTA05000248.1' 'AGTA05000249.1' 'AGTA05000250.1'
)

# mkdir -p $OUTPUT_DIR

for CHR in '${CHROMOSOMES[@]}'; do
  echo 'Processing $CHR...'
  beagle -Xmx16g gt=$INPUT_VCF out=$OUTPUT_DIR/vcf_phased_$CHR nthreads=16 chrom=$CHR
  bcftools index -t $OUTPUT_DIR/vcf_phased_$CHR.vcf.gz
done

VCF_FILES=()
for CHR in '${CHROMOSOMES[@]}'; do
  VCF_FILES+=('$OUTPUT_DIR/vcf_phased_${CHR}.vcf.gz')
done

# echo 'Concatenating files...'
bcftools concat '${VCF_FILES[@]}' -Oz -o $FINAL_OUTPUT

# echo 'Done! Final file: $FINAL_OUTPUT'

# %% [markdown]
# ## extract vcf

# %%
import subprocess
import os

def vcftools_keep_sample(vcf_file, output_prefix, sample_file):
        
    cmd = f'vcftools --gzvcf {vcf_file} --maf 0.0001 --keep {sample_file} --min-alleles 2 --recode --recode-INFO-all --out {output_prefix}'
    print(f'Executing command: {cmd}')
    
    subprocess.run(cmd, shell=True, check=True)

    cmd_zip = f'bgzip -@ 16 {output_prefix}.recode.vcf'
    subprocess.check_output(cmd_zip, shell=True)

vcf_file = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05.vcf.recode.vcf.gz'

output_prefix_Hchi = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05_Hchi'
sample_Hchi_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/vcf_filter/species_pair/Hchi.txt'
vcftools_keep_sample(vcf_file, output_prefix_Hchi, sample_Hchi_path)

output_prefix_Hsau = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05_Hsau'
sample_Hsau_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/vcf_filter/species_pair/Hsau.txt'
vcftools_keep_sample(vcf_file, output_prefix_Hsau, sample_Hsau_path)

output_prefix_Hmic = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05_Hmic'
sample_Hmic_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/vcf_filter/species_pair/Hmic.txt'
vcftools_keep_sample(vcf_file, output_prefix_Hmic, sample_Hmic_path)

output_prefix_Lruf = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05_Lruf'
sample_Lruf_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/vcf_filter/species_pair/Lruf.txt'
vcftools_keep_sample(vcf_file, output_prefix_Lruf, sample_Lruf_path)

output_prefix_Ppun = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05_Ppun'
sample_Ppun_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/vcf_filter/species_pair/Ppun.txt'
vcftools_keep_sample(vcf_file, output_prefix_Ppun, sample_Ppun_path)

output_prefix_Pnye = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05_Pnye'
sample_Pnye_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/vcf_filter/species_pair/Pnye.txt'
vcftools_keep_sample(vcf_file, output_prefix_Pnye, sample_Pnye_path)
