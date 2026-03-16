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
# # 250401 Calculate LD

# %%
import gzip
import os

def filter_monomorphic_sites(input_vcf, output_vcf):

    with gzip.open(input_vcf, 'rt') as in_f, open(output_vcf, 'w') as out_f:
        for line in in_f:
            if line.startswith('#'):
                out_f.write(line)
                continue
            
            fields = line.strip().split('\t')
            genotypes = fields[9:]
            
            alleles = []
            for gt in genotypes:
                alleles.extend(gt.split(':')[0].replace('|', '/').split('/'))
            
            unique_alleles = set(alleles) - {'.'}
            
            if len(unique_alleles) > 1:
                out_f.write(line)

vcf_files = [
    'vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Hmic.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Hsau.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Lruf.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Pnye.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Ppun.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hmic_Hsau.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hmic_Lruf.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hmic_Pnye.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hmic_Ppun.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hsau_Lruf.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hsau_Pnye.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Hsau_Ppun.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Lruf_Pnye.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Lruf_Ppun.recode.vcf.gz',
    'vcf_phased_dataset1_missing0.9_maf0.05_Pnye_Ppun.recode.vcf.gz'
]

base_path = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged'

for vcf_file in vcf_files:
    output_file = vcf_file.replace('.gz', '').replace('.vcf', '_removeNonSNPs.vcf')
    print(f'処理中: {vcf_file}')
    filter_monomorphic_sites(f'{base_path}/{vcf_file}', output_file)
    print(f'完了: {output_file}')


# %%
# bash
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Hmic.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Hsau.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Lruf.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Pnye.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hchi_Ppun.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hmic_Hsau.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hmic_Lruf.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hmic_Pnye.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hmic_Ppun.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hsau_Lruf.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hsau_Pnye.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Hsau_Ppun.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Lruf_Pnye.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Lruf_Ppun.recode_removeNonSNPs.vcf -@ 128
bgzip /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/vcf_phased_dataset1_missing0.9_maf0.05_Pnye_Ppun.recode_removeNonSNPs.vcf -@ 128

# %% [markdown]
# ## Preprocessing

# %%
import pandas as pd

def merge_adjacent_regions(df):
    merged_regions = []
    
    for chrom in df['CHROM'].unique():
        chrom_data = df[df['CHROM'] == chrom].sort_values('POS_START')
        
        current_start = chrom_data.iloc[0]['POS_START']
        current_end = chrom_data.iloc[0]['POS_END']
        
        for i in range(1, len(chrom_data)):
            next_start = chrom_data.iloc[i]['POS_START']
            next_end = chrom_data.iloc[i]['POS_END']
            
            if next_start <= current_end:
                current_end = max(current_end, next_end)
            else:
                merged_regions.append({
                    'CHROM': chrom,
                    'POS_START': current_start,
                    'POS_END': current_end
                })
                current_start = next_start
                current_end = next_end
        
        merged_regions.append({
            'CHROM': chrom,
            'POS_START': current_start,
            'POS_END': current_end
        })
    return pd.DataFrame(merged_regions)

hdr_top1_dup4_dataset1_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/result/hdr_top1%_dup4_dataset1.csv'
df_hdr_top1_dup4_dataset1 = pd.read_csv(hdr_top1_dup4_dataset1_path)

df_hdr_top1_dup4_dataset1 = df_hdr_top1_dup4_dataset1[['CHROM', 'POS_START', 'POS_END']]
df_hdr_top1_dup4_dataset1['POS_START'] = df_hdr_top1_dup4_dataset1['POS_START']-100000
df_hdr_top1_dup4_dataset1['POS_END'] = df_hdr_top1_dup4_dataset1['POS_END']+100000

merged_df = merge_adjacent_regions(df_hdr_top1_dup4_dataset1)

print('結合後の領域:')
print(merged_df)


# %% [markdown]
# Filter short chromosome

# %%

fai_path = '/Volumes/research_backup/genome/Genome4VariantCall/M_zebra/Maylandia_zebra.M_zebra_UMD2a.dna_sm.toplevel.masked.fa.fai'
df_fai= pd.read_csv(fai_path, sep='\t', header=None)

# merge df_fai
merged_df = pd.merge(merged_df, df_fai, how='left', left_on='CHROM', right_on=0)

merged_df = merged_df[(merged_df['POS_START']>0) & (merged_df['POS_END']<merged_df[1])]
merged_df = merged_df[['CHROM', 'POS_START', 'POS_END']]

merged_df.to_csv('/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/LDcal_HighFst_vcftools_with100000_merged_regions.bed', sep='\t', index=False)

# %% [markdown]
# ## Calculate LD

# %%
import subprocess
import itertools
from concurrent.futures import ThreadPoolExecutor
import os

# set path
base_path = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/'
bed_file = f'{base_path}/LDcal_HighFst_vcftools_with100000_merged_regions.bed'
vcf_path = '/Volumes/Transcend/250514_vcf_filtered/vcf_merged'

N_THREADS = 16

species = ['Hchi', 'Hmic', 'Hsau', 'Lruf', 'Pnye', 'Ppun']

def run_vcftools(combination):
    sp1, sp2 = combination
    vcf_file = f'{vcf_path}/vcf_phased_dataset1_missing0.9_maf0.05_{sp1}_{sp2}.recode_removeNonSNPs.vcf.gz'
    output = f'LD_{sp1}_{sp2}_region'
    
    cmd = [
        'vcftools',
        '--gzvcf', vcf_file,
        '--bed', bed_file,
        '--hap-r2',
        '--ld-window-bp', '150000',
        '--out', output
    ]
    print(cmd)
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(f'Completed: {sp1}-{sp2}')
        return result
    except Exception as e:
        print(f'Error processing {sp1}-{sp2}: {str(e)}')
        return None

combinations = list(itertools.combinations(species, 2))

print(f'Starting analysis with {N_THREADS} threads')

with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
    results = list(executor.map(run_vcftools, combinations))

for combination, result in zip(combinations, results):
    if result is None:
        print(f'Failed: {combination[0]}-{combination[1]}')
    elif result.returncode != 0:
        print(f'Error in {combination[0]}-{combination[1]}: {result.stderr}')


# %%
# bash

for file in /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/*.hap.ld; do
    bgzip -f "$file" -@ 128
done

# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import gc

# リージョン情報の読み込み
path_region = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250331_dataset1_realanysis/dataset1/hdr_top1%_dup4_dataset1.csv'
df_region = pd.read_csv(path_region)

# 解析する種ペアのリスト
species = ['Hchi_Hmic', 'Hchi_Hsau', 'Hchi_Pnye', 'Hchi_Ppun','Hchi_Lruf',
           'Hmic_Hsau', 'Hmic_Pnye', 'Hmic_Ppun', 'Hmic_Lruf', 'Hsau_Lruf', 
           'Hsau_Pnye', 'Hsau_Ppun', 'Lruf_Pnye', 'Lruf_Ppun', 'Pnye_Ppun']

# マージンのリスト
list_margin = [0, 1000, 5000, 10000, 25000, 50000, 100000]

# 各種ペアについて処理
for sp in species:
    print(f"\nProcessing species pair: {sp}")

    # 各マージン値について個別のファイルを作成
    result_files = {margin: open(f'LD_{sp}_margin{margin}_vcftools.csv', 'w') for margin in list_margin}
    
    # ヘッダーを書き込み
    for f in result_files.values():
        f.write('species_pair,CHROM,HDR_START,HDR_END,LDcal_START,LDcal_END,LD_mean\n')


    # LDファイルの読み込み（1回のみ）
    file_path = f"/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/LD_{sp}_region.hap.ld.gz"

    try:
        print(f"Reading file: {file_path}")
        df_tmp = pd.read_csv(file_path, sep="\t")
        df_tmp = df_tmp.dropna(subset=['R^2'])
        
        # 各リージョンについて処理
        for tup in df_region.itertuples():
            chr = tup.CHROM
            base_start = int(tup.POS_START)
            base_end = int(tup.POS_END)
            
            # 各マージン値について処理し、直接ファイルに書き込み
            for margin in list_margin:
                start = base_start - margin
                end = base_end + margin
                
                region_df_tmp = df_tmp[
                    (df_tmp['CHR'] == chr) & 
                    (df_tmp['POS1'] >= start) & 
                    (df_tmp['POS1'] <= end) & 
                    (df_tmp['POS2'] >= start) & 
                    (df_tmp['POS2'] <= end)
                ]
                
                if not region_df_tmp.empty:
                    ld_mean = region_df_tmp['R^2'].mean()
                    result_files[margin].write(f'{sp},{chr},{base_start},{base_end},{start},{end},{ld_mean}\n')
                else:
                    result_files[margin].write(f'{sp},{chr},{base_start},{base_end},{start},{end},nan\n')

        del df_tmp
        gc.collect()
        
    except Exception as e:
        print(f"Error processing {sp}: {str(e)}")
        continue

# ファイルを閉じる
for f in result_files.values():
    f.close()


# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# 保存先ディレクトリの作成（存在しない場合）
output_dir = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/figures'
os.makedirs(output_dir, exist_ok=True)

# merginのリスト
mergin_list = [0, 1000, 5000, 10000, 25000, 50000, 100000]

# 種のリスト
species_list = ['Hchi_Hmic', 'Hchi_Hsau', 'Hchi_Pnye', 'Hchi_Ppun', 'Hchi_Lruf', 'Hmic_Hsau', 'Hmic_Pnye', 'Hmic_Ppun', 'Hmic_Lruf', 'Hsau_Lruf', 'Hsau_Pnye', 'Hsau_Ppun', 'Lruf_Pnye', 'Lruf_Ppun', 'Pnye_Ppun']

list_df_LD_mean = []

for mergin in mergin_list:
    for sp in species_list:
        # データの読み込み
        file_path = f'/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/LD_mean_output/LD_{sp}_margin{mergin}_vcftools.csv'
        df_tmp = pd.read_csv(file_path)
    
        df_tmp['mergin'] = mergin
        list_df_LD_mean.append(df_tmp)

df_LD_mean_concat = pd.concat(list_df_LD_mean, axis=0)

# 新しい列 'Lip_tag' を作成し、デフォルト値を 'Normal_pair' に設定
df_LD_mean_concat['Lip_tag'] = 'Normal_pair'

# Species列に 'Hchi' が含まれている行を 'Lip_pair' に設定
df_LD_mean_concat.loc[df_LD_mean_concat['species_pair'].str.contains('Hchi'), 'Lip_tag'] = 'Lip_pair'

df_LD_mean_concat = df_LD_mean_concat.dropna()

df_LD_mean_concat

color_list = ['#F2668B', '#5E8B84']

fig = plt.figure(figsize=(7, 6))
sns.boxenplot(data=df_LD_mean_concat, x='mergin', y='LD_mean', hue='Lip_tag', gap=.1, width_method='exponential', palette=color_list, flier_kws=dict(linewidth=.2))

plt.rcParams['svg.fonttype'] = 'none'
plt.savefig('/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/figures/LD_cal_mean_with_mergin.svg', dpi=300, bbox_inches='tight')

plt.show()

# %%
df_LD_mean_concat.to_csv('/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/figures/df_LD_mean_concat.csv')

# %%
from scipy.stats import mannwhitneyu

mergin_list = [0, 1000, 5000, 10000, 25000, 50000, 100000]

for mergin in mergin_list:

    print(f'test for {mergin}')
    df_LD_mean_concat_mergin = df_LD_mean_concat[df_LD_mean_concat['mergin']==mergin]
    
    # 各グループのデータを取得
    group1 = df_LD_mean_concat_mergin[df_LD_mean_concat_mergin['Lip_tag'] == 'Lip_pair']['LD_mean']
    group2 = df_LD_mean_concat_mergin[df_LD_mean_concat_mergin['Lip_tag'] == 'Normal_pair']['LD_mean']
    
    # Mann-Whitney U検定を実行
    statistic, p = mannwhitneyu(group1, group2, method='asymptotic', alternative='two-sided')
    
    print(f'bonfferoni P: {p*len(mergin_list)}')

# %%
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

# merginのリスト
mergin_list = [0, 1000, 5000, 10000, 25000, 50000, 100000]

# 空のデータフレームを作成
df_Fisher_result = pd.DataFrame(columns=['mergin', 'top500_contingency', 'top500_odds_ratio', 
                                       'top500_p_value', 'top500_adj_p_value', 
                                       'propotion_Normal', 'propotion_Lip'])

# 各merginについて処理を実行
for mergin in mergin_list:
    print(f'test for {mergin}')
    df_LD_mean_concat_mergin = df_LD_mean_concat[df_LD_mean_concat['mergin']==mergin]
    
    # 全体の上位500を抽出
    df_top500 = df_LD_mean_concat_mergin.nlargest(500, 'LD_mean')    
    
    # 全体の上位500におけるNormal_pairとLip_pairのカウント
    top500_normal_count = sum(df_top500['Lip_tag'] == 'Normal_pair')
    top500_lip_count = sum(df_top500['Lip_tag'] == 'Lip_pair')
    
    # データセット全体でのNormal_pairとLip_pairの総数
    total_normal_count = sum(df_LD_mean_concat_mergin['Lip_tag'] == 'Normal_pair')
    total_lip_count = sum(df_LD_mean_concat_mergin['Lip_tag'] == 'Lip_pair')
    
    # 全体の比率の分析のための2x2分割表
    top500_contingency = np.array([
        [top500_lip_count, total_lip_count],
        [top500_normal_count, total_normal_count]
    ])
    
    # Fisherの正確確率検定を実行
    top500_odds_ratio, top500_p_value = fisher_exact(top500_contingency)
    
    # 新しい行のデータを辞書として作成
    new_row = {
        'mergin': mergin,
        'top500_contingency': [top500_contingency],
        'top500_odds_ratio': top500_odds_ratio,
        'top500_p_value': top500_p_value,
        'top500_adj_p_value': min(top500_p_value * len(mergin_list), 1),
        'propotion_Normal': (top500_normal_count/500)*100,
        'propotion_Lip': (top500_lip_count/500)*100
    }
    
    # データフレームに新しい行を追加
    df_Fisher_result = pd.concat([df_Fisher_result, pd.DataFrame([new_row])], ignore_index=True)
    
    print('Overall Top 500 Analysis:')
    print(f'Normal pairs in top 500: {top500_normal_count}')
    print(f'Lip pairs in top 500: {top500_lip_count}')
    print(f'Total Normal pairs: {total_normal_count}')
    print(f'Total Lip pairs: {total_lip_count}')
    print(f'Odds ratio: {top500_odds_ratio:.3f}')
    print(f'P-value: {top500_p_value:.2e}')
    print('\n')

# 全行のデータを辞書として作成
top500_normal_count = sum(df_LD_mean_concat_mergin['Lip_tag'] == 'Normal_pair')
top500_lip_count = sum(df_LD_mean_concat_mergin['Lip_tag'] == 'Lip_pair')

total_row = {
    'mergin': 'total',
    'top500_contingency': [total_lip_count, total_normal_count],
    'propotion_Normal': (top500_normal_count/(total_lip_count+total_normal_count))*100,
    'propotion_Lip': (top500_lip_count/(total_lip_count+total_normal_count))*100
}

df_Fisher_result = pd.concat([df_Fisher_result, pd.DataFrame([total_row])], ignore_index=True)

# 結果の確認
print('\nFinal DataFrame:')
print(df_Fisher_result)

# データを整形（ロングフォーマットに変換）
df_long = pd.melt(
    df_Fisher_result,
    id_vars=['mergin'],
    value_vars=['propotion_Normal', 'propotion_Lip'],
    var_name='Type',
    value_name='Proportion'
)

# Typeの名前を整形（'propotion_'を削除）
df_long['Type'] = df_long['Type'].str.replace('propotion_', '')

# プロットの設定
plt.figure(figsize=(6, 4))
sns.barplot(
    data=df_long,
    x='mergin',
    y='Proportion',
    hue='Type',
    hue_order=['Lip', 'Normal'],  # 順番を指定
    palette=['#F2668B', '#5E8B84']  # カラーパレットを設定
)

# グラフの装飾
plt.title('Proportion of Normal and Lip Pairs by Margin', fontsize=8, pad=15)
plt.xlabel('Margin', fontsize=8)
plt.ylabel('Proportion (%)', fontsize=8)
plt.xticks(rotation=45)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.legend(title='Type')

# レイアウトの調整
plt.tight_layout()

plt.rcParams['svg.fonttype'] = 'none'
plt.savefig('/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/figures/overall_odds.svg', dpi=300, bbox_inches='tight')

plt.show()

# %%
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns

mergin = 25000

# merginごとに区切る
df_mergin = df_LD_mean_concat[df_LD_mean_concat['mergin']==mergin]

# 染色体ごとの解析
normal_top500 = df_mergin[df_mergin['Lip_tag'] == 'Normal_pair'].nlargest(500, 'LD_mean')
lip_top500 = df_mergin[df_mergin['Lip_tag'] == 'Lip_pair'].nlargest(500, 'LD_mean')

def get_chrom_counts(data):
    return data['CHROM'].value_counts().to_dict()

normal_counts = get_chrom_counts(normal_top500)
lip_counts = get_chrom_counts(lip_top500)

all_chroms = sorted(list(set(list(normal_counts.keys()) + list(lip_counts.keys()))))

results = []
for chrom in all_chroms:
    normal_count = normal_counts.get(chrom, 0)
    lip_count = lip_counts.get(chrom, 0)
    
    contingency_table = np.array([
        [lip_count, total_lip_count],
        [normal_count, total_normal_count]
    ])
    
    odds_ratio, p_value = fisher_exact(contingency_table)
    
    results.append({
        'Chromosome': chrom,
        'Normal_count': normal_count,
        'Lip_count': lip_count,
        'Odds_ratio': odds_ratio,
        'P_value': p_value
    })

results_df = pd.DataFrame(results)

# Bonferroni補正
n_tests = len(all_chroms)
results_df['P_value_adjusted'] = results_df['P_value'] * n_tests
results_df['P_value_adjusted'] = results_df['P_value_adjusted'].clip(upper=1)
results_df = results_df.sort_values('Odds_ratio', ascending=False)


# 結果のDataFrameからinfを除外
results_df = results_df[results_df['Odds_ratio'] != float('inf')]

# 散布図をプロット
fig, ax1 = plt.subplots(1, 1, figsize=(6, 5))

# x軸の目盛りを設定
x_positions = np.arange(len(results_df['Chromosome']))
ax1.set_xticks(x_positions)
ax1.set_xticklabels(results_df['Chromosome'], rotation=45, ha='right')

# 点のサイズを調整（合計カウントに基づく）
sizes = (results_df['Normal_count'] + results_df['Lip_count']) * 2

# P値に基づいて色を設定
colors = ['#FF0000' if p < 0.01 else  # 赤
          '#FFB2B2' if p < 0.05 else  # 灰色
          '#808080' for p in results_df['P_value_adjusted']]

# 散布図をプロット
ax1.scatter(x_positions, results_df['Odds_ratio'], 
           c=colors, s=sizes, alpha=0.6)

ax1.set_title('Chromosome Enrichment Analysis')
ax1.set_ylabel('Odds Ratio')

# 凡例を作成（色とサイズの両方）
# 色の凡例
legend_elements_color = [
    plt.scatter([], [], c='#FF0000', label='p < 0.01', s=100),
    plt.scatter([], [], c='#FFB2B2', label='p < 0.05', s=100),
    plt.scatter([], [], c='#808080', label='p ≥ 0.05', s=100)
]

# サイズの凡例
# 例として50, 100, 200のカウントに対応するサイズを表示
size_examples = [10, 50, 100, 200]
legend_elements_size = [
    plt.scatter([], [], c='gray', alpha=0.6, s=count*2,
               label=f'n = {count}') for count in size_examples
]

# 2つの凡例を別々に配置
# 色の凡例を左に
first_legend = ax1.legend(handles=legend_elements_color, 
                         title='Significance', 
                         loc='upper right')
# サイズの凡例を右に
ax1.legend(handles=legend_elements_size, 
          title='Total counts',
          loc='center right')
# 最初の凡例を追加し直す
ax1.add_artist(first_legend)

# より広いマージンを確保するために tight_layout の調整
plt.tight_layout(rect=[0, 0, 0.85, 1])

plt.rcParams['svg.fonttype'] = 'none'
plt.savefig(f'/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/figures/LD_plot_fisher_{mergin}.svg', dpi=300, bbox_inches='tight')
plt.show()

# 結果をCSVファイルに保存
results_df.to_csv(f'/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/figures/enrichment_analysis_results_{mergin}.csv', index=False)

# %% [markdown]
# ## make boxen plot

# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from IPython.display import display

output_dir = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/figures'
os.makedirs(output_dir, exist_ok=True)

species = ['Hchi_Hmic', 'Hchi_Hsau', 'Hchi_Pnye', 'Hchi_Ppun', 'Hchi_Lruf', 'Hmic_Hsau', 'Hmic_Pnye', 'Hmic_Ppun', 'Hmic_Lruf', 'Hsau_Lruf', 'Hsau_Pnye', 'Hsau_Ppun', 'Lruf_Pnye', 'Lruf_Ppun', 'Pnye_Ppun']

list_df_LD_tmp = []

for sp in species:
    file_path = f'/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/LD_{sp}_region.hap.ld.gz'
    df_LD_tmp = pd.read_csv(file_path, sep='\t', compression='gzip')
    df_LD_tmp['species_pair'] = sp
    list_df_LD_tmp.append(df_LD_tmp)
    
df_LD_concat = pd.concat(list_df_LD_tmp, axis=0)

display(df_LD_concat)

regions = ['LG8_16830001_16900001']

color_list_LD = ['#F2668B', '#03A688', '#9FC131', '#40E0D0', '#FFB347', '#b3b3b3', '#b3b3b3', '#b3b3b3', '#b3b3b3', '#b3b3b3', '#b3b3b3', '#b3b3b3', '#b3b3b3', '#b3b3b3', '#b3b3b3']

for region in regions:
    chrom, start, end = region.split('_')
    start, end = int(start), int(end)
    
    df_LD_region = df_LD_concat[(df_LD_concat['CHR'] == chrom) & (df_LD_concat['POS1'] >= start) & (df_LD_concat['POS1'] <= end) & (df_LD_concat['POS2'] >= start) & (df_LD_concat['POS1'] <= end)]
    
    fig, ax = plt.subplots(figsize = (7, 5))
        
    # plot violin_plot
    sns.boxenplot(ax=ax, x='species_pair', y='R^2', data=df_LD_region, palette=color_list_LD)
    
    # Graph parameters
    plt.tight_layout()
    plt.xticks(rotation=-90)
    
    output_file = f'{output_dir}/LD_violin_vlnplot_{region}.svg'
    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

# %% [markdown]
# ## make heatmap

# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

output_dir = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/heatmap_output"
os.makedirs(output_dir, exist_ok=True)

species = ['Hchi_Hmic', 'Hchi_Hsau', 'Hchi_Pnye', 'Hchi_Ppun', 'Hchi_Lruf', 'Hmic_Hsau', 'Hmic_Pnye', 'Hmic_Ppun', 'Hmic_Lruf', 'Hsau_Lruf', 'Hsau_Pnye', 'Hsau_Ppun', 'Lruf_Pnye', 'Lruf_Ppun', 'Pnye_Ppun']

info_path = "/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/LD_info.txt"
df_info = pd.read_csv(info_path, sep="\t")

regions = ['LG8_16830001_16900001']

for sp in species:
    file_path = f"/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/LD/LD_{sp}_region.hap.ld.gz"
    df = pd.read_csv(file_path, sep="\t", compression="gzip")
    
    print(f"Processing species: {sp}")

    # make heatmap in set regions
    for region in regions:
        print(region)
        chr, start, end = region.split('_')
        start, end = int(start), int(end)
        
        # filtering
        region_df = df[(df['CHR'] == chr) & (df['POS1'] >= start) & (df['POS1'] <= end) & 
                       (df['POS2'] >= start) & (df['POS2'] <= end)]
        
        if region_df.empty:
            print(f"リージョン {region} にはデータがありません。スキップします。")
            continue
        
        # make pivot df
        pivot_df = region_df.pivot(index="POS1", columns="POS2", values="R^2")
        
        # set mask
        mask = np.tril(np.ones_like(pivot_df, dtype=bool))
        
        # draw heatmap
        plt.figure(figsize=(12, 10))
        heatmap = sns.heatmap(pivot_df, cmap="YlOrRd", square=True, mask=mask)

        # graph parameters
        plt.xlabel("POS2")
        plt.ylabel("POS1")
        plt.title(f"Linkage Disequilibrium (R^2) Heatmap for {sp} - {region}")
        plt.tight_layout()
        
        filename = f"LD_heatmap_{region}_{sp}.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        
        plt.close()
        
        print(f"ヒートマップが保存されました: {filepath}")

print(f"すべてのヒートマップが {output_dir} ディレクトリに保存されました。")
