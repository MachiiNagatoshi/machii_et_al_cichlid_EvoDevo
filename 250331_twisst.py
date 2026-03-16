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
# # 250331 twisst

# %% [markdown]
# ## make geno file

# %% [markdown]
# 以下のコードで、.genoファイルを作成した。
# > python3 /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/twisst/genomics_general/VCF_processing/parseVCF.py -i /Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05_twisst.recode.vcf.gz --skipIndels -o /Volumes/Transcend/250514_vcf_filtered/vcf_merged/vcf_phased_dataset1_missing0.9_maf0.05_twisst.geno.gz

# %% [markdown]
# | #CHROM | POS  | HC_2020_1 | HC_2020_3 | HC_2020_4 | SAMD00269291 | SAMD00269292 | SAMD00269293 | SAMD00269294 | SAMD00269295 | SAMD00269296 | SRR12394623 |
# |--------|------|-----------|-----------|-----------|--------------|--------------|--------------|--------------|--------------|--------------|-------------|
# | LG1    | 9130 | A\|A     | A\|A      | A\|A      | A\|A         | A\|A         | A\|A         | A\|A         | A\|A         | A\|A         | A\|A        |
# | LG1    | 9196 | A\|A     | A\|A      | A\|A      | A\|A         | A\|A         | A\|A         | A\|A         | A\|A         | A\|A         | A\|A        |
# | LG1    | 9227 | A\|A     | A\|A      | A\|A      | A\|A         | A\|A         | A\|A         | A\|A         | A\|A         | A\|A         | A\|A        |
# | LG1    | 9277 | A\|A     | A\|G      | A\|A      | A\|A         | A\|G         | A\|A         | A\|A         | A\|G         | A\|A         | G\|G        |
# | LG1    | 9288 | G\|G     | G\|G      | G\|G      | G\|G         | G\|G         | G\|G         | G\|G         | G\|G         | G\|G         | G\|G        |
# | LG1    | 9449 | G\|G     | G\|G      | G\|G      | G\|G         | G\|G         | G\|G         | G\|G         | G\|G         | G\|G         | G\|G        |
#

# %% [markdown]
# ## run raxml_sliding_window

# %%
# Importing libraries
import pandas as pd
import os
import subprocess
import shutil
import multiprocessing as mp

def mtrees_generator(inputs):
    # set variables
    windSize, name, path = inputs
    
    # Change the current working directory to the specified path
    os.chdir(path)

    # set input (.geno), output and raxml path
    geno_file='%s/%s.geno.gz'%(path,name)
    outname='%s_%sSNP'%(name,windSize)
    raxml_path = '/home/nakaharu/anaconda3/bin/raxmlHPC'

    # set command
    cmd="python /home/nakaharu/HDD4-2/machii/raxml_sliding_window_maf005/raxml_sliding_windows.py --windType sites --windSize %s --genoFile %s --prefix %s -M 2 --raxml %s"%(windSize,geno_file,outname,raxml_path)
    print(cmd)

    # Execute the command using subprocess
    p = subprocess.check_output(cmd, shell=True)


# Set the base path
base_path = '/home/nakaharu/HDD5/machii/vcf'

# Create a multiprocessing pool with 4 processes
p = mp.Pool(6)

# Create an empty list to store parameters
params = []

# Iterate over the specified window sizes
for windSize in [50]: #, 100

    # set input .geno file name
    name = 'vcf_phased_dataset1_missing0.9_maf0.05_twisst'

    # Append parameters to the list using f-string
    params.append([windSize, name, base_path])

# Print the list of parameters
print(params)  

# Execute mtrees_generator function with multiple parameters using multiprocessing
result = p.map(mtrees_generator, params)

# Close the multiprocessing pool
p.close()


# %%
$ conda activate raxml_sliding_window
$ python raxml_sliding_twisst.py

# %% [markdown]
# ## Run twisst

# %%
import pandas as pd
import os
import subprocess
import shutil
import multiprocessing as mp

## Running additional twisst groupings for 50, 100SNPs windSize

# generating parameters for twisst_run by adding list"params"
def param_gen(prefix, group_prefix, suffix_list, out):
    params = []
    for windSize in [50, 100]:
        base_path = "/Volumes/Transcend/250514_vcf_filtered/twisst/"
        trees_file = "%s/%s_%sSNP.trees.gz" % (base_path, prefix, windSize)
    
        for suffix in suffix_list:
            outname = "%s/%s_%sSNP_%s" % (base_path, prefix, windSize, suffix)
            group_tsv = "%s/%s_%s.txt" % (base_path, group_prefix, suffix)
            params.append([trees_file, outname, group_tsv, out])
    return params

def twisst_run(inputs):
    trees_file, outname, group_tsv, out = inputs
    
    df = pd.read_table(group_tsv, names=["indID", "group"])
    g_command = ""
    
    for group in list(set(df["group"])):
        g_command = g_command + "-g %s " % (group)
    
    cmd = "python /Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/twisst/twisst.py -t %s -w %s.weights.tsv.gz %s --method complete --groupsFile %s --outgroup %s" % (trees_file, outname, g_command, group_tsv, out)
    print(cmd)
    p = subprocess.check_output(cmd, shell=True)
    return p


if __name__ == "__main__":
    out="Nbri"
    params = param_gen(
        "vcf_phased_dataset1_missing0.9_maf0.05_twisst",
        "twisst",
        ["Victoria_Malawi", "Victoria_Tanganyika"],
        out)
    
    with mp.Pool(8) as p:
        result = p.map(twisst_run, params)


# %% [markdown]
# ## Visualize

# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
import os

def manhattan_plot_twisst(df_merged_data, chrom="scaffold", pos="mid", topo_columns=['%_topo1', '%_topo2', '%_topo3']):
    fig = plt.figure(figsize=(25, 6))
    ax = fig.add_subplot(111)
    
    df_merged_data[chrom] = df_merged_data[chrom].astype(str)
    
    last_xpos = 0
    xs_by_id = []
    x = []
    colors = ['#F2668B', '#9FC131', 'white']

    for seqid, group_data in df_merged_data.groupby(by=chrom, sort=False):
        for site in group_data[pos]:
            x.append(last_xpos + site)
        xs_by_id.append([seqid, last_xpos + (group_data[pos].iloc[0] + group_data[pos].iloc[-1]) / 2])
        last_xpos = x[-1]

    bottom = np.zeros(len(x))
    for i, topo in enumerate(topo_columns):
        values = df_merged_data[topo].values
        ax.bar(x, values, bottom=bottom, width=max(x)/len(x)/2, color=colors[i], 
               align='center', label=topo)
        bottom += values

    ax.set_xticks([v for c, v in xs_by_id])
    ax.set_xticklabels([c for c, v in xs_by_id], rotation=45, ha='right')
    
    ax.set_xlim(0, x[-1])
    ax.set_ylim(0, 100)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.set_ylabel('Topology Percentage')
    ax.set_xlabel('Genomic Position')
    ax.legend(title='Topologies')

    plt.tight_layout()
    return fig

def topology_violin_plot(df_merged_data, topo_columns=['%_topo1', '%_topo2', '%_topo3']):
    df_long = pd.melt(df_merged_data, value_vars=topo_columns, var_name='Topology', value_name='Percentage')
    fig, ax = plt.subplots(figsize=(7, 6))
    
    # 色を指定
    colors = {'%_topo1': '#4C9FD8', '%_topo2': '#727171', '%_topo3': '#D8D8D8'}
    
    # violinplotを描画し、色を指定
    sns.violinplot(x='Topology', y='Percentage', data=df_long, ax=ax,
                   palette=colors, order=topo_columns, cut=0)
    
    # %_topo1の上位0.5%の閾値を計算
    topo1_threshold = df_merged_data['%_topo1'].quantile(0.995)  # 上位0.5% = 99.5%パーセンタイル
    
    # 閾値を超えるデータを特定
    high_values_mask = df_merged_data['%_topo1'] >= topo1_threshold
    
    # 閾値を超えるデータを%_topo1のみでswarmplotで重ねて描画
    if high_values_mask.any():  # 該当データが存在する場合のみ
        # %_topo1の閾値を超えるデータのみを抽出
        df_high_topo1 = df_merged_data[high_values_mask][['%_topo1']]
        df_long_high = pd.melt(df_high_topo1, value_vars=['%_topo1'], var_name='Topology', value_name='Percentage')
        
        # 赤色のswarmplotを重ねて描画（%_topo1のみ）
        sns.swarmplot(x='Topology', y='Percentage', data=df_long_high, ax=ax,
                     color='red', size=2, alpha=0.8)
    
    # 閾値の位置に水平なグレーの点線を描画
    ax.axhline(y=topo1_threshold, color='gray', linestyle='--', alpha=0.7, linewidth=1)
    
    # y軸との交差点に閾値を表示
    ax.text(-0.05, topo1_threshold, f'{topo1_threshold:.1f}', 
            horizontalalignment='right', verticalalignment='center', 
            transform=ax.get_yaxis_transform(), 
            color='black', fontsize=10)
    
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Topology')
    
    # 統計検定
    stat, p_value = stats.mannwhitneyu(df_merged_data['%_topo1'], df_merged_data['%_topo2'])
    ax.text(0.5, 1.05, f'Mann-Whitney U test (topo1 vs topo2):\np-value = {p_value:.4e}',
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        
    plt.tight_layout()
    return fig

# ファイルのリスト
twisst_files = [
    'vcf_phased_dataset1_missing0.9_maf0.05_twisst_50SNP_Victoria_Malawi.weights.tsv',
    'vcf_phased_dataset1_missing0.9_maf0.05_twisst_50SNP_Victoria_Tanganyika.weights.tsv'
] #'vcf_phased_dataset1_missing1_maf0.05_twisst_100SNP_Victoria_Malawi.weights.tsv','vcf_phased_dataset1_missing1_maf0.05_twisst_100SNP_Victoria_Tanganyika.weights.tsv',


base_path = '/Volumes/Transcend/250514_vcf_filtered/twisst'

for twisst_file in twisst_files:
    # Twisst結果ファイルを読み込む
    df_twisst_results = pd.read_csv(os.path.join(base_path, twisst_file), sep='\t') #, comment='#'
    
    # ファイル名から対応する系統樹情報ファイルを決定
    snp_count = twisst_file.split('_')[6].split('SNP')[0]
    tree_info_file = f'vcf_phased_dataset1_missing0.9_maf0.05_twisst_{snp_count}SNP.data.tsv'
    
    # 系統樹の情報ファイルを読み込む
    df_tree_info = pd.read_csv(os.path.join(base_path, tree_info_file), sep='\t')
    
    # データフレームをマージ
    df_merged_data = pd.merge(df_tree_info, df_twisst_results, left_index=True, right_index=True)

    print(df_merged_data)
    
    # 各樹形の割合を計算
    total_topos = df_merged_data['topo1'] + df_merged_data['topo2'] + df_merged_data['topo3']
    df_merged_data['%_topo1'] = df_merged_data['topo1'] / total_topos * 100
    df_merged_data['%_topo2'] = df_merged_data['topo2'] / total_topos * 100
    df_merged_data['%_topo3'] = df_merged_data['topo3'] / total_topos * 100

    # 必要に応じて結果を保存
    # df_merged_data.to_csv(f'/Users/machiinagatoshi/Desktop/{twisst_file}_percentages.csv', index=False)
    
    # マンハッタンプロットを描画
    # fig_manhattan = manhattan_plot_twisst(df_merged_data)
    # fig_manhattan.savefig(os.path.join(base_path, f'manhattan_plot_{twisst_file.split(".weights")[0]}.png'), dpi=300, bbox_inches='tight')
    # plt.close(fig_manhattan)
    
    # バイオリンプロットを描画
    fig_violin = topology_violin_plot(df_merged_data)
    plt.rcParams["svg.fonttype"] = "none"
    fig_violin.savefig(os.path.join(base_path, f'violin_plot_{twisst_file.split(".weights")[0]}.svg'), bbox_inches='tight')
    plt.close(fig_violin)
    
    # 統計情報を表示
    print(f"\nResults for {twisst_file}:")
    print("Descriptive statistics:")
    print(df_merged_data[['%_topo1', '%_topo2', '%_topo3']].describe())
    
    print("\nMann-Whitney U test results (topo1 vs topo2):")
    stat, p_value = stats.mannwhitneyu(df_merged_data['%_topo1'], df_merged_data['%_topo2'])
    print(f"Statistic: {stat}")
    print(f"p-value: {p_value}")
    print("\n" + "="*50 + "\n")

print("Analysis complete. Check the output directory for plot images.")


# %%
import gffutils

# GFFファイルからデータベースを作成
gff_file = '/Volumes/research_backup/genome/Genome4VariantCall/M_zebra/Maylandia_zebra.M_zebra_UMD2a.109.gff3'
db_file = '/Volumes/research_backup/genome/Genome4VariantCall/M_zebra/Maylandia_zebra.M_zebra_UMD2a.109.gff3.db'  # データベースを保存するパスを指定

print(f"GFFファイルからデータベースを作成しています: {gff_file}")

# データベースがすでに存在する場合は読み込み、存在しない場合は新規作成
if not os.path.exists(db_file):
    print(f"データベースを作成して保存します: {db_file}")
    db = gffutils.create_db(gff_file, dbfn=db_file, force=True, keep_order=True,
                          merge_strategy='merge', sort_attribute_values=True)
else:
    print(f"既存のデータベースを読み込みます: {db_file}")
    db = gffutils.FeatureDB(db_file)

# %%
import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def process_twisst_file(weights_file, data_file, comparison_tag):
    # Twisst結果ファイルを読み込む（数値型を明示的に指定）
    df_weights = pd.read_csv(weights_file, sep='\t')
    # print(df_weights.shape)
    
    # 位置情報を含むファイルを読み込む
    df_data = pd.read_csv(data_file, sep='\t')
    # print(df_data.shape)
    
    # データフレームをマージ
    df_merged = pd.concat([df_data, df_weights], axis=1)
    
    # 各樹形の割合を計算
    total_topos = df_merged['topo1'] + df_merged['topo2'] + df_merged['topo3']
    df_merged['%_topo1'] = df_merged['topo1'] / total_topos * 100
    df_merged['%_topo2'] = df_merged['topo2'] / total_topos * 100
    df_merged['%_topo3'] = df_merged['topo3'] / total_topos * 100
    
    df_merged['Comparison_Tag'] = comparison_tag
    
    return df_merged


def plot_twisst_region(df, region, base_path, db, mergin):
    
    # フォント設定
    plt.rcParams['svg.fonttype'] = 'none'
    
    # 領域情報の解析
    chrom, start, end = region.split('_')
    start = max(0, int(start)-int(mergin))
    end = int(end) + int(mergin)
    
    # 指定された領域のデータを抽出
    df_region = df[(df['scaffold'] == chrom) & 
                   (df['start'] >= start) & 
                   (df['end'] <= end)]

    df_region = df_region.sort_values(['start'], ascending=True)
    
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(15, 8), sharex=True, gridspec_kw={'height_ratios': [0.5, 4, 4]})
    
    colors_mal = ['#4c9fd8', '#717071', '#D9D9D9']
    colors_tan = ['#e36485', '#717071', '#D9D9D9']

    # 遺伝子情報を描画
    for gene in db.region(seqid=chrom, start=start, end=end, featuretype='gene'):
        gene_start = max(gene.start, start)
        gene_end = min(gene.end, end)
        rect = Rectangle((gene_start, 0), gene_end - gene_start, 1, 
                        fill=True, facecolor='lightgrey', edgecolor='none')
        ax0.add_patch(rect)
        
        # エキソンを濃いグレーで表示
        for exon in db.children(gene, featuretype='exon'):
            exon_start = max(exon.start, start)
            exon_end = min(exon.end, end)
            if exon_start < end and exon_end > start:
                rect = Rectangle((exon_start, 0), exon_end - exon_start, 1,
                               fill=True, facecolor='darkgrey', edgecolor='none')
                ax0.add_patch(rect)
        
        gene_name = gene.attributes.get('Name', [gene.id.split(':')[1]])[0]
        ax0.text((gene_start + gene_end) / 2, 0.5, gene_name, 
                 ha='center', va='center', fontsize=8)
    
    ax0.set_ylim(0, 1)
    ax0.set_yticks([])
    ax0.set_xlim(start, end)
    ax0.axis('off')
    
    for ax, comparison, colors in zip([ax1, ax2], ['Victoria_Malawi', 'Victoria_Tanganyika'], [colors_mal, colors_tan]):
        df_comp = df_region[df_region['Comparison_Tag'] == comparison]

        # print(df_comp)
        
        ax.bar(df_comp['start'], df_comp['%_topo1'], width=df_comp['end'] - df_comp['start'],
               color=colors[0], label='Topo1', align='edge')
        ax.bar(df_comp['start'], df_comp['%_topo2'], width=df_comp['end'] - df_comp['start'],
               bottom=df_comp['%_topo1'], color=colors[1], label='Topo2', align='edge')
        ax.bar(df_comp['start'], df_comp['%_topo3'], width=df_comp['end'] - df_comp['start'],
               bottom=df_comp['%_topo1'] + df_comp['%_topo2'], color=colors[2], label='Topo3', align='edge')
        
        ax.set_title(f'{comparison}')
    
    ax2.set_xlim(start, end)
    ax2.set_xticks(np.linspace(start, end, 6))
    ax2.set_xticklabels([f'{int(x):,}' for x in np.linspace(start, end, 6)])
    
    plt.suptitle(f'Twisst Results for {region}')
    plt.tight_layout()
    
    output_file = f'{base_path}/twisst_plot_{region}.svg'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved as: {output_file}")


# ファイルパスを設定
base_path = '/Volumes/Transcend/250514_vcf_filtered/twisst/'
tree_info_file = os.path.join(base_path, 'vcf_phased_dataset1_missing0.9_maf0.05_twisst_50SNP.data.tsv')
vic_mal_file = os.path.join(base_path, 'vcf_phased_dataset1_missing0.9_maf0.05_twisst_50SNP_Victoria_Malawi.weights.tsv',)
vic_tan_file = os.path.join(base_path, 'vcf_phased_dataset1_missing0.9_maf0.05_twisst_50SNP_Victoria_Tanganyika.weights.tsv')


# データの処理
vic_mal = process_twisst_file(vic_mal_file, tree_info_file, 'Victoria_Malawi')
vic_tan = process_twisst_file(vic_tan_file, tree_info_file, 'Victoria_Tanganyika')

# データの結合
df_merged = pd.concat([vic_mal, vic_tan], ignore_index=True)

regions = ['LG8_16840034_16878542']

for region in regions:
    plot_twisst_region(df_merged, region, base_path, db,100000)
