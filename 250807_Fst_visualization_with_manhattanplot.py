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
# # 250807 Visualization

# %% [markdown]
# ## Visualize upset plot

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from upsetplot import plot, from_memberships
import os

def upsetplot_fst_with_violin(data, chrom="CHROM", pos="POS", fst="WEIR_AND_COCKERHAM_FST", top_threshold=0.99, max_value_of_yaxis=1):
    fig = plt.figure(figsize=(25, 3))
    gs = fig.add_gridspec(1, 2, width_ratios=[22, 3])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    
    data[[chrom]] = data[[chrom]].astype(str)
    
    last_xpos = 0
    xs_by_id = []
    x, y, c = [], [], []
    colors = ["#C0C0C0", "#717375"]
    color_index = 0
    for seqid, group_data in data.groupby(by=chrom, sort=False):
        color = colors[color_index]
        for site, fst_value in zip(group_data[pos], group_data[fst]):
            x.append(last_xpos + site)
            y.append(fst_value)
            c.append(color)
        xs_by_id.append([seqid, last_xpos + (group_data[pos].iloc[0] + group_data[pos].iloc[-1]) / 2]) 
        last_xpos = x[-1]
        color_index = (color_index + 1) % 2
    
    fst_one_ratio = (data[fst] == 1).sum() / len(data)
    if fst_one_ratio >= (1 - top_threshold):
        c = np.where(np.array(y) == 1, "r", c)
    else:
        threshold = data[fst].quantile(top_threshold)
        c = np.where(np.array(y) >= threshold, "r", c)
    
    ax1.scatter(x, y, c=c, edgecolors="none", s=0.5)
    
    ax1.set_xticks([v for c, v in xs_by_id]) 
    ax1.set_xticklabels([c for c, v in xs_by_id])
    
    ax1.set_xlim(0, x[-1])
    ax1.set_ylim(ymin=-0.2, ymax=1)
    
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    sns.violinplot(y=data[fst], ax=ax2, color="gray", orient="v", split=False, inner=None, cut=0)
    
    if fst_one_ratio >= (1 - top_threshold):
        y_data = data[data[fst] == 1][fst]
    else:
        y_data = data[data[fst] >= threshold][fst]
    
    x_data = [0] * len(y_data)
    ax2.scatter(x_data, y_data, color="red", s=5)
    ax2.set_ylim(ax1.get_ylim())
    ax2.axis("off")
    # Fstの平均値を表示
    mean_fst = data[fst].mean()
    ax2.text(0.5, 0.5, f"Mean: {mean_fst:.3f}", ha='center', va='center', fontsize=20)
    # 上位0.5％のFst値を表示
    if fst_one_ratio >= (1 - top_threshold):
        ax2.text(0.5, 0.8, f"Fst=1: {fst_one_ratio*100:.1f}%", ha='center', va='center', fontsize=20)
    else:
        ax2.text(0.5, 0.8, f"Top{(1-top_threshold)*100:.1f}%: {threshold:.3f}", ha='center', va='center', fontsize=20)
    plt.tight_layout()
    
    return ax1, ax2


# %%
### pi raw file の加工
# データをDataFrameに読み込む
vcftools_result_dir = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/vcftools_result_window10000'

pi_files = [
    'pi_10000_Hchi.txt.windowed.pi.gz', 'pi_10000_Hsau.txt.windowed.pi.gz', 'pi_10000_Lruf.txt.windowed.pi.gz',
    'pi_10000_Hmic.txt.windowed.pi.gz', 'pi_10000_Pnye.txt.windowed.pi.gz', 'pi_10000_Ppun.txt.windowed.pi.gz'
]

df_list_pi = []
df_out_pi = pd.DataFrame()

for pi_file in pi_files:

    # pathとprefixの定義
    prefix = os.path.splitext(pi_file)[0].split('.')[0].replace('_10000','')
    pi_path = os.path.join(vcftools_result_dir, pi_file)
    
    df_tmp_pi = pd.read_csv(pi_path, sep='\t')

    # CHROM_POS列の作成
    df_tmp_pi['CHROM_POS'] = df_tmp_pi['CHROM'] + '_' + df_tmp_pi['BIN_START'].astype(str) + '_' + df_tmp_pi['BIN_END'].astype(str)
    df_tmp_pi = df_tmp_pi.set_index('CHROM_POS')

    # カラムを整形
    df_tmp_pi = df_tmp_pi.rename(columns={'PI': prefix})
    df_tmp_pi = df_tmp_pi[prefix]

    df_list_pi.append(df_tmp_pi)


# concat関数を使用して横軸に連結
df_out_pi = pd.concat(df_list_pi, axis=1)

# 重複する列名を処理
df_out_pi = df_out_pi.loc[:,~df_out_pi.columns.duplicated()]
df_out_pi = df_out_pi[['pi_Hchi', 'pi_Hsau', 'pi_Lruf', 'pi_Hmic', 'pi_Ppun', 'pi_Pnye']]

# CHROM, POS列の作成
df_out_pi = df_out_pi.reset_index()
df_out_pi['CHROM'] = df_out_pi['CHROM_POS'].apply(lambda x: x.split('_')[0])
df_out_pi['POS_START'] = df_out_pi['CHROM_POS'].apply(lambda x: x.split('_')[1]).astype(int)
df_out_pi['POS_END'] = df_out_pi['CHROM_POS'].apply(lambda x: x.split('_')[2]).astype(int)
df_out_pi = df_out_pi.dropna()

print(df_out_pi)


## fst raw file の加工
# ファイルのパスを指定
vcftools_result_dir = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/vcftools_result_window10000'

fst_files = [
    'fst_10000_Hchi_Hmic.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Hsau.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Lruf.txt.windowed.weir.fst.gz',
    'fst_10000_Hchi_Pnye.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Ppun.txt.windowed.weir.fst.gz'
]

top_threshold = 0.99

# 空のリストを作成してデータフレームを格納
df_list_fst = []
df_out_fst = pd.DataFrame()
    
for fst_file in fst_files:
    prefix = os.path.splitext(fst_file)[0].split('.')[0].replace('_10000','')
    
    fst_path = os.path.join(vcftools_result_dir, fst_file)
    df_tmp_fst = pd.read_csv(fst_path, sep='\t')

    # CHROM_POS列の作成し、indexに設定
    df_tmp_fst['CHROM_POS'] = df_tmp_fst['CHROM'] + '_' + df_tmp_fst['BIN_START'].astype(str) + '_' + df_tmp_fst['BIN_END'].astype(str)
    df_tmp_fst = df_tmp_fst.set_index('CHROM_POS')

    # fstの値が1を取るものが上位1％以上存在するか確認
    fst_one_ratio = (df_tmp_fst['WEIGHTED_FST'] == 1).sum() / len(df_tmp_fst)
    if fst_one_ratio >= (1 - top_threshold):
        df_tmp_fst[prefix+'_top1%'] = (df_tmp_fst['WEIGHTED_FST'] == 1).astype(int)
    else:
        # 上位1％の閾値を計算
        threshold = df_tmp_fst['WEIGHTED_FST'].quantile(top_threshold)
        df_tmp_fst[prefix+'_top1%'] = (df_tmp_fst['WEIGHTED_FST'] >= threshold).astype(int)

    # カラムを整形
    df_tmp_fst = df_tmp_fst.rename(columns={'WEIGHTED_FST': prefix})
    df_tmp_fst = df_tmp_fst[[prefix, prefix+'_top1%']]
    
    df_list_fst.append(df_tmp_fst)

df_out_fst = pd.concat(df_list_fst, axis=1)

# 重複する列名を処理
df_out_fst = df_out_fst.loc[:,~df_out_fst.columns.duplicated()]
df_out_fst = df_out_fst[['fst_Hchi_Hmic', 'fst_Hchi_Hsau', 'fst_Hchi_Lruf', 'fst_Hchi_Pnye', 'fst_Hchi_Ppun', 
                         'fst_Hchi_Hmic_top1%', 'fst_Hchi_Hsau_top1%', 'fst_Hchi_Lruf_top1%', 'fst_Hchi_Pnye_top1%', 'fst_Hchi_Ppun_top1%']]

# CHROM, POS列の作成
df_out_fst = df_out_fst.reset_index()
df_out_fst['CHROM'] = df_out_fst['CHROM_POS'].apply(lambda x: x.split('_')[0])
df_out_fst['POS_START'] = df_out_fst['CHROM_POS'].apply(lambda x: x.split('_')[1]).astype(int)
df_out_fst['POS_END'] = df_out_fst['CHROM_POS'].apply(lambda x: x.split('_')[2]).astype(int)
df_out_fst = df_out_fst.dropna()

print(df_out_fst)


# %%
fst_files = [
    'fst_10000_Hchi_Hmic.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Hsau.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Lruf.txt.windowed.weir.fst.gz',
    'fst_10000_Hchi_Pnye.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Ppun.txt.windowed.weir.fst.gz'
]


for fst_file in fst_files:

    prefix = fst_file.split('.')[0].replace('_10000','')
    ax1, ax2 = upsetplot_fst_with_violin(data=df_out_fst, chrom="CHROM", pos="POS_START", fst=f"{prefix}", max_value_of_yaxis=1)
    ax1.set_title(prefix,fontsize=20)

    output_file = f'/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/result/upset/fst_manhattan_plot_dataset1_{prefix}_10000.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    plt.show()


# %%
top_columns = ['fst_Hchi_Hmic_top1%', 'fst_Hchi_Hsau_top1%', 'fst_Hchi_Lruf_top1%', 'fst_Hchi_Pnye_top1%', 'fst_Hchi_Ppun_top1%']

# 各比較におけるtop 1%のwindowsを示すデータフレームを作成
indicators = df_out_fst[top_columns]

# 全てのtop1%カラムで0のものを除外
indicators = indicators[indicators.sum(axis=1) > 0]

# カラム名を簡略化
indicators.columns = [col.replace('fst_Hchi_', '').replace('_top1%', '') for col in indicators.columns]

# 各行について、どのカテゴリに属しているかをリストとして取得し、重複を除去
membership_dict = {}
for _, row in indicators.iterrows():
    # その行が属するカテゴリのリストを作成
    categories = tuple(sorted([col for col, val in row.items() if val == 1]))
    if categories in membership_dict:
        membership_dict[categories] += 1
    else:
        membership_dict[categories] = 1

# membership_dictからリストを作成
memberships = []
counts = []
for categories, count in membership_dict.items():
    memberships.append(list(categories))
    counts.append(count)

# upset plotを作成
plt.figure(figsize=(10, 6))
plot(from_memberships(memberships, data=counts), 
     sort_by='-degree',
     show_counts=True,
     subset_size='sum')

plt.title('Overlap of top 1% Fst windows across comparisons')
plt.rcParams["svg.fonttype"] = "none"
plt.savefig('/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/result/upset/upset_HighFst_dataset1_10000.svg', dpi=300, bbox_inches='tight')


# データの要約を表示（オプション）
print("\nNumber of windows in each comparison:")
for col in indicators.columns:
    print(f"{col}: {(indicators[col] == 1).sum()}")

# 除外した行数を表示
print(f"\nNumber of windows excluded (all zeros): {len(df_out_fst) - len(indicators)}")

# %%

# top 1%のカラムとFst値のカラムを選択
top_columns = [col for col in df_out_fst.columns if 'top1%' in col]
fst_columns = [col for col in df_out_fst.columns if 'fst' in col and 'top1%' not in col and 'CHROM' not in col and 'POS' not in col]

# 全てのtop1%カラムで0のものを除外
indicators = df_out_fst[top_columns]
df_filtered = df_out_fst[indicators.sum(axis=1) > 0]
indicators = indicators[indicators.sum(axis=1) > 0]

# カテゴリの組み合わせと対応するFst値を収集
combinations = []
fst_values = []

for idx, row in indicators.iterrows():
    # カテゴリの組み合わせを取得
    categories = tuple(sorted([col.replace('fst_Hchi_', '').replace('_top1%', '') 
                             for col, val in row.items() if val == 1]))
    
    # 対応するFst値の平均を計算
    fst_vals = [df_filtered.loc[idx, 'fst_Hchi_Hmic':'fst_Hchi_Ppun']]
    avg_fst = np.mean(fst_vals)
    
    combo_str = '&'.join(categories)
    combinations.append(combo_str)
    fst_values.append(avg_fst)

# データフレームを作成
plot_df = pd.DataFrame({
    'Combination': combinations,
    'Average_Fst': fst_values
})

order = [
    'Hmic&Hsau&Lruf&Pnye&Ppun',
    'Hmic&Lruf&Pnye&Ppun',
    'Hmic&Hsau&Lruf&Ppun',
    'Hmic&Hsau&Pnye&Ppun',
    'Hmic&Hsau&Lruf&Pnye',
    'Hsau&Lruf&Pnye&Ppun',
    'Hmic&Lruf&Ppun',
    'Hmic&Pnye&Ppun',
    'Hmic&Hsau&Ppun',
    'Hmic&Lruf&Pnye',
    'Hmic&Hsau&Lruf',
    'Hmic&Hsau&Pnye',
    'Lruf&Pnye&Ppun',
    'Hsau&Lruf&Ppun',
    'Hsau&Pnye&Ppun',
    'Hsau&Lruf&Pnye',
    'Hmic&Ppun',
    'Hmic&Lruf',
    'Hmic&Pnye',
    'Hmic&Hsau',    
    'Lruf&Ppun',
    'Pnye&Ppun',
    'Hsau&Ppun',
    'Lruf&Pnye',
    'Hsau&Lruf',
    'Hsau&Pnye',
    'Hmic',
    'Ppun',
    'Lruf',
    'Pnye',
    'Hsau'
]

color = ['#D6617F','#D6617F','#D6617F','#D6617F','#D6617F','#D6617F',
         '#cccccc','#cccccc','#cccccc','#cccccc','#cccccc',
         '#cccccc','#cccccc','#cccccc','#cccccc','#cccccc',
         '#cccccc','#cccccc','#cccccc','#cccccc','#cccccc',
         '#cccccc','#cccccc','#cccccc','#cccccc','#cccccc',
         '#cccccc','#cccccc','#cccccc','#cccccc','#cccccc']

# プロットの作成
plt.figure(figsize=(10, 4))
sns.violinplot(data=plot_df, x='Combination', y='Average_Fst', order=order, palette=color, inner=None)
plt.xticks(rotation=90, ha='right')
plt.title('Distribution of Average Fst Values by Category Combination')
plt.xlabel('Category Combinations')
plt.ylabel('Average Fst')
plt.tight_layout()
plt.rcParams["svg.fonttype"] = "none"
plt.savefig('/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/result/upset/violin_HighFst_forupset_dataset1_10000.svg', dpi=300, bbox_inches='tight')
plt.show()

# 各組み合わせの要約統計量を表示
print("\nSummary statistics for each combination:")
summary_stats = plot_df.groupby('Combination')['Average_Fst'].describe()
# 定義した順序でソート
summary_stats = summary_stats.reindex(order)
print(summary_stats)

# %% [markdown]
# ## Visualize manhattan plot using data of fst, pi, and RNAseq coverage

# %%
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import gffutils
import numpy as np
import os
from matplotlib.lines import Line2D
import math

def plot_combined_region_with_coverage(df_fst, df_pi, df_coverage, region, species_list_fst, species_list_pi, 
                        color_list_fst, color_list_pi, color_list_hch_RNAseq, color_list_hsa_RNAseq, db, mergin):
    """
    Parameters:
    -----------
    df_fst : pandas.DataFrame
    df_pi : pandas.DataFrame
    df_coverage : pandas.DataFrame
    region : str
    species_list_fst : list
    species_list_pi : list
    color_list_fst : list
    color_list_pi : list
    color_list_hch_RNAseq : list
    color_list_hsa_RNAseq : list
    db : gffutils.FeatureDB
    """
    # フォント設定
    plt.rcParams['svg.fonttype'] = 'none'
    
    # 領域情報の解析
    chrom, start, end = region.split('_')
    start, end = max(0, int(start)-mergin), int(end)+mergin
    
    # POSを整数型に変換
    df_fst['POS_START'] = pd.to_numeric(df_fst['POS_START'], errors='coerce')
    df_pi['POS_START'] = pd.to_numeric(df_pi['POS_START'], errors='coerce')
    
    # 指定された領域のデータを抽出
    df_region_fst = df_fst[(df_fst['CHROM'] == chrom) & 
                          (df_fst['POS_START'] >= start) & 
                          (df_fst['POS_START'] <= end)]
    df_region_pi = df_pi[(df_pi['CHROM'] == chrom) & 
                        (df_pi['POS_START'] >= start) & 
                        (df_pi['POS_START'] <= end)]
    
    # カバレッジデータを抽出（filter_region関数を利用）
    df_region_coverage = filter_region(df_coverage, chrom, start, end)
    
    # カバレッジデータの位置情報の中間点を計算
    if not df_region_coverage.empty:
        df_region_coverage['position'] = (pd.to_numeric(df_region_coverage['start']) + 
                                         pd.to_numeric(df_region_coverage['end'])) / 2
    
    # Hch列とHsa列を特定
    hch_cols = [col for col in df_region_coverage.columns if col.startswith('Hch')]
    hsa_cols = [col for col in df_region_coverage.columns if col.startswith('Hsa')]

    
    # すべてのカバレッジデータの最大値を見つけて共通のy軸スケールを設定
    if not df_region_coverage.empty:
        all_data_cols = hch_cols + hsa_cols
        all_max_values = []
        
        for col in all_data_cols:
            if col in df_region_coverage.columns:
                max_val = pd.to_numeric(df_region_coverage[col]).max()
                all_max_values.append(max_val)
        
        if all_max_values:
            global_max = np.ceil(max(all_max_values))
            
            # y軸目盛りを設定
            if global_max <= 10:
                yticks = np.linspace(0, global_max, 3)
            elif global_max <= 50:
                yticks = np.linspace(0, global_max, 4)
            else:
                yticks = np.linspace(0, global_max, 5)
            
            # 整数値に丸める
            yticks = np.round(yticks).astype(int)
    
    # プロットの作成（6段組み）
    fig, (ax1, ax_top1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(7, 1, figsize=(15, 12), 
                                                      height_ratios=[1, 1.5, 4, 4, 4, 4, 4], 
                                                      sharex=True)
    
    # 上段: 遺伝子座の情報
    for gene in db.region(seqid=chrom, start=start, end=end, featuretype='gene'):
        gene_start = max(gene.start, start)
        gene_end = min(gene.end, end)
        rect = Rectangle((gene_start, 0), gene_end - gene_start, 1, 
                        fill=True, facecolor='lightgrey', edgecolor='none')
        ax1.add_patch(rect)
        
        # エキソンを濃いグレーで表示
        for exon in db.children(gene, featuretype='exon'):
            exon_start = max(exon.start, start)
            exon_end = min(exon.end, end)
            if exon_start < end and exon_end > start:
                rect = Rectangle((exon_start, 0), exon_end - exon_start, 1,
                               fill=True, facecolor='darkgrey', edgecolor='none')
                ax1.add_patch(rect)
        
        gene_name = gene.attributes.get('Name', [gene.id.split(':')[1]])[0]
        ax1.text((gene_start + gene_end) / 2, 0.5, gene_name, 
                 ha='center', va='center', fontsize=8)
    
    ax1.set_ylim(0, 1)
    ax1.set_yticks([])
    ax1.axis('off')

    # 二段目: Fst top 1%の領域を四角で表示
    for i, species in enumerate(species_list_fst):
        top1_column = f'{species}_top1%'
        top1_regions = df_region_fst[df_region_fst[top1_column] == 1]

        for _, row in top1_regions.iterrows():
            rect = Rectangle((row['POS_START']-5000, i), 10000, 1,
                           fill=True, facecolor='black', alpha=0.2)
            ax_top1.add_patch(rect)

    ax_top1.set_ylim(0, len(species_list_fst))
    ax_top1.set_yticks(np.arange(len(species_list_fst)) + 0.5)
    ax_top1.set_yticklabels(species_list_fst)
    ax_top1.tick_params(axis='y', which='both', left=False, right=False)
    ax_top1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_top1.spines['top'].set_visible(False)
    ax_top1.spines['right'].set_visible(False)
    ax_top1.spines['bottom'].set_visible(False)
    ax_top1.spines['left'].set_visible(False)
    
    # 三段目: Fstマンハッタンプロット
    for species, color in zip(species_list_fst, color_list_fst):
        x = df_region_fst['POS_START']
        y = df_region_fst[species]
        ax2.scatter(x, y, c=color, s=10, label=species)
    
    ax2.set_ylabel('FST')
    ax2.grid(True, linestyle='--', alpha=0.3)
    ax2.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    
    # 四段目: Piマンハッタンプロット
    for species, color in zip(species_list_pi, color_list_pi):
        x = df_region_pi['POS_START']
        y = df_region_pi[species]
        ax3.scatter(x, y, c=color, s=10, label=species)
    
    ax3.set_ylabel('π')
    ax3.grid(True, linestyle='--', alpha=0.3)
    ax3.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    
    # 五段目: Hchサンプルのカバレッジプロット
    if not df_region_coverage.empty:
        for i, col in enumerate(hch_cols):
            if col in df_region_coverage.columns:
                y_data = pd.to_numeric(df_region_coverage[col])
                ax4.plot(df_region_coverage['position'], y_data, color=color_list_hch_RNAseq[i], 
                        label=col, linewidth=0.8)
        
        ax4.set_ylabel('Coverage (Hch)')
        ax4.set_yticks(yticks)
        ax4.set_yticklabels([str(int(y)) for y in yticks])
        ax4.set_ylim(0, global_max * 1.05)
        ax4.legend(loc='upper right', bbox_to_anchor=(1.15, 1), fontsize=8)
        ax4.grid(True, linestyle='--', alpha=0.3)
    
    # 六段目: Hsaサンプルのカバレッジプロット
    if not df_region_coverage.empty:
        for i, col in enumerate(hsa_cols):
            if col in df_region_coverage.columns:
                y_data = pd.to_numeric(df_region_coverage[col])
                ax5.plot(df_region_coverage['position'], y_data, color=color_list_hsa_RNAseq[i], 
                        label=col, linewidth=0.8)
        
        ax5.set_ylabel('Coverage (Hsa)')
        ax5.set_yticks(yticks)
        ax5.set_yticklabels([str(int(y)) for y in yticks])
        ax5.set_ylim(0, global_max * 1.05)
        ax5.legend(loc='upper right', bbox_to_anchor=(1.15, 1), fontsize=8)
        ax5.grid(True, linestyle='--', alpha=0.3)

    if not df_region_coverage.empty:
        # 発現差のプロット
        df_region_coverage['log2FC'] = df_region_coverage[hch_cols].mean(axis=1).apply(lambda x: math.log2(x+1)) - df_region_coverage[hsa_cols].mean(axis=1).apply(lambda x: math.log2(x+1))
        ax6.bar(df_region_coverage['position'], df_region_coverage['log2FC'], width=100, color='black')    
        ax6.set_xlabel('Genomic Position')
        ax6.set_ylabel('log2 Fold Change')
        ax6.grid(True, linestyle='--', alpha=0.3)

    
    # x軸の設定
    for ax in [ax2, ax3, ax4, ax5]:
        ax.set_xlim(start, end)
        ax.set_xticks(np.linspace(start, end, 6))
        ax.set_xticklabels([f'{int(x):,}' for x in np.linspace(start, end, 6)])
    
    plt.tight_layout()
    
    # プロットの保存
    output_file = f'/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/result/vcftools_fst_pi_coverage_dataset1_{region}.svg'
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

def filter_region(df, chrom, start, end):
    """
    Parameters:
    -----------
    df : DataFrame
    chrom : str
    start : int
    end : int
        
    Returns:
    --------
    DataFrame: filter vcf
    """
    # 列名を確認
    chr_col = 'chr'
    start_col = 'start'
    end_col = 'end'
        
    # データ型を変換して確実に比較できるようにする
    df[start_col] = pd.to_numeric(df[start_col], errors='coerce')
    df[end_col] = pd.to_numeric(df[end_col], errors='coerce')
    
    # フィルタリング
    region_data = df[(df[chr_col] == chrom) & 
                     (df[start_col] >= start) & 
                     (df[end_col] <= end)]
    
    # データが存在しない場合は空のDataFrameを返す
    if region_data.empty:
        print(f"警告: 領域 {chrom}:{start}-{end} にデータが見つかりませんでした。")
    
    return region_data


# %%
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

### pi raw file の加工
# データをDataFrameに読み込む
vcftools_result_dir = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/vcftools_result_window10000'

pi_files = [
    'pi_10000_Hchi.txt.windowed.pi.gz', 'pi_10000_Hsau.txt.windowed.pi.gz', 'pi_10000_Lruf.txt.windowed.pi.gz',
    'pi_10000_Hmic.txt.windowed.pi.gz', 'pi_10000_Pnye.txt.windowed.pi.gz', 'pi_10000_Ppun.txt.windowed.pi.gz'
]

df_list_pi = []
df_out_pi = pd.DataFrame()

for pi_file in pi_files:

    # pathとprefixの定義
    prefix = os.path.splitext(pi_file)[0].split('.')[0].replace('_10000','')
    pi_path = os.path.join(vcftools_result_dir, pi_file)
    
    df_tmp_pi = pd.read_csv(pi_path, sep='\t')

    # CHROM_POS列の作成
    df_tmp_pi['CHROM_POS'] = df_tmp_pi['CHROM'] + '_' + df_tmp_pi['BIN_START'].astype(str) + '_' + df_tmp_pi['BIN_END'].astype(str)
    df_tmp_pi = df_tmp_pi.set_index('CHROM_POS')

    # カラムを整形
    df_tmp_pi = df_tmp_pi.rename(columns={'PI': prefix})
    df_tmp_pi = df_tmp_pi[prefix]

    df_list_pi.append(df_tmp_pi)


# concat関数を使用して横軸に連結
df_out_pi = pd.concat(df_list_pi, axis=1)

# 重複する列名を処理
df_out_pi = df_out_pi.loc[:,~df_out_pi.columns.duplicated()]
df_out_pi = df_out_pi[['pi_Hchi', 'pi_Hsau', 'pi_Lruf', 'pi_Hmic', 'pi_Ppun', 'pi_Pnye']]

# CHROM, POS列の作成
df_out_pi = df_out_pi.reset_index()
df_out_pi['CHROM'] = df_out_pi['CHROM_POS'].apply(lambda x: x.split('_')[0])
df_out_pi['POS_START'] = df_out_pi['CHROM_POS'].apply(lambda x: x.split('_')[1]).astype(int)
df_out_pi['POS_END'] = df_out_pi['CHROM_POS'].apply(lambda x: x.split('_')[2]).astype(int)
df_out_pi = df_out_pi.fillna(0) # 重要：集団内で固定しているサイトはpiが計算できない。そのためここを0で埋めた

## fst raw file の加工
# ファイルのパスを指定
vcftools_result_dir = '/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/vcftools_result_window10000'

fst_files = [
    'fst_10000_Hchi_Hmic.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Hsau.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Lruf.txt.windowed.weir.fst.gz',
    'fst_10000_Hchi_Pnye.txt.windowed.weir.fst.gz', 'fst_10000_Hchi_Ppun.txt.windowed.weir.fst.gz'
]

top_threshold = 0.99

# 空のリストを作成してデータフレームを格納
df_list_fst = []
df_out_fst = pd.DataFrame()
    
for fst_file in fst_files:
    prefix = os.path.splitext(fst_file)[0].split('.')[0].replace('_10000','')
    
    fst_path = os.path.join(vcftools_result_dir, fst_file)
    df_tmp_fst = pd.read_csv(fst_path, sep='\t')

    # CHROM_POS列の作成し、indexに設定
    df_tmp_fst['CHROM_POS'] = df_tmp_fst['CHROM'] + '_' + df_tmp_fst['BIN_START'].astype(str) + '_' + df_tmp_fst['BIN_END'].astype(str)
    df_tmp_fst = df_tmp_fst.set_index('CHROM_POS')

    # fstの値が1を取るものが上位1％以上存在するか確認
    fst_one_ratio = (df_tmp_fst['WEIGHTED_FST'] == 1).sum() / len(df_tmp_fst)
    if fst_one_ratio >= (1 - top_threshold):
        df_tmp_fst[prefix+'_top1%'] = (df_tmp_fst['WEIGHTED_FST'] == 1).astype(int)
    else:
        # 上位1％の閾値を計算
        threshold = df_tmp_fst['WEIGHTED_FST'].quantile(top_threshold)
        df_tmp_fst[prefix+'_top1%'] = (df_tmp_fst['WEIGHTED_FST'] >= threshold).astype(int)

    # カラムを整形
    df_tmp_fst = df_tmp_fst.rename(columns={'WEIGHTED_FST': prefix})
    df_tmp_fst = df_tmp_fst[[prefix, prefix+'_top1%']]    
    df_list_fst.append(df_tmp_fst)

df_out_fst = pd.concat(df_list_fst, axis=1)

# 重複する列名を処理
df_out_fst = df_out_fst.loc[:,~df_out_fst.columns.duplicated()]
df_out_fst = df_out_fst[['fst_Hchi_Hmic', 'fst_Hchi_Hsau', 'fst_Hchi_Lruf', 'fst_Hchi_Pnye', 'fst_Hchi_Ppun',
                         'fst_Hchi_Hmic_top1%', 'fst_Hchi_Hsau_top1%', 'fst_Hchi_Lruf_top1%', 'fst_Hchi_Pnye_top1%', 'fst_Hchi_Ppun_top1%']]

# CHROM, POS列の作成
df_out_fst = df_out_fst.reset_index()
df_out_fst['CHROM'] = df_out_fst['CHROM_POS'].apply(lambda x: x.split('_')[0])
df_out_fst['POS_START'] = df_out_fst['CHROM_POS'].apply(lambda x: x.split('_')[1]).astype(int)
df_out_fst['POS_END'] = df_out_fst['CHROM_POS'].apply(lambda x: x.split('_')[2]).astype(int)
# df_out_fst = df_out_fst.dropna()

# RNA-seqカバレッジデータを読み込む
coverage_file = '/Volumes/Transcend/250304_RNAseq/bigwig/RNAseq_Lip_Hch_Hsa_TMM_handedit.tsv'
print(f"RNA-seqカバレッジデータを読み込んでいます: {coverage_file}")
df_coverage = pd.read_csv(coverage_file, sep='\t')

# データの最初の数行を表示して確認
print("Fstデータサンプル:")
print(df_out_fst.head())
print("\nPiデータサンプル:")
print(df_out_pi.head())
print("\nRNA-seqカバレッジデータサンプル:")
print(df_coverage.head())

# %%
df_out_fst[df_out_fst[['fst_Hchi_Hmic_top1%', 'fst_Hchi_Hsau_top1%', 'fst_Hchi_Lruf_top1%', 'fst_Hchi_Pnye_top1%', 'fst_Hchi_Ppun_top1%']].sum(axis=1) >3].to_csv('/Users/machiinagatoshi/Desktop/Research_Desk/genome_analysis/250219_for_paper/fst/250519_dataset1_realanysis/result/hdr_top1%_dup4_dataset1_10000.csv')

# %%
# 種リストと色リストの設定
species_list_fst = ['fst_Hchi_Hmic', 'fst_Hchi_Hsau', 'fst_Hchi_Lruf', 'fst_Hchi_Pnye', 'fst_Hchi_Ppun']
species_list_pi = ['pi_Hsau', 'pi_Lruf', 'pi_Hmic', 'pi_Ppun', 'pi_Pnye', 'pi_Hchi']
color_list_fst = ['#03A688', '#9FC131', '#99B4BF', '#F2668B', '#253C59']
color_list_pi = ['#809A94', '#A7AB9E', '#A7ADB2', '#6B7580', '#B3A89C', '#F2668B']
hch_colors = ['#FF0000', '#FF4D6D', '#FF99AC', '#FFC2D1'] 
hsa_colors = ['#2E2E2E', '#5A5A5A', '#878787', '#C1C1C1']

regions = ['LG8_16830001_16900001']

mergin = 0

for region in regions:
    
    print(f"領域をプロットしています: {region}")
    plot_combined_region_with_coverage(
        df_out_fst, df_out_pi, df_coverage, region, 
        species_list_fst, species_list_pi, 
        color_list_fst, color_list_pi, 
        hch_colors, hsa_colors, db, mergin
    )
    print(f"プロットが完了しました: {region}")
