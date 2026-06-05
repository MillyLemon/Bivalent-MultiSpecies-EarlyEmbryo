import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# ===================================================
# 统一字体大小
# ===================================================
R = 1.2
base_font_size = 10 * R
title_font_size = 14 * R
axis_label_font_size = 12 * R
tick_label_font_size = 10 * R
legend_font_size = 11 * R
suptitle_font_size = 16 * R

plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({
    'font.size': base_font_size,
    'axes.titlesize': title_font_size,
    'axes.labelsize': axis_label_font_size,
    'xtick.labelsize': tick_label_font_size,
    'ytick.labelsize': tick_label_font_size,
    'legend.fontsize': legend_font_size,
    'figure.titlesize': suptitle_font_size,
})

# ===================================================
# 工具函数
# ===================================================
def extract_gene_id(bed_gene_id):
    if '-' in bed_gene_id:
        return bed_gene_id.rsplit('-', 1)[0]
    return bed_gene_id

def read_tpm_file(tpm_file):
    tpm_dict = {}
    try:
        df = pd.read_csv(tpm_file, sep='\t')
        if 'Name' in df.columns and 'TPM' in df.columns:
            for _, row in df.iterrows():
                tpm_dict[row['Name']] = row['TPM']
    except:
        pass
    return tpm_dict

def count_unique_genes_and_get_tpm(bed_file, tpm_dict):
    try:
        df = pd.read_csv(bed_file, sep='\t', header=None, usecols=[3], names=['gene_id'])
        df['base_gene_id'] = df['gene_id'].apply(extract_gene_id)
        unique_genes = df['base_gene_id'].unique()
        tpm_values = [np.log2(tpm_dict.get(g, 0) + 1) for g in unique_genes]
        return tpm_values
    except:
        return []

def process_sample(sample_name):
    sample_dir = sample_name
    tpm_file = os.path.join(sample_dir, 'quant.genes.sf')
    tpm_dict = read_tpm_file(tpm_file) if os.path.exists(tpm_file) else {}

    regions = {
        'Promoter': 'promoter_genes.trans.bed',
        'TSS': 'tss_genes.trans.bed',
        'Broad': 'broad_genes.trans.bed',
    }
    res = {}
    for r, f in regions.items():
        path = os.path.join(sample_dir, f)
        res[r] = count_unique_genes_and_get_tpm(path, tpm_dict)
    return res

def read_sample_names(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip()]

# ===================================================
# 显著性检验
# ===================================================
def sig_stars(p):
    if p < 0.001: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return 'ns'

def mw_test(a, b):
    if len(a) == 0 or len(b) == 0:
        return 1.0
    return mannwhitneyu(a, b, alternative='two-sided')[1]

# ===================================================
# 绘图：3个柱子 + 3组两两比较
# ===================================================
def plot_with_stats(ax, data_dict, title):
    pro = data_dict['Promoter']
    tss = data_dict['TSS']
    bro = data_dict['Broad']

    plot_data = [pro, tss, bro]
    labels = ['Promoter', 'TSS', 'Broad']
    colors = ['#00D4FF', '#FFE600', '#000080']

    bp = ax.boxplot(plot_data, labels=labels, patch_artist=True, widths=0.4,
                    medianprops=dict(color='black', linewidth=2*R))
    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)

    p1 = mw_test(pro, tss)
    p2 = mw_test(pro, bro)
    p3 = mw_test(tss, bro)

    all_vals = pro + tss + bro
    y_max = max(all_vals) if all_vals else 1
    h = y_max * 0.08

    def add_sig(i1, i2, pval, lvl):
        y = y_max + lvl * h
        ax.plot([i1, i2], [y, y], c='black', lw=1)
        ax.text((i1+i2)/2, y, sig_stars(pval), ha='center', fontsize=10*R, weight='bold')

    add_sig(1, 2, p1, 1)
    add_sig(1, 3, p2, 2)
    add_sig(2, 3, p3, 3)

    ax.set_title(title, weight='bold')
    ax.set_ylabel('log2(TPM+1)')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(top=y_max + 4*h)

# ===================================================
# 主函数：两行布局
# 第一行：H3K27me3
# 第二行：H3K4me3
# ===================================================
def main():
    samples = read_sample_names("sample.name")
    
    # 自动分成两组
    k27_samples = [s for s in samples if "K27" in s or "H3K27" in s]
    k4_samples = [s for s in samples if "K4" in s or "H3K4" in s]
    
    n_cols = max(len(k27_samples), len(k4_samples))
    fig, axes = plt.subplots(2, n_cols, figsize=(3.2*n_cols, 10), sharey=True)

    # 第一行：K27
    for i, smp in enumerate(k27_samples):
        dat = process_sample(smp)
        plot_with_stats(axes[0, i], dat, smp)

    # 第二行：K4
    for i, smp in enumerate(k4_samples):
        dat = process_sample(smp)
        plot_with_stats(axes[1, i], dat, smp)

    # 隐藏空白子图
    for row in range(2):
        for col in range(n_cols):
            if (row == 0 and col >= len(k27_samples)) or (row == 1 and col >= len(k4_samples)):
                axes[row, col].set_visible(False)

    # 行标题
    axes[0,0].text(0.02, 0.98, 'H3K27me3', transform=axes[0,0].transAxes,
                   weight='bold', fontsize=12*R, va='top')
    axes[1,0].text(0.02, 0.98, 'H3K4me3', transform=axes[1,0].transAxes,
                   weight='bold', fontsize=12*R, va='top')

    fig.suptitle('TPM Distribution (log2(TPM+1), with zero)', weight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('combined_sample_boxplots_log2_plus1_with_zero.fix.pdf', bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()
