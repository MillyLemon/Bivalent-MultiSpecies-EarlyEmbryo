import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def count_unique_genes(bed_file):
    """统计bed文件中唯一的基因数量（基于第4列）"""
    try:
        # 读取bed文件，使用制表符分隔，没有表头
        df = pd.read_csv(bed_file, sep='\t', header=None, 
                        usecols=[3], names=['gene_id'])
        # 统计唯一基因数量
        return df['gene_id'].nunique()
    except Exception as e:
        print(f"读取文件 {bed_file} 时出错: {e}")
        return 0

def process_sample(sample_name, mode='gene'):
    """
    处理单个样本，根据mode选择文件类型: 'gene' 或 'trans'
    - mode='gene'  : 读取 *_genes.bed 文件
    - mode='trans' : 读取 *_genes.trans.bed 文件
    """
    sample_dir = sample_name
    suffix = '.trans' if mode == 'trans' else ''
    
    broad_file = os.path.join(sample_dir, f'broad_genes{suffix}.bed')
    promoter_file = os.path.join(sample_dir, f'promoter_genes{suffix}.bed')
    tss_file = os.path.join(sample_dir, f'tss_genes{suffix}.bed')
    
    counts = {
        'Region': ['Promoter', 'TSS','Broad'],
        'Gene_Count': [],
        'Sample': sample_name,
        'Mode': mode
    }
    
    bed_files = [promoter_file, tss_file,broad_file]
    regions = ['Promoter', 'TSS','Broad']
    
    for bed_file, region in zip(bed_files, regions):
        if os.path.exists(bed_file):
            gene_count = count_unique_genes(bed_file)
            counts['Gene_Count'].append(gene_count)
            print(f"{sample_name} ({mode}) - {region}: {gene_count}")
        else:
            counts['Gene_Count'].append(0)
            print(f"警告: {bed_file} 不存在")
    
    return counts

def plot_single_sample(ax, counts, sample_name):
    """在指定的axes上绘制单个样本的条形图"""
    colors = ['#00D4FF','#FFE600','#000080']  # 为三个区域设置不同的颜色
    
    bars = ax.bar(counts['Region'], counts['Gene_Count'], 
                 color=colors, alpha=0.8, edgecolor='black')
    
    # 添加数值标签
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height):,}', ha='center', va='bottom', fontsize=9)
    
    # 设置标题和标签
    ax.set_title(sample_name, fontsize=11, fontweight='bold')
    ax.set_ylabel('Gene Count', fontsize=10)
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    
    # 添加网格线
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    
    # 设置y轴从0开始
    max_count = max(counts['Gene_Count']) if counts['Gene_Count'] else 1
    ax.set_ylim(0, max_count * 1.15)
    
    return ax

def create_gene_count_plot(all_data, output_filename):
    """
    根据所有样本的数据创建多子图网格，并保存为PDF
    all_data: 由process_sample返回的字典列表
    output_filename: 输出PDF文件名
    """
    n_samples = len(all_data)
    # 固定2行6列布局（假设最多12个样本，原脚本风格）
    n_rows, n_cols = 2, 6
    
    fig = plt.figure(figsize=(24, 8))
    gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.3)
    
    # 绘制每个样本的子图
    for i, data in enumerate(all_data):
        row = i // n_cols
        col = i % n_cols
        ax = fig.add_subplot(gs[row, col])
        plot_single_sample(ax, data, data['Sample'])
    
    # 隐藏多余的子图（如果样本数不足12）
    for i in range(len(all_data), n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        fig.add_subplot(gs[row, col]).set_visible(False)
    
    # 添加整体图例
    colors = ['#00D4FF','#FFE600','#000080']
    labels = ['Promoter', 'TSS','Broad']
    legend_elements = [plt.Rectangle((0,0), 1, 1, facecolor=color, edgecolor='black', alpha=0.8) 
                      for color, label in zip(colors, labels)]
    
    fig.legend(legend_elements, labels, loc='upper center', 
               ncol=3, fontsize=11, bbox_to_anchor=(0.5, 1.02))
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_filename, format='pdf', bbox_inches='tight')
    print(f"图像已保存为: {output_filename}")
    plt.show()

def read_sample_names(filename):
    """从文件中读取样本名称"""
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def main():
    # 从samples.name文件中读取样本名称
    sample_file = "sample.name"
    if not os.path.exists(sample_file):
        print(f"错误: 找不到样本文件 {sample_file}")
        print("请确保samples.name文件存在于当前目录")
        return
    
    all_samples = read_sample_names(sample_file)
    print(f"从 {sample_file} 中读取到 {len(all_samples)} 个样本")
    
    # 按H3K27me3和H3K4me3分组并排序
    h3k27me3_samples = sorted([s for s in all_samples if "H3K27me3" in s])
    h3k4me3_samples = sorted([s for s in all_samples if "H3K4me3" in s])
    sorted_samples = h3k27me3_samples + h3k4me3_samples
    
    # 收集两种模式的数据
    print("\n正在处理基因模式 (gene) 数据...")
    gene_all_data = []
    for sample_name in sorted_samples:
        gene_data = process_sample(sample_name, mode='gene')
        gene_all_data.append(gene_data)
    
    print("\n正在处理转录本模式 (trans) 数据...")
    trans_all_data = []
    for sample_name in sorted_samples:
        trans_data = process_sample(sample_name, mode='trans')
        trans_all_data.append(trans_data)
    
    if not gene_all_data or not trans_all_data:
        print("错误: 没有找到任何样本数据")
        return
    
    # 生成两个PDF图像
    create_gene_count_plot(gene_all_data, 'combined_gene_counts_plot.fix.pdf')
    create_gene_count_plot(trans_all_data, 'combined_trans_counts_plot.fix.pdf')

if __name__ == "__main__":
    main()
