#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
提取 HOMER annotatePeaks.pl 输出的 Motif 靶基因文件中的唯一基因名，
并保存为独立的基因列表文件（每行一个基因）。
"""

import pandas as pd
import os

# ================= 配置区域（按实际路径修改） =================
# 文件路径字典：转录因子名称 -> 对应的 HOMER 输出文件
FILE_PATHS = {
    'NR4A2':   'NR4A2_targets.txt',
    'PRDM1':   'PRDM1_targets.txt',
    'ZBTB18':  'ZBTB18_targets.txt',
    'ZNF135':  'ZNF135_targets.txt'
}

# 输出目录（如果不指定，则保存在当前目录）
OUTPUT_DIR = '.'   # '.' 表示当前目录
# =============================================================

def extract_unique_genes(file_path):
    """
    从 HOMER annotatePeaks.pl 输出文件中提取唯一基因名。
    文件应为制表符分隔，含有 'Gene Name' 列。
    """
    # HOMER 输出的第一行通常是注释（以 # 开头），跳过这些行
    df = pd.read_csv(file_path, sep='\t', comment='#')
    
    # 检查必需的列
    if 'Gene Name' not in df.columns:
        raise ValueError(f"文件 {file_path} 中未找到 'Gene Name' 列")
    
    # 提取基因名，去除空值和 '-' 占位符
    genes = df['Gene Name'].dropna().unique()
    genes = [g for g in genes if g != '-' and pd.notna(g)]
    
    # 去重并排序
    unique_genes = sorted(set(genes))
    return unique_genes


def save_gene_list(gene_list, output_path):
    """将基因列表保存到文本文件，每行一个基因"""
    with open(output_path, 'w', encoding='utf-8') as f:
        for gene in gene_list:
            f.write(gene + '\n')


def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    for tf_name, file_path in FILE_PATHS.items():
        if not os.path.exists(file_path):
            print(f"警告：文件 {file_path} 不存在，跳过 {tf_name}")
            continue
        
        print(f"正在处理 {tf_name} ({file_path}) ...")
        genes = extract_unique_genes(file_path)
        
        output_file = os.path.join(OUTPUT_DIR, f"{tf_name}_gene_list.txt")
        save_gene_list(genes, output_file)
        
        print(f"  -> 提取到 {len(genes)} 个唯一基因，已保存至 {output_file}")
        print(f"     前10个基因示例: {genes[:10]}")
        print("")

    print("所有文件处理完成！")

if __name__ == "__main__":
    main()
