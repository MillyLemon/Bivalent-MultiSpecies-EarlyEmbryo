#!/usr/bin/env python3
# create_transcript_bed.py

import gzip
import re
from collections import defaultdict
import sys

def extract_gene_id(attributes):
    """从GTF属性字符串中提取gene_id"""
    # 使用正则表达式匹配gene_id
    match = re.search(r'gene_id\s+"([^"]+)"', attributes)
    if match:
        return match.group(1)
    
    # 备用方法
    for item in attributes.split(';'):
        item = item.strip()
        if item.startswith('gene_id'):
            # 提取引号内的内容
            match = re.search(r'"([^"]+)"', item)
            if match:
                return match.group(1)
    
    return None

def main():
    input_gtf = sys.argv[1] if len(sys.argv) > 1 else "gencode.vM23.annotation.gtf.gz"
    output_bed = sys.argv[2] if len(sys.argv) > 2 else "geneWithDiff_trans.bed"
    
    print(f"处理文件: {input_gtf}")
    
    # 存储每个基因的转录本
    gene_transcripts = defaultdict(list)
    
    line_count = 0
    with gzip.open(input_gtf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type != 'transcript':
                continue
            
            line_count += 1
            if line_count % 10000 == 0:
                print(f"已处理 {line_count} 行...")
            
            # 提取gene_id
            gene_id = extract_gene_id(fields[8])
            if not gene_id:
                print(f"警告: 第{line_count}行无法提取gene_id")
                continue
            
            # 提取其他信息
            chrom = fields[0]
            start = int(fields[3]) - 1  # GTF 1-based → BED 0-based
            end = int(fields[4])
            strand = fields[6]
            
            # 存储转录本信息
            gene_transcripts[gene_id].append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand
            })
    
    print(f"处理完成，共找到 {len(gene_transcripts)} 个基因")
    
    # 写入BED文件
    with open(output_bed, 'w') as out_f:
        for gene_id, transcripts in sorted(gene_transcripts.items()):
            # 按起始位置排序
            transcripts.sort(key=lambda x: x['start'])
            
            # 为每个转录本编号
            for i, transcript in enumerate(transcripts, 1):
                transcript_name = f"{gene_id}-{i}"
                bed_line = f"{transcript['chrom']}\t{transcript['start']}\t{transcript['end']}\t{transcript_name}\t0\t{transcript['strand']}"
                out_f.write(bed_line + '\n')
    
    total_transcripts = sum(len(t) for t in gene_transcripts.values())
    print(f"生成文件: {output_bed}")
    print(f"总转录本数: {total_transcripts}")
    
    # 显示示例
    print("\n示例输出:")
    for gene_id in list(gene_transcripts.keys())[:3]:
        print(f"\n基因 {gene_id} 有 {len(gene_transcripts[gene_id])} 个转录本:")
        transcripts = gene_transcripts[gene_id]
        for i, t in enumerate(transcripts[:3], 1):
            print(f"  {gene_id}-{i}: {t['chrom']}:{t['start']}-{t['end']} ({t['strand']})")

if __name__ == "__main__":
    main()
