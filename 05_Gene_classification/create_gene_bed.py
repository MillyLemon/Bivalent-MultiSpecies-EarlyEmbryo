import re
import gzip

def parse_gtf_attributes(attr_str):
    """解析GTF属性字符串"""
    attributes = {}
    # 移除结尾的分号，然后按分号分割
    attr_str = attr_str.strip(';')
    # 查找所有键值对
    pattern = r'(\w+)\s+"([^"]*)"'
    matches = re.findall(pattern, attr_str)
    
    for key, value in matches:
        attributes[key] = value
    
    return attributes

def gtf_to_gene_bed(gtf_file, output_bed):
    """从GTF文件中提取gene信息并保存为BED格式"""
    
    processed_genes = set()  # 记录已处理的gene，避免重复
    
    # 尝试用gzip打开，如果失败则用普通方式打开
    try:
        if gtf_file.endswith('.gz'):
            # 使用gzip打开
            with gzip.open(gtf_file, 'rt') as infile, open(output_bed, 'w') as outfile:
                process_file(infile, outfile, processed_genes)
        else:
            # 普通文本文件
            with open(gtf_file, 'r') as infile, open(output_bed, 'w') as outfile:
                process_file(infile, outfile, processed_genes)
    except UnicodeDecodeError:
        # 如果是压缩文件但没有.gz后缀
        with gzip.open(gtf_file, 'rt') as infile, open(output_bed, 'w') as outfile:
            process_file(infile, outfile, processed_genes)
    
    print(f"成功提取 {len(processed_genes)} 个基因")
    print(f"输出文件: {output_bed}")

def process_file(infile, outfile, processed_genes):
    """处理文件内容"""
    for line in infile:
        # 跳过注释行
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        
        # 确保有足够的字段
        if len(fields) < 9:
            continue
        
        # 只处理gene类型的行
        if fields[2] != 'gene':
            continue
        
        # 获取基本信息
        chrom = fields[0]
        start = int(fields[3]) - 1  # GTF是1-based，BED是0-based
        end = int(fields[4])
        strand = fields[6]
        
        # 解析属性字段
        attributes = parse_gtf_attributes(fields[8])
        
        # 提取gene_id（移除版本号）
       # gene_id = attributes.get('gene_id', '').split('.')[0]
       
	#保留版本号
        gene_id = attributes.get('gene_id', '')
	# 检查是否已处理过这个gene
        if gene_id in processed_genes:
            continue
        
        # 写入BED格式（6列）
        bed_line = f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n"
        outfile.write(bed_line)
        
        processed_genes.add(gene_id)

# 使用示例
if __name__ == "__main__":
    gtf_file = "gencode.vM23.annotation.gtf.gz"  # 注意：可能是压缩文件
    output_bed = "gencode.vM23.annotation.genes.bed"         # 输出BED文件
    
    gtf_to_gene_bed(gtf_file, output_bed)
