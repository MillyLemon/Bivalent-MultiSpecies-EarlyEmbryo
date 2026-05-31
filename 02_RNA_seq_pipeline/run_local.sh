#!/bin/bash
# run_local.sh - 本地运行版本

echo "======================================"
echo "Snakemake工作流 - 本地模式"
echo "开始时间: $(date)"
echo "======================================"

# 设置环境
export PYTHONUNBUFFERED=1

# 创建目录结构
mkdir -p {logs/snakemake,logs/cluster,fastp_fq,fastp_report,salmon_output}

# 检测系统资源
AVAILABLE_CORES=$(grep -c ^processor /proc/cpuinfo)
USED_CORES=$((AVAILABLE_CORES - 2))
if [ $USED_CORES -lt 4 ]; then
    USED_CORES=4  # 最小使用4个核心
fi
if [ $USED_CORES -gt 16 ]; then
    USED_CORES=16  # 最大不超过16个核心
fi

echo "系统信息:"
echo "- CPU核心: $AVAILABLE_CORES (使用: $USED_CORES)"
echo "- 内存: $(free -h | awk '/^Mem:/{print $2}')"
echo "- 磁盘空间: $(df -h . | awk 'NR==2{print $4}') 可用"

# 运行snakemake
echo "开始执行Snakemake工作流..."
snakemake \
    --snakefile Snakefile \
    --configfile config.json \
    --cores $USED_CORES \
    --jobs $USED_CORES \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    --printshellcmds \
    2>&1 | tee logs/snakemake/run_$(date +%Y%m%d_%H%M%S).log

# 检查执行状态
EXIT_CODE=${PIPESTATUS[0]}
echo "======================================"
echo "工作流执行完成"
echo "结束时间: $(date)"
echo "退出代码: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ 所有任务成功完成！"
else
    echo "⚠ 部分任务失败或中断"
fi
echo "======================================"
