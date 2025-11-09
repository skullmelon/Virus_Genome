#下载获取样本测序文件（本步骤文件可替换为样本测序文件）
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10971381/SRR10971381 -O SRR10971381.sra

#将测序文件拆分为fastq格式并压缩
fasterq-dump SRR10971381.sra
nohup gzip *.fastq

#下载centrifuge的hvc数据库（包含人类病毒和新冠病毒）并解压缩
mkdir hvc
wget https://zenodo.org/record/3732127/files/h+v+c.tar.gz?download=1
tar -zxvf h+v+c.tar.gz ./

#centrifuge比对
centrifuge -x ~/04.ncov/28.data/database/hvc/hvc -1 ~/04.ncov/28.data/sra/SRR10971381.sra_1.fastq.gz -2 ~/04.ncov/28.data/sra/SRR10971381.sra_2.fastq.gz -S result.tsv --report-file report.tsv -p 12 >centrifuge.log
#centrifuge -x 索引数据库 -1 样本的fastq（或其压缩包）第一条 -2 样本的fastq（或其压缩包）第二条 -S 结果输入 --report-file report.tsv -p 12 >centrifuge.log

# 初始化日志文件
echo "开始处理文件: $(date)" > workflow.log
echo "----------------------------------------" >> workflow.log

# 检查输入文件是否存在
if [ ! -f "report.tsv" ]; then
    echo "错误: 文件 report.tsv 不存在!" >> workflow.log
    echo "错误: 文件 report.tsv 不存在! 请检查文件是否存在。"
    exit 1
fi

# 验证第6列列名
header=$(head -n 1 report.tsv)
col6_name=$(echo "$header" | awk -F'\t' '{print $6}')
if [ "$col6_name" != "numUniqueReads" ]; then
    echo "错误: 第6列不是numUniqueReads，实际是:$col6_name" >> workflow.log
    echo "错误: 列名错位，请检查表头"
    exit 1
fi
echo "列名验证通过：第6列确为numUniqueReads" >> workflow.log

# 验证列数
header_cols=$(echo "$header" | awk -F'\t' '{print NF}')
if [ "$header_cols" -ne 7 ]; then
    echo "错误: 表头列数异常，预期7列，实际$header_cols列" >> workflow.log
    echo "错误: 表头格式错误，列数不符合预期"
    exit 1
fi
echo "表头列数检查通过，共7列" >> workflow.log

# 提取表头
echo "$header" > header.tmp
echo "已提取表头信息" >> workflow.log

# 处理数据行（清理+验证）
tail -n +2 report.tsv | awk -F'\t' '
{
    # 清理第6列非数字字符（防御性处理）
    gsub(/[^0-9]/, "", $6)
    if ($6 ~ /^[0-9]+$/) {
        print $0
    } else {
        print "警告: 行"NR+1"的numUniqueReads无效，已跳过" >> "workflow.log"
    }
}' OFS='\t' > data_cleaned.tmp

# 关键修复：明确指定制表符分隔符，强制数值排序
# 筛选 numUniqueReads < 100 并排序
cat header.tmp > less100.tsv
awk -F'\t' '$6 + 0 < 100' data_cleaned.tmp | sort -t$'\t' -k6,6nr >> less100.tsv
echo "已生成 less100.tsv（按numUniqueReads降序）" >> workflow.log

# 筛选 numUniqueReads >= 100 并排序
cat header.tmp > more100.tsv
awk -F'\t' '$6 + 0 >= 100' data_cleaned.tmp | sort -t$'\t' -k6,6nr >> more100.tsv
echo "已生成 more100.tsv（按numUniqueReads降序）" >> workflow.log

# 检查特定病毒名称
echo "开始检查 more100.tsv 中的特定名称..." >> workflow.log
if grep -q "Severe acute respiratory syndrome-related coronavirus" more100.tsv; then
    echo "阳性"
    echo "检查结果: 阳性" >> workflow.log
else
    echo "阴性"
    echo "检查结果: 阴性" >> workflow.log
fi

# 清理临时文件
rm -f header.tmp data_cleaned.tmp

echo "----------------------------------------" >> workflow.log
echo "处理完成: $(date)" >> workflow.log
echo "处理详情已记录到 workflow.log"
