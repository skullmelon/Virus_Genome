#本脚本为新冠病毒基因组拼接自动流程脚本

#1.下载新冠病毒基因组原始测序数据数据
mkdir illumina_data
cd illumina_data
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10971381/SRR10971381 -O SRR10971381.sra    #下载链接可能发生改变，在NCBI官网搜索SRR10971381获取最新URL

#2.将illumina测序文件sra拆分为两条正反义链的fastq文件并压缩
fasterq-dump SRR10971381.sra
gzip *.fastq

#3.对测序数据进行第一次质量检测
mkdir qc_1
fastqc -f fastq -o ./qc_1 SRR10971381.sra_1.fastq.gz SRR10971381.sra_2.fastq.gz
sz ./qc_1
        #运行完毕后可在本地端查看文件，判断质量。
while true; do
    echo -e "\n请选择后续执行路径："
    echo "1) 进行质量控制"
    echo "2) 直接进行下一步"
    echo -n "请输入数字 1 或 2（默认选择1）：" 
    read -r user_choice

    # case语句匹配输入（支持多字符匹配，比如1或one都触发路径1）
    case "$user_choice" in
        1)  
            echo -e "\n执行路径1：进行质量控制"
            
            #用户确定各参数
            while true; do
            echo -n "请输入严格程度（参数-z ，默认4）："
            read -r param_z
            # 若用户未输入，使用默认值4
            param_z=${param_z:-4}

            echo -n "请输入质量阈值（参数-q ，默认20）："
            read -r param_q
            # 若用户未输入，使用默认值20
            param_q=${param_q:-20}

            echo -n "请输入允许的低质量碱基比例，低于此比例将丢弃（参数-u ，默认30）："
            read -r param_u
            # 若用户未输入，使用默认值30
            param_u=${param_u:-30}

            echo -n "请输入允许的未知碱基N数量（-n 参数，默认10）："
            read -r param_n
            # 若用户未输入，使用默认值10
            param_n=${param_n:-10}

            # 确认命令
            echo -e "\n各参数设定如下：-z=$param_z, -q=$param_q, -u=$param_u, -n=$param_n"
            echo -n "请选择：继续执行(y/Y) 或 重新设定(r/R)："

            while true; do
            read -r user_choice
            case "$user_choice" in
                y|Y)  # 选择继续，跳出外层循环
                    echo -e "\n用户确认继续，开始进行质量控制."
                    break 2
                    ;;
                r|R)  # 选择重新设定，重新执行外层循环（重新询问参数）
                    echo -e "\n重新重新设定参数"
                    break 1
                    ;;
                *)  # 无效输入，提示后重新让用户选择
                    echo -n "输入无效，请输入 y/Y（继续）或 r/R（重新设定）："
                    ;;
            esac
        done
    done

            # 最终命令（固定参数 -i、-I、-o、-O、-h 保持不变）
            echo -n "fastp -i SRR10971381.sra_1.fastq.gz -I SRR10971381.sra_2.fastq.gz -o clean.1.fq.gz -O clean.2.fq.gz -z $param_z -q $param_q -u $param_u -n $param_n -h clean.html"
            fastp -i SRR10971381.sra_1.fastq.gz -I SRR10971381.sra_2.fastq.gz -o clean.1.fq.gz -O clean.2.fq.gz -z $param_z -q $param_q -u $param_u -n $param_n -h clean.html

            #将过滤后的序列移动至新文件夹
            mkdir filter
            mv -f clean.1.fq.gz clean.2.fq.gz clean.html fastp.json ./filter

            #进行第二次质量检测
            cd filter
            mkdir qc_2
            fastqc -f fastq -o qc_2 clean.1.fq.gz clean.2.fq.gz
            sz ./qc_2
            break
            ;;
        
        2)  
            mv SRR10971381.sra_1.fastq.gz clean.1.fq.gz
            mv SRR10971381.sra_2.fastq.gz clean.2.fq.gz
            echo -e "\n执行路径2：直接进行下一步"
            break
            ;;

        *)  # 无效输入
            echo -e "\n错误：输入必须是 1 或 2。"
            echo -e "重新选择执行路径"
            ;;
    esac

echo -e "质控过滤结束，开始拼接"

#该步骤有两个选择，选择1：通过refgenie直接下载bwa库中的已经建好的索引文件，方便性能较低的服务器。选择2：:直接下载人基因组并重新构建索引库。
#二选一即可，默认为选择1.
#选择1
#安装refgenie软件：
mamba install -y refgenie
#创建refgenie软件配置文件存放目录：
mkdir -p ~/refgenie_config
#初始refgenie化软件并指定配置文件目录：
refgenie init -c ~/refgenie_config/genome_config.yaml
#显示软件库中所有索引文件：
echo -n "显示软件库中所有索引文件"
refgenie listr
echo -n -e "\n请输入索引文件名称（bwa格式对应的基因组名，默认hg38）："
read -r bwa_index
# 若用户未输入，使用默认值hg38
bwa_index=${bwa_index:-hg38} 
# 执行下载（修复原代码变量引用错误：加$引用变量）
echo -e "\n开始下载索引：$bwa_index/bwa_index"
refgenie pull "$bwa_index/bwa_index" -c ~/refgenie_config/genome_config.yaml
# 简单判断下载结果
if [ $? -eq 0 ]; then
    echo -e "\n索引下载完成，进行下一步"
else
    echo -e "\n错误：索引下载失败，请检查输入的索引名称是否正确！"
    exit 1
fi
#将所有索引文件转移至项目文件夹中：下载完成后返回信息中标注有文件的存储位置。
mkdir humangeome_index
while true; do
    # 提示用户输入绝对路径
    echo -n "请输入下载文件所在绝对路径："
    read -r file_path
    # 检查路径是否存在
    if [ -d "$file_path" ]; then
        # 路径存在，移动路径内所有文件到当前文件夹
        echo -e "\n路径有效，正在移动文件至当前文件夹"
        mv -f "$file_path"/* ./humangeome_index 2>/dev/null
        
        # 简单判断移动是否成功（检查是否有文件被移动）
        if [ $? -eq 0 ] || [ -n "$(ls -A "$file_path" 2>/dev/null)" ]; then
            echo -e "文件移动完成。"
            break
        else
            echo -e "路径内无有效文件，请重新输入。"
        fi
    else
        # 路径不存在，提示错误并重新输入
        echo -e "路径无效，请重新输入。\n"
    fi
done

#选择2（适用高性能服务器）
#直接下载人类基因组数据 若下载失败，可用网址：https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
#~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz ./
#解压文件
#gunzip human_g1k_v37.fasta.gz
#建立索引文件
#bwa-mem2 index human_g1k_v37.fasta


#序列比对：样本和人类基因组之间比对
bwa mem -t 16 ./humangeome_index/*.fa clean.1.fq.gz clean.2.fq.gz >ncov.sam
#正反链均比对不上的reads保存至ncov.1.fq和ncov.2.fq并压缩
samtools fastq -G 2 ncov.sam -1 ncov.1.fq -2 ncov.2.fq
gzip ncov.1.fq ncov.2.fq
#去除人类基因组前后的数据统计比较
seqkit stat clean.1.fq.gz clean.2.fq.gz >humanfilter_stat
seqkit stat ncov.1.fq.gz ncov.2.fq.gz >>humanfilter_stat
cat humanfilter_stat

# 询问用户输入线程数，默认值10
echo -n "即将开始基因组拼接"
echo -n "请输入线程数（-t 参数，默认为10）："
read -r thread_num
# 用户未输入时，默认使用10
thread_num=${thread_num:-10}

# 显示当前选择的线程数（可选，提升体验）
echo -e "\n已选择线程数：$thread_num，开始拼接。"

# 执行 megahit 命令
megahit -t "$thread_num" -o megahit/ -1 ncov.1.fq -2 ncov.2.fq 1>megahit.log 2>megahit.err
# 简单判断命令执行结果（可选）
if [ $? -eq 0 ]; then
    echo -e "拼接完成，日志文件：megahit.log，错误日志：megahit.err"
else
    echo -e "拼接失败,！请查看错误日志：megahit.err"
    exit 1
fi

#查看拼接结果
cd ./megahit
seqkit stat final.contigs.fa >final.contigs_stat
seqkit fx2tab final.contigs.fa  |awk -F"\t" '{print $1,length($2)}' >>final.contigs_stat
cat final.contigs_stat
 
#下载新冠病毒数据
# 下载目标目录（不存在则创建）
DOWNLOAD_DIR="./database"
[ ! -d "$DOWNLOAD_DIR" ] && echo "创建下载目录: $DOWNLOAD_DIR" && mkdir -p "$DOWNLOAD_DIR"
# ascp 基础命令（固定参数部分）
ASCP_BASE_CMD="~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l 6000m"
# 远程文件前缀（NCBI Betacoronavirus数据库路径）
REMOTE_PREFIX="anonftp@ftp.ncbi.nlm.nih.gov:blast/db/Betacoronavirus."

#校验md5完整性并显示结果
echo -e "MD5完整性校验"
md5sum -c *.md5 > md5_test.log 2>&1
cat md5_test.log
echo -e "校验完成，请查看出错文件"

# 补充出错的下载文件
while true; do
    # 读取用户输入（允许空输入，提示更清晰）
    read -p $'\n请输入错误/缺失的文件编号（如 00 01 02），直接回车结束下载：' INPUT

    # 处理输入：去除前后空格，判断是否为空
    TRIMMED_INPUT=$(echo "$INPUT" | xargs)  # xargs自动去除前后空格和空行
    if [ -z "$TRIMMED_INPUT" ]; then
        echo -e "文件完整，即将进行比对流程。"
        break
    fi

    # 构建输入的编号数组
    NUMBERS=($TRIMMED_INPUT)  # 按空格分割为数组
    echo -e "\n补充下载文件编号：${NUMBERS[*]}"

    # 构造完整的ascp下载命令
    ASCP_FULL_CMD="$ASCP_BASE_CMD"
    for NUM in "${NUMBERS[@]}"; do
        ASCP_FULL_CMD+=" $REMOTE_PREFIX$NUM.*"
    done
    ASCP_FULL_CMD+=" $DOWNLOAD_DIR"

    #下载数据
    echo -e "开始下载：$ASCP_FULL_CMD"
    $ASCP_FULL_CMD

    # 6. 下载完成后，重新执行MD5校验
    echo -e "\n下载完成，重新进行MD5校验"
    md5sum -c *.md5 > md5_test.log 2>&1
    cat md5_test.log
    echo -e "重新校验完成"
done

#解压下载文件后删除压缩包和md5文件
cd ./database
for i in *.tar.gz
do
    tar -zxvf "$i"
done
echo -n "数据库解压完成"
#rm -rf ./database/*.tar.gz *.md5
#echo -n "下载完成，压缩包已删除"
cd ../ #回到megahnit文件

#将拼接后的序列与新冠病毒数据库比对
# 定义默认参数值
default_evalue="1e-5"
default_num_threads="10"
# 询问用户evalue值
read -p "请输入evalue值（默认：$default_evalue）：" user_evalue
# 若用户未输入，使用默认值
evalue=${user_evalue:-$default_evalue}
# 询问用户num_threads值
read -p "请输入线程数num_threads（默认：$default_num_threads）：" user_threads
# 若用户未输入，使用默认值
num_threads=${user_threads:-$default_num_threads}
# 执行blastn命令（保留原命令所有其他参数，仅替换用户指定的两个值）
echo "开始执行blastn命令，参数：evalue=$evalue，num_threads=$num_threads"
blastn -query final.contigs.fa \
       -db ./database \
       -out blast.out \
       -outfmt 6 \
       -evalue "$evalue" \
       -num_threads "$num_threads" \
       > blast.log

echo "比对完成，输出日志已保存到 blast.log"
#将比对结果保存到桌面
sz blast.log

#情况1（默认）：保存一个目标文件
while true; do
    # 提示用户输入序列ID
    echo -n "请输入序列ID号："
    read seq_id
    # 处理空输入（用户直接回车的情况）
    if [ -z "$seq_id" ]; then
        echo "错误：ID号不能为空，请重新输入！"
        continue
    fi
    # 验证ID是否存在于文件中
    seqkit grep -q -p "$seq_id" final.contigs.fa
    # 根据seqkit的退出状态码判断ID是否有效
    if [ $? -eq 0 ]; then
        # ID有效，提取序列
        echo "提取序列 $seq_id ..."
        seqkit grep -p "$seq_id" final.contigs.fa > "${seq_id}.fa"
        break
    else
        # ID无效，提示错误并重新循环
        echo "ID号错误：$seq_id 不存在于 final.contigs.fa 中，请重新输入！"
    fi
done
#查看文件内容
echo -n "是否查看生成的文件内容？(y/n，默认n)："
read view_choice
# 若用户输入 y/Y，执行 cat 查看文件
if [[ "$view_choice" == "y" || "$view_choice" == "Y" ]]; then
    echo -e "\n===== 文件 $output_file 内容 ====="
    cat "${seq_id}.fa"
    echo -e "===== 查看结束 ====="
fi

#情况2：需要输入多个文件则可去除if false
if false; then
    while true; do
        # 提示用户输入ID（多个用空格分隔）
        read -p "请输入序列ID号：" ids

        # 拆分ID为数组（空格分隔）
        id_array=($ids)

        # 标记是否所有ID都有效
        all_valid=1

        # 验证每个ID是否存在
        for id in "${id_array[@]}"; do
            # 用seqkit检查ID是否存在（--quiet只返回状态码，不输出内容）
            if ! seqkit grep -p "$id" --quiet "$INPUT_FA"; then
                echo "ID号错误：$id 不存在"
                all_valid=0  # 标记存在无效ID
            fi
        done

        # 所有ID都有效则退出循环，否则重新输入
        [ $all_valid -eq 1 ] && break || echo "请重新输入！"
    done
    # 逐个执行命令（每个ID生成对应.fa文件）
    for id in "${id_array[@]}"; do
        seqkit grep -p "$id" "$INPUT_FA" > "${id}.fa"
    done

    echo "所有ID处理完成！输出文件：${id_array[@]/%/.fa}"
fi
