#!/bin/bash
#set -x  # 启用调试
output_prefix="./outputdatas/target_data"  # 修改保存路径
# Specify the CSV file and model.props file
csv_file="samestep1.csv"
props_file="props/model.props"
line_count=$(wc -l < samestep1.csv)
echo "Number of lines in the CSV file: $line_count"

# 初始化计数器
counter=1
# Process each line of the CSV file
#i=1
#while [ "$i" -lt 5 ]; do
while IFS=',' read -r corn1 corn2 corn3 corn4 corn5; do

                                                                           
    # 检查 model.props 是否有写入权限
    # if [ ! -w "$props_file" ]; then
    #     echo "Error: No write permission for '$props_file'."
    #     exit 1
    # fi                                                                      
    # Update the model.props file

     echo "Updating model.props with values: $corn1, $corn2, $corn3, $corn4, $corn5"

     sed -i "s/^\(max\.fission\.age\s*=\s*\).*\$/\1$corn1/" $props_file
     #sed -i "s/^\(min\.fission\.age\s*=\s*\).*\$/\1$corn1/" $props_file
    
     sed -i "s/^\(max\.death\.age\s*=\s*\).*\$/\1$corn2/" $props_file
     #sed -i "s/^\(min\.death\.age\s*=\s*\).*\$/\1$corn2/" $props_file

     sed -i "s/^\(annual\.variance\s*=\s*\).*\$/\1$corn3/" $props_file

     sed -i "s/^\(fertility\.prop\s*=\s*\).*\$/\1$corn4/" $props_file

     sed -i "s/^\(harvest\.adj\s*=\s*\).*\$/\1$corn5/" $props_file

    
    echo "准备运行 mpirun 命令..."
    mpirun -n 1 bin/main.exe props/config.props props/model.props < /dev/null
    #exit_status=$?
    #echo "mpirun 命令${counter}执行完毕，退出状态：$exit_status"



    # 重命名输出文件为 target_data+counter.csv，并移动到新路径
    output_file="${output_prefix}${counter}.csv"
    cp "NumberOfHousehold.csv" "$output_file"

    echo "计数器现在的值: $counter"
    echo "采样矩阵第${counter}行处理完成。"
    # 等待 mpirun 进程完成
    ((counter++))
    #i=$((i+1))
done < "$csv_file"