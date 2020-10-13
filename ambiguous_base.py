# _*_ coding: utf-8 _*_
import pandas as pd
import numpy as np
import os
import sys

print(
'''
====================================================================================
██╗   ██╗██╗   ██╗███████╗██╗   ██╗
╚██╗ ██╔╝██║   ██║██╔════╝╚██╗ ██╔╝
 ╚████╔╝ ██║   ██║███████╗ ╚████╔╝ 
  ╚██╔╝  ██║   ██║╚════██║  ╚██╔╝  
   ██║   ╚██████╔╝███████║   ██║   
   ╚═╝    ╚═════╝ ╚══════╝   ╚═╝   
                                   
Time: 2020年10月2日
Author: Yusy
Version: V0.1
File: ambiguous_base.py
Describe: The only requirement is that the input file should be a aligined FASTA file
Usage: python ambiguous_base.py
=====================================================================================
''')

if __name__ == '__main__':
    # 模糊碱基列表
    ambiguous_bases = {
        'R': ['A', 'G'],
        'Y': ['C', 'T', 'U'],
        'M': ['A', 'C'],
        'K': ['G', 'T'],
        'S': ['C', 'G'],
        'W': ['A', 'T'],
        'H': ['A', 'C', 'T'],
        'B': ['C', 'G', 'T'],
        'V': ['A', 'C', 'G'],
        'D': ['A', 'G', 'T'],
        'N': ['A', 'C', 'G', 'T']
    }

    #filename = input('请输入需要处理的文件名：')
    filename = sys.argv[1]
    names = []
    samples = []

    with open(filename,'r') as fr:
        seq_list = fr.read().split(">")
        for seq in seq_list:
            if seq:
                seq_name = seq.split("\n")[0]
                seq_fa = "".join(seq.split("\n")[1:])
                seq_fa = seq_fa.upper()
                name = ">" + seq_name
                names.append(name)
                samples.append(seq_fa)
    reads = samples

    for read in reads:
        for i in range(len(read)):
            with open('base.txt', 'a') as f:
                f.write(str(ord(read[i])))
                f.write('\t')
        with open('base.txt','a') as f:
            f.write('\n')

    base_table = pd.read_table("base.txt",sep='\t',header=None)
    os.remove('base.txt')
    base_table = base_table.loc[:,list(base_table.columns)[:-1]]#去除最后一列多出的tab
    #print(base_table)
    #计算每一列的平均数
    means = []
    for col in base_table.columns:
        means.append(base_table[col].mean())
    #print("每一列的平均数",means)

    #匹配模糊碱基列表
    ambiguous_base_ascii = []
    for ambiguous_base in ambiguous_bases.keys():
        ambiguous_base_ascii.append(ord(ambiguous_base))
    #print("模糊碱基的ascii：",ambiguous_base_ascii)

    # 数据框->列表
    data_array = np.array(base_table)
    data_list =data_array.tolist()
    #print(len(data_list[0]))

    # 寻找模糊碱基，输出所有模糊碱基和其所在的列
    # 三个/四个可选择序列选择的做法有点冗余,根据二分法可以精简,改了一个不想改了
    all_ambiguous_bases = []
    all_ambiguous_bases_cols = []
    for i in range(len(data_list)):
        for j in range(len(data_list[0])):
            if data_list[i][j] in ambiguous_base_ascii:
                # 所有序列中匹配到的模糊碱基
                all_ambiguous_bases.append(data_list[i][j])
                all_ambiguous_bases_cols.append(j)
                # 确定填充哪个
                for k in range(len(all_ambiguous_bases)):
                    if all_ambiguous_bases[k] == 82:
                        if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('G') - means[all_ambiguous_bases_cols[k]]):
                            data_list[i][j] = ord('A')#print('此处应该改为A')
                        else:
                            data_list[i][j] = ord('G')#print('此处应该改成G')
                    elif all_ambiguous_bases[k] == 89:
                        if abs(ord('C') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                            if abs(ord('C') - means[all_ambiguous_bases_cols[k]]) < abs(ord('U') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('C')#print('此处应该改成C')
                            else:
                                data_list[i][j] = ord('U')#print('此处应该改成U')
                        else:
                            if abs(ord('T') - means[all_ambiguous_bases_cols[k]]) < abs(ord('U') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('T')#print('此处应该改成T')
                            else:
                                data_list[i][j] = ord('U')#print('此处应该改成U')
                    elif all_ambiguous_bases[k] == 77:
                        if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('C') - means[all_ambiguous_bases_cols[k]]):
                            data_list[i][j] = ord('A')#print('此处应该改为A')
                        else:
                            data_list[i][j] = ord('C')#print('此处应该改成C')
                    elif all_ambiguous_bases[k] == 75:
                        if abs(ord('G') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                            data_list[i][j] = ord('G')#print('此处应该改为G')
                        else:
                            data_list[i][j] = ord('T')#print('此处应该改成T')
                    elif all_ambiguous_bases[k] == 83:
                        if abs(ord('C') - means[all_ambiguous_bases_cols[k]]) < abs(ord('G') - means[all_ambiguous_bases_cols[k]]):
                            data_list[i][j] = ord('C')#print('此处应该改为C')
                        else:
                            data_list[i][j] = ord('G')#print('此处应该改成G')
                    elif all_ambiguous_bases[k] == 87:
                        if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                            data_list[i][j] = ord('A')#print('此处应该改为A')
                        else:
                            data_list[i][j] = ord('T')#print('此处应该改成T')
                    elif all_ambiguous_bases[k] == 72:
                        if abs(ord('C') - means[all_ambiguous_bases_cols[k]]) < abs(ord('A') - means[all_ambiguous_bases_cols[k]]):
                            if abs(ord('C') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('C')#print('此处应该改成C')
                            else:
                                data_list[i][j] = ord('T')#print('此处应该改成T')
                        else:
                            if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('A')#print('此处应该改成A')
                            else:
                                data_list[i][j] = ord('T')#print('此处应该改成T')
                    elif all_ambiguous_bases[k] == 66:
                        if abs(ord('C') - means[all_ambiguous_bases_cols[k]]) < abs(ord('G') - means[all_ambiguous_bases_cols[k]]):
                            # print("data_list{}{}".format(i,j))
                            data_list[i][j] = ord('C')
                            # print('此处应该改成C')
                        else:
                            if abs(ord('G') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('G')#print('此处应该改成G')
                            else:
                                data_list[i][j] = ord('T')#print('此处应该改成T')
                    elif all_ambiguous_bases[k] == 86:
                        if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('C') - means[all_ambiguous_bases_cols[k]]):
                            if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('A')#print('此处应该改成A')
                            else:
                                data_list[i][j] = ord('T')#print('此处应该改成T')
                        else:
                            if abs(ord('C') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('C')#print('此处应该改成C')
                            else:
                                data_list[i][j] = ord('T')#print('此处应该改成T')
                    elif all_ambiguous_bases[k] == 68:
                        if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('G') - means[all_ambiguous_bases_cols[k]]):
                            if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('A')#print('此处应该改成A')
                            else:
                                data_list[i][j] = ord('T')#print('此处应该改成T')
                        else:
                            if abs(ord('G') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                data_list[i][j] = ord('G')#print('此处应该改成G')
                            else:
                                data_list[i][j] = ord('T')#print('此处应该改成T')
                    elif all_ambiguous_bases[k] == 78:
                        if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('G') - means[all_ambiguous_bases_cols[k]]):
                            if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                if abs(ord('A') - means[all_ambiguous_bases_cols[k]]) < abs(ord('C') - means[all_ambiguous_bases_cols[k]]):
                                    data_list[i][j] = ord('A')#print('此处应该改成A')
                                else:
                                    data_list[i][j] = ord('C')#print('此处应该改成C')
                            else:
                                if abs(ord('T') - means[all_ambiguous_bases_cols[k]]) < abs(ord('C') - means[all_ambiguous_bases_cols[k]]):
                                    data_list[i][j] = ord('T')#print("此处应该改成T")
                                else:
                                    data_list[i][j] = ord('C')#print("此处应该改成C")
                        else:
                            if abs(ord('G') - means[all_ambiguous_bases_cols[k]]) < abs(ord('T') - means[all_ambiguous_bases_cols[k]]):
                                if abs(ord('G') - means[all_ambiguous_bases_cols[k]]) < abs(ord('C') - means[all_ambiguous_bases_cols[k]]):
                                    data_list[i][j] = ord('G')#print("此处应该改成G")
                                else:
                                    if abs(ord('T') - means[all_ambiguous_bases_cols[k]]) < abs(ord('C') - means[all_ambiguous_bases_cols[k]]):
                                        data_list[i][j] = ord('T')#print("此处应该改成T")
                                    else:
                                        data_list[i][j] = ord('C')#print("此处应该改成C")

    for i in range(len(data_list)):
        for j in range(len(data_list[0])):
            data_list[i][j] = chr(data_list[i][j])
    for i in range(len(data_list)):
        with open('clean_data.fasta','a') as file:
            file.write(names[i])
            file.write('\n')
            for j in range(len(data_list[i])):
                file.write(data_list[i][j])
            file.write('\n')



















