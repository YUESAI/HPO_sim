import numpy as np
import struct
import sys

# SLPA结果文件路径
GROUP_FILE = '../data/PPI/SLPAw_ppi_run1_r0.01_v3_T100.icpm'

def read_group():
	gdict = []
	for line in open(GROUP_FILE, 'r'):
		if line[0] == '#':
			continue
		line_data = line.strip().split()
		if len(line_data) == 0:
			continue
		gdict.append([int(i) for i in line_data])
	return gdict

# 获取第一个连通子图模块划分
gdict = read_group()

import PPI
# 读取蛋白质网络文件
G = PPI.net("../data/PPI/PPI.info")
# PPI网络中的所有基因
genes = G.genes

# 将其他子图的每个子图作为一个模块加入
PPI_dict = [set() for i in range(0, len(G._Sub))]
for i in genes:
	sub_id = G._Map[G._Dict[i]]
	PPI_dict[sub_id].add(i)

for i in PPI_dict[1 : len(PPI_dict)]:
	gdict.append(i)

# 保存划分情况到文件
l_tax = len(gdict)

DICT_FILE = '../data/PPI/dict.txt'

dict_file = open(DICT_FILE, 'wb')

dict_file.write(struct.pack("i", l_tax))

for i in range(0, l_tax):
	dict_file.write(struct.pack("i", i))

for i in range(0, l_tax):
	dict_file.write(struct.pack("i", len(gdict[i])))

for i in range(0, l_tax):
	for j in gdict[i]:
		dict_file.write(struct.pack("i", j))

