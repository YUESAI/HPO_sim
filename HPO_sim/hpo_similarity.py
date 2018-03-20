import networkx as nx
import sys
import numpy as np
import os

# 数据集路径
DIR = 'data/HPO_GENESET'

# =============================================================================

def read_gene_list(gene_file):
	genes_set = set()
	# 获取文件中所有基因
	for line in open(DIR + '\\' + gene_file, 'r'):
		if line[0] == '#':
			continue
		genes_set.add(int(line.strip().split('\t')[0]))
	return genes_set

# =============================================================================

def get_mods(T, gene_set):
	mods = set()
	cnt = 0
	for i in gene_set:
		try:
			mods.update(T[i])
		except:
			mods.add(-1)
			cnt = cnt + 1
	# 对超过50%的基因不在字典中的情况拒判
	if cnt >= len(gene_set) * 0.5:
		raise '未能识别集合'
	return mods

def jac_sim(mods_1, mods_2):
	# Jaccard Index
	return len(mods_1 & mods_2) / len(mods_1 | mods_2)

def BP_sim(mods_1, mods_2):
	# 交集的大小
	l_and = len(mods_1 & mods_2)
	# Bin Prediction Accuracy
	return (l_and / len(mods_1) +  l_and / len(mods_2)) / 2

def aff_sim(mods_1, mods_2):
	# Affinity Score
	return np.square(len(mods_1 & mods_2)) / len(mods_1) / len(mods_2)

# =============================================================================
def main():
	# 处理网络
	sys.path.append("tool")
	import plotter
	import FM

	# 读取功能模块
	T = FM.mod("data/PPI/dict.txt").belones
	genes = set(T.keys())

	# 读取测试集
	pheno_list = list()
	hpo_mods = dict()
	for i in os.listdir(DIR):
		genes_set = read_gene_list(i)
		pid = int(i.split('.txt')[0])
		# 只保留网络中存在的基因
		try:
			hpo_mods[pid] = get_mods(T, genes_set)
		except:
			continue
		pheno_list.append(pid)
	
	print("done reading files")
		
	# 计算相似度

	similarity = dict()
	mk = list()

	# 对测试集中数据两两比较
	for i in pheno_list:
		mk.append(i)
		mods_A = hpo_mods[i]
		for j in pheno_list:
			if j in mk:
				continue
			mods_B = hpo_mods[j]
			d_AB = jac_sim(mods_A, mods_B)
			#d_AB = BP_sim(mods_A, mods_B)
			#d_AB = aff_sim(mods_A, mods_B)
			similarity[i, j] = d_AB
		#print('done',i)
	print("done processing")
	# 输出结果
	plotter.plot(similarity, 'data/GO/GO_BP.txt', 0.05)
# =============================================================================
if __name__ == '__main__':
	main()



