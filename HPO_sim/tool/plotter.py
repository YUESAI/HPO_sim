import numpy as np
import matplotlib.pyplot as plt


def read_GO_similarity(file_name):
	go_similarity = {}
	for line in open(file_name,'r'):
			line_data = line.strip().split('\t')
			go_similarity[(int(line_data[0]), int(line_data[1]))] = float(line_data[2])
	return go_similarity

def correlation(p, min_value , para):
	medium = []
	y = []
	for i in range(len(p)):
		if len(p[i])  == 0:
			continue
		medium.append(float(np.percentile(p[i], 50)))
		y.append(float('%.3f' % (i * para + min_value + 0.5 * para)))
	result = np.corrcoef(medium, y)
	return result

def plot(hpo_similarity, file_name, para):
	min_value = float('%.2f' % min(hpo_similarity.values()))
	max_value = float('%.2f' % max(hpo_similarity.values()))

	length = int(np.ceil((max_value - min_value) / para))

	p = [ [] for i in range(length)]

	go_similarity = read_GO_similarity(file_name)

	for hpo_pair in hpo_similarity.keys():
		similarity = hpo_similarity[hpo_pair]
		for i in range(length):
			if similarity >= (i * para + min_value) and similarity < (i * para + min_value + para):
				p[i].append(go_similarity[hpo_pair])

	x = [float('%.2f' % (i * para + min_value)) for i in range(length)]

	plt.boxplot(p, showcaps = False, showfliers = False, positions = x, widths = 0.4 * para)
	
	plt.xlabel('similarity')
	plt.ylabel('GO_based similarity')
	plt.title('HPO similarity')
	Correlation_coefficient = correlation(p, min_value, para)[0][1]
	print(" the correlation coefficient between HPO_based and GO_based method is %s" % abs(Correlation_coefficient))
	plt.show()
	#return correlation(p, min_value , length )[0][1]









