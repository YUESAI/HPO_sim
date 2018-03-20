import networkx as nx
import numpy as np
import sys

def main():

	similarity = {}
	simi = []
	with open(r'hpo_similarity.txt','r') as f:
		for line in f:
			line_data = line.strip().split('\t')
			simi.append(float(line_data[2]))
			similarity[(int(line_data[0]), int(line_data[1]))] = float(line_data[2])

	sys.path.append("tool")
	import plotter
	plotter.plot(similarity, 'data/GO/GO_CC.txt', 0.1)


if __name__ == '__main__':
	main()


