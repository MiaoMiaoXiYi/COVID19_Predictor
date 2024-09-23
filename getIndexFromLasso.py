import geneCiteIndexFunction as gcifunction
import time
from datetime import datetime
import matplotlib.pyplot as plt


continent = 'All'

fileInput = './result/lasso_output.txt'

inputList = []
parseOutput = []

with open(fileInput, 'r', encoding='utf-8') as f:
	for line in f:
		line = line.strip('\n').split(' ')
		inputList.append(line)

print('loading input file finished:',fileInput)

filename = './indiciesToKeepFromLasso.txt'
fileFeatureSelected = open(filename, 'w')

for i in range(len(inputList)-14):  # the last three features are age gender and clade ...
	input = inputList[i][1]
	coefficient = float(inputList[i][2])
	gene_name, gene_cite = gcifunction.geneSiteInputParse2(input)
	if gene_name in gcifunction.transGenesStart:
		index_ = gcifunction.geneSiteToIndex(gene_name,gene_cite) 
		index_ = index_ - 1
	else:
		index_ = -1
	# print(input,index_)
	if index_ != -1 and coefficient != 0.0:
		fileFeatureSelected.write(str(index_) + '\n')

