import glob
import os
import pandas as pd
from Bio import SeqIO
import csv
import numpy

reffastaPath = './data/refAA.fasta'
refFasta = SeqIO.parse(reffastaPath, "fasta")
for seq in refFasta:
    idRef = seq.id
    seqRef = seq.seq


sequence_file = './data_input/Input_AAseq.fasta'

# ID2Clade = dict()
#
# with open(lineage_file, 'r') as f:
# 	next(f)
# 	for line in f:
# 		cladeTypeIndex = 17
# 		line = line.strip('\n').split('\t')
# 		ID2Clade[line[0]] = line[cladeTypeIndex]

seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(sequence_file, "fasta")}

mutation_frequency_res_file = open('./result/mutation_rate.txt', 'w')

feature_site = 0
for site in range(0,len(seqRef)):
	diff_num = 0
	total_num = 0
	for key in seq_dict.keys():
		queryline = seq_dict[key]
		querysite = queryline[site].upper()
		refsite   = seqRef[site].upper()
		total_num += 1
		if querysite != refsite and refsite != '*' and querysite != '*':
			diff_num += 1
	mutation_frequency_res_file.write(str(site) + ' ' + str(diff_num / total_num) + ' ' + str(refsite) + '\n')
	if diff_num/total_num >= 0.0001:
		feature_site += 1
		print(site,' ',refsite,' ', diff_num/total_num)

mutation_frequency_res_file.close()
print('feature_site: ',feature_site)
print('done')
