import geneCiteIndexFunction as gcifunction
import time
from datetime import datetime
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np

# ----------------------------------------------------- step 1: -------------------------------------------------------
x = ord('A')
print(x)
# v1.0
# input : S_614
# output: S_D614G
input = 'S_455' # 'S_614'
gene_name, gene_cite = gcifunction.geneSiteInputParse2(input)
index_ = gcifunction.geneSiteToIndex(gene_name,gene_cite) 
index_ = index_ - 1 

print('index:',index_)
print('index of shap:',index_+1)


referenceFile = './data_input/refAA.fasta'

refFasta = SeqIO.parse(referenceFile, "fasta")
for seq in refFasta:
	idRef = seq.id
	seqRef = seq.seq
refAA = seqRef[index_]


loadFile = './data_input/AASeqInput.fasta'
seq_dict = {rec.id: rec.seq for rec in SeqIO.parse(loadFile, "fasta")}
idToAA = dict()
for key in seq_dict.keys():
	seq = seq_dict[key]
	idToAA[key] = seq[index_]

aaTypedict = {}
for key in idToAA:
	aaTypedict[idToAA[key]] = aaTypedict.get(idToAA[key], 0) + 1
print('\n'+'****************aa Type distribution******************')
print(aaTypedict)
print('aa type number: ', len(aaTypedict))

target = ' '
t_num_max = 0
for key in aaTypedict:
	if key != refAA:
		num_tmp = aaTypedict[key]
		if num_tmp > t_num_max:
			target = key
			t_num_max = num_tmp
output = gene_name + '_' + refAA + str(gene_cite) + target
output2 = gene_name + ':' + refAA + str(gene_cite) + target
print('mutation type: ',output)
print('refAA:',refAA,ord(refAA))
print('mutAA:',target,ord(target))

print('\nAA char to number:')
input = 'R'
print(input,ord(input))
input2 = '-'
print(input2,ord(input2))

print('\nAA number to char:')
input_num = 80
print(input_num,chr(input_num))

# ----------------------------------------------------- step 2: -------------------------------------------------------
if 0:
	# Statistics of variants in nextclade
	loadpath = 'meta.tsv'
	tsvList = []
	listColumn = []
	# seqName 0
	# clade   1
	# aaSubstitutions 30
	idToClade = dict()
	idToAAMutation = dict()
	idToDeletion = dict()
	with open(loadpath, 'r', encoding='utf-8') as f:
		firstLine = True
		for line in f:
			line = line.strip('\n').split('\t')
			if firstLine:
				listColumn = line
				firstLine = False
			else:
				tsvList.append(line)
				idToClade[line[0]] = line[1]
				idToAAMutation[line[0]] = line[30]
				idToDeletion[line[0]] = line[31]

	# substitution
	# substitutionName = 'S:N501Y'
	substitutionName = output2
	substitutionName = 'N:R32-'
	substitutionToClade = dict()
	for id in idToAAMutation:
		AAMutations = idToAAMutation[id]
		if AAMutations.find(substitutionName) != -1:
			substitutionToClade[id] = idToClade[id]

	Cladedict = {}
	for key in substitutionToClade:
		Cladedict[substitutionToClade[key]] = Cladedict.get(substitutionToClade[key], 0) + 1
	print('\n'+'****************Clade Type distribution******************')
	print('substitutionName:',substitutionName)
	print(Cladedict)
	print('Clade type number: ', len(Cladedict))

	# deletion
	delectionName = 'N:G30-'
	delectionToClade = dict()
	for id in idToDeletion:
		AADelections = idToDeletion[id]
		if AADelections.find(delectionName) != -1:
			delectionToClade[id] = idToClade[id]

	Cladedict2 = {}
	for key in delectionToClade:
		Cladedict2[delectionToClade[key]] = Cladedict2.get(delectionToClade[key], 0) + 1
	print('\n'+'****************Clade Type distribution del******************')
	print('delectionName:',delectionName)
	print(Cladedict2)
	print('del Clade type number: ', len(Cladedict2))

