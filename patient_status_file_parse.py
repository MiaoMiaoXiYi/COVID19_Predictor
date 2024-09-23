import glob
import os
import pandas as pd
from Bio import SeqIO
import csv
import numpy



file = []
dir = []
dir_res = []
readPath = './data_input/tsv/'
def list_dir(start_dir):
    dir_res = os.listdir(start_dir)
    for path in dir_res:
        temp_path = start_dir + '/' + path
        if os.path.isfile(temp_path):
            file.append(temp_path)
        if os.path.isdir(temp_path):
            dir.append(temp_path)
            list_dir(temp_path)

list_dir(readPath)
# print("file：", file)
# print("dir：", dir)

# rewriteFastaSampleNumber = 5000000


# idList = []
# sequencesList = []
tsvList = []
saveSeqId = 0

info_size = 16  #status 16 ; tech:14   for 20240114--0910: 17


for file1 in file:
	tsvPath = file1

	with open(tsvPath, 'r', encoding='utf-8') as f:
		next(f)  
		for line in f:
			line = line.strip('\n').split('\t')
			# line = line[:-1]
			if len(line) != info_size:
				continue
			tsvList.append(line)

print('loading tsv files finished!')


rewriteTsvFile = readPath + 'PatientStatus.tsv'
with open(rewriteTsvFile, 'w', newline='') as f:
	tsv_w = csv.writer(f, delimiter='\t')
	if info_size == 16:
		tsv_w.writerow(
			['Virus name', 'Accession ID', 'Collection date', 'Location', 'Host', 'Additional location information', 'Sampling strategy', 'Gender', \
			 'Patient age', 'Patient status', 'Last vaccinated', 'Passage', 'Specimen', 'Additional host information', 'Lineage', \
			 'Clade'])

	tsv_w.writerows(numpy.array(tsvList).tolist()) 

print('done!')
