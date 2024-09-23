import glob
import os
import pandas as pd
from Bio import SeqIO
import csv
import numpy


file = []
dir = []
dir_res = []
readPath = './data_input/fasta/'
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
print("file：", file)
# print("dir：", dir)

rewriteFastaSampleNumber = 1000000


idList = []
sequencesList = []
# tsvList = []
saveSeqId = 0

for file1 in file:
	fastaPath = file1
	# tsvPath = readPath + str(fileIndex) + '.tsv'
	inputFasta = SeqIO.parse(fastaPath, "fasta")

	for seq in inputFasta:
		idList.append(seq.id)
		sequencesList.append(seq.seq)
		saveSeqId += 1
		if saveSeqId%10000==0:
			print('sequence number: ',saveSeqId)


print('loading fasta files finished!')

rewriteFastaIndex = 1
rewriteFastaFileName = readPath + str(rewriteFastaIndex) + '.fasta'
rewriteId = 0
fa_out = open(rewriteFastaFileName, 'w', newline='')
for seq in sequencesList:
    fa_out.write(">" + idList[rewriteId] + "\n")
    fa_out.write(str(seq) + "\n")
    rewriteId += 1
    if rewriteId%rewriteFastaSampleNumber == 0:
        print('rewriteFastaNumber: ',rewriteFastaIndex)
        fa_out.close()
        rewriteFastaIndex += 1
        rewriteFastaFileName = readPath + str(rewriteFastaIndex) + '.fasta'
        fa_out = open(rewriteFastaFileName, 'w', newline='')


print('done!')
