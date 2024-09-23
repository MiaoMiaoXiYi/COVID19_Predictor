import sys
import Bio
import os
import time
from Bio import GenBank
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
import json
from gencoder import translate_dna
from Bio.Seq import MutableSeq


with open('SARS-CoV-2.json', 'r') as fr:
	json_file = json.loads(fr.read())
genome = json_file['genome']
genes = json_file['genes']

ORF1a = genes['ORF1a']
ORF1a_start = ORF1a['coordinates']['from']
ORF1a_to = ORF1a['coordinates']['to']

ORF1b = genes['ORF1b']
ORF1b_start = ORF1b['coordinates']['from']
ORF1b_to = ORF1b['coordinates']['to']

S = genes['S']
S_start = S['coordinates']['from']  
S_to = S['coordinates']['to']

ORF3a = genes['ORF3a']
ORF3a_start = ORF3a['coordinates']['from']
ORF3a_to = ORF3a['coordinates']['to']

E = genes['E']
E_start = E['coordinates']['from']
E_to = E['coordinates']['to']

M = genes['M']
M_start = M['coordinates']['from']
M_to = M['coordinates']['to']

ORF6 = genes['ORF6']
ORF6_start = ORF6['coordinates']['from']
ORF6_to = ORF6['coordinates']['to']

ORF7a = genes['ORF7a']
ORF7a_start = ORF7a['coordinates']['from']
ORF7a_to = ORF7a['coordinates']['to']

ORF7b = genes['ORF7b']
ORF7b_start = ORF7b['coordinates']['from']
ORF7b_to = ORF7b['coordinates']['to']

ORF8 = genes['ORF8']
ORF8_start = ORF8['coordinates']['from']
ORF8_to = ORF8['coordinates']['to']

N = genes['N']
N_start = N['coordinates']['from']  # N:28274-29533
N_to = N['coordinates']['to']

ORF10 = genes['ORF10']
ORF10_start = ORF10['coordinates']['from']
ORF10_to = ORF10['coordinates']['to']


lenthAccum = 0
transGenesStart = dict()
transGenesStart['ORF1a'] = 1
lenthAccum += (genes['ORF1a']['coordinates']['to'] - genes['ORF1a']['coordinates']['from'] + 1) / 3
transGenesStart['ORF1b'] = lenthAccum + 1
lenthAccum += (genes['ORF1b']['coordinates']['to'] - genes['ORF1b']['coordinates']['from'] + 1) / 3
transGenesStart['S'] = lenthAccum + 1
lenthAccum += (genes['S']['coordinates']['to'] - genes['S']['coordinates']['from'] + 1) / 3
transGenesStart['ORF3a'] = lenthAccum + 1
lenthAccum += (genes['ORF3a']['coordinates']['to'] - genes['ORF3a']['coordinates']['from'] + 1) / 3
transGenesStart['E'] = lenthAccum + 1
lenthAccum += (genes['E']['coordinates']['to'] - genes['E']['coordinates']['from'] + 1) / 3
transGenesStart['M'] = lenthAccum + 1
lenthAccum += (genes['M']['coordinates']['to'] - genes['M']['coordinates']['from'] + 1) / 3
transGenesStart['ORF6'] = lenthAccum + 1
lenthAccum += (genes['ORF6']['coordinates']['to'] - genes['ORF6']['coordinates']['from'] + 1) / 3
transGenesStart['ORF7a'] = lenthAccum + 1
lenthAccum += (genes['ORF7a']['coordinates']['to'] - genes['ORF7a']['coordinates']['from'] + 1) / 3
transGenesStart['ORF7b'] = lenthAccum + 1
lenthAccum += (genes['ORF7b']['coordinates']['to'] - genes['ORF7b']['coordinates']['from'] + 1) / 3
transGenesStart['ORF8'] = lenthAccum + 1
lenthAccum += (genes['ORF8']['coordinates']['to'] - genes['ORF8']['coordinates']['from'] + 1) / 3
transGenesStart['N'] = lenthAccum + 1
lenthAccum += (genes['N']['coordinates']['to'] - genes['N']['coordinates']['from'] + 1) / 3
transGenesStart['ORF10'] = lenthAccum + 1
lenthAccum += (genes['ORF10']['coordinates']['to'] - genes['ORF10']['coordinates']['from'] + 1) / 3

target1 = transGenesStart['S'] - 1 + 614  
target1 = transGenesStart['ORF1b'] - 1 + 2612  
# print('target1:',target1)
target2 = transGenesStart['ORF3a'] - 1 + 1

def geneSiteToIndex(gene, site):  
    index_ = transGenesStart[gene] - 1 + site 
    index_ = int(index_)

    return index_


def indexToGeneSite(index_): 
    last_key = 'none'
    last_site = 0
    first_loop = True
    for key in transGenesStart:
        if index_ <= (genes['ORF1a']['coordinates']['to'] - genes['ORF1a']['coordinates']['from'] + 1) / 3:
            last_key = key
            last_site = index_
            break
        elif first_loop:
            last_key = key
            first_loop = False
        else:
            if transGenesStart[key] > index_:
                last_site = index_ + 1 - transGenesStart[last_key]
                break
            else:
                last_key = key
        if key == 'ORF10' and index_ >= transGenesStart[key]:
            last_key = key
            last_site = index_+1-transGenesStart[key]
    geneSiteName = last_key+'_'+str(int(last_site))
    return geneSiteName


# for key in transGenesStart:
#     print(key)
# gene_input = 'S'
# aa_site = 614
def geneSiteInputParse(genesiteInpt):
    # genesiteInpt = 'S:D614'   # S:D614G or S:D614   ('gene'+':'+'refAA'+'Site')
    split_site = genesiteInpt.find(':')
    gene_input = genesiteInpt[0:split_site]
    last_str = genesiteInpt[-1]
    if genesiteInpt[split_site + 1:split_site + 5] != 'STOP':
        if genesiteInpt.find('STOP') != -1:
            gene_site = int(genesiteInpt[split_site + 2:-4])
        else:
            if last_str.isdigit():
                gene_site = int(genesiteInpt[split_site + 2:])
            else:
                gene_site  = int(genesiteInpt[split_site+2:-1])
        return gene_input, gene_site
    else:
        return 'none', -1

def geneSiteInputParse2(genesiteInpt):
    # genesiteInpt = 'ORF1a_610'
    split_site = genesiteInpt.find('_')
    gene_input = genesiteInpt[0:split_site]
    gene_site = int(genesiteInpt[split_site + 1:])
    return gene_input, gene_site



# # test
# genesiteInpt = 'S:D614'
# gene_input, gene_site = geneSiteInputParse(genesiteInpt)
# print('Input: ',gene_input + ':' + str(gene_site))
#
# print('Index (start from 1): ', geneSiteToIndex(gene_input, gene_site))
# print('GeneStieName: ', indexToGeneSite(geneSiteToIndex(gene_input, gene_site)))


