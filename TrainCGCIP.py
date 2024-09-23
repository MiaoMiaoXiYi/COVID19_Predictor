import pandas as pd
from pandas.core.ops import invalid
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, \
	confusion_matrix, make_scorer
from datetime import datetime
import joblib
import sys
import os
# os.environ['CUDA_VISIBLE_DEVICES'] = '0'
import time
from sklearn.model_selection import cross_val_score
from Bio import SeqIO
import numpy as np
from sklearn.tree import plot_tree
import matplotlib.pyplot as plt  
from matplotlib import rcParams  

from sklearn.svm import SVC
from sklearn.ensemble import AdaBoostClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import MultinomialNB

import numpy as np
from catboost import CatBoostClassifier, Pool
import matplotlib.pyplot as plt
import pylab as pl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import torch
from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets
from torchvision.transforms import ToTensor, Lambda, Compose
import gc
import seaborn as sns
from sklearn.metrics import confusion_matrix

import xgboost as xgb
from xgboost import XGBClassifier
import csv
from sklearn.model_selection import KFold

from lightgbm.sklearn import LGBMClassifier
import gpboost
import shap

import StatusEncoding as se
import re
from scipy import stats

from scipy.stats import chisquare
from scipy.stats import chi2
from scipy.stats import chi2_contingency
import tsvFilesParse

from sklearn.metrics import roc_auc_score
from hi_lasso.hi_lasso import HiLasso
from random_lasso import random_lasso

from feature_selection import Feature_Selection as fs

from sklearn.feature_selection import chi2
from sklearn.feature_selection import SelectKBest

import geneCiteIndexFunction as gcifunction
from sklearn.model_selection import GridSearchCV

file_ps = './patientStatus.tsv'


target_continent = 'Asia'

shap_flag = False

Feature_Selection_flag = False

OffIndiciesFilePath = "./indiciesToKeep.txt"


k_fold_train = 1  
k_set = 10  

testRF = 0  # better than decision tree
testXGBoost = 0
testLightGBM = 0
testGPBoost = 0




idToContinent = dict()
idToAge = dict()
idToGender = dict()
idToStatus = dict()
idToVacc = dict()
idToReinfected = dict()
idToAddition = dict()
idToLineage = dict()
idToGisaid = dict()
statusOutputs = []
genderOutputs = []
ageOutputs = []

ageToCoding = dict()

# add by rrk 20221026
idToSelected = dict()
idToLabel = dict()

vacc_num = 0
unvacc_num = 0
rein_num = 0

with open(file_ps, 'r') as f:
	next(f)
	for line in f:
		line = line.strip('\n').split('\t')
		id_ = line[0].replace(" ", "")  
		idToGender[id_] = se.gender_encoding(line[7])  
		idToAge[id_] = se.age_encoding(line[8]) 
		idToStatus[id_] = se.status_encoding(line[9], mild_tight_flag) 
		idToVacc[id_] = line[10]
		idToAddition[id_] = line[13]
		idToLineage[id_] = line[14]
		idToGisaid[id_] = line[15]

		idToContinent[id_] = line[3].split('/')[0].replace(" ", "")


		genderOutputs.append(idToGender[id_])
		ageOutputs.append(idToAge[id_])
		ageToCoding[line[8]] = idToAge[id_]

		status_ = idToStatus[id_]
		vacc_ = idToVacc[id_]
		addition_ = idToAddition[id_]

		vaccFromStatus = se.status_encoding_vacc(status_)
		vaccFromVacc = se.vacc_encoding_vacc(vacc_)
		vaccFromAddi = se.addition_encoding_vacc(addition_)

		vaccOutput = 'unknown'
		if vaccFromStatus == 'Vacc' or vaccFromVacc == 'Vacc' or vaccFromAddi == 'Vacc':
			vaccOutput = 'Vacc'
			vacc_num += 1
		elif vaccFromStatus == 'Unvacc' or vaccFromVacc == 'Unvacc' or vaccFromAddi == 'Unvacc':
			vaccOutput = 'Unvacc'
			unvacc_num += 1
		idToVacc[id_] = vaccOutput

		reInfectFromStatus = se.status_encoding_reinfection(status_)
		reInfectFromAddi = se.addition_encoding_reinfection(addition_)

		reInfectOutput = 'unknown'
		if reInfectFromStatus == 'Reinfection' or reInfectFromAddi == 'Reinfection':
			reInfectOutput = 'Reinfection'
			rein_num += 1
		idToReinfected[id_] = reInfectOutput


dataModetest = 1  # 1--n*1     0--n*5     


referenceFile = './data_input/refAA.fasta'  


newSeqId = []

idToPangoType = dict()
idToGisaidType = dict()
# data storage
dataList = []
# dict for lookup efficiency
indiciesToKeep = dict()
indiciesToKeep_ref = dict()

with open(OffIndiciesFilePath, "r") as f2:
	data2 = f2.readlines()
for i in data2:
	indiciesToKeep_ref[int(i)] = True
# add by rrk 20220817
indicies_ref = list(indiciesToKeep_ref.keys())
indicies_ref.sort()

refFasta = SeqIO.parse(referenceFile, "fasta")
for seq in refFasta:
	idRef = seq.id
	seqRef = seq.seq

finalLineRef = []

for index in indicies_ref:
	finalLineRef.append(seqRef[index])

# referenceId = "Wuhan/WH04/2020"
referenceId = "Wuhan/WH01/2019"
referenceSeq = ""

idToLineage = dict()
idToSeq = dict()

mustKeepIds = []
mustKeepLineages = []

tsvListSub = []
idSub = dict()
FastaIdListSub = []
FastaSeqListSub = []


# function for handling weird sequence characters
def clean(x, loc):
	x = x.upper()

	if x == 'T' or x == 'A' or x == 'G' or x == 'C' or x == '-':
		return x

	if x == 'U':
		return 'T'

	# otherwise return value from reference
	return referenceSeq[loc]


def findReferenceSeq():
	with open(referenceFile) as f:
		currentSeq = ""

		for line in f:
			if ">" not in line:
				currentSeq = currentSeq + line.strip()

	f.close()
	return currentSeq


def getDataLine(seqId, seq):
	dataLine = []
	dataLine.append(seqId)

	newSeq = ""

	# for each character in the sequence
	if FormatDataFileFlag == 1:  # add by rrk 20220228
		for index in range(len(seq)):
			newSeq = newSeq + clean(seq[index], index)
		dataLine.append(newSeq)
	else:
		dataLine.append(str(seq))  # add by rrk 20220228

	return dataLine


def readInAndFormatData():
	global nextstrainFailedNum
	# add the data line for the reference seq
	# idToLineage[referenceId] = "A"
	# dataList.append(getDataLine(referenceId, referenceSeq))

	# create a dictionary of sequence ids to their assigned lineages
	# add by rrk 20220220
	# with open(lineage_file, 'r') as f:
	# 	for line in f:
	# 		line = line.strip()
	# 		split = line.split(",")
	# 		idToLineage[split[0]] = split[1]

	# clades = []
	# index_tsv = 0
	# index_ = 0
	# with open(lineage_file, 'r') as f:
	# 	next(f)
	# 	for line in f:
	# 		index_ += 1
	# 		if index_ % test_seq_freq != 0 and loadTrainedModelFlag == 1:
	# 			continue
	# 		line = line.strip('\n').split('\t')
	# 		# add by rrk 20220313
	# 		idtmp = line[0]
	# 		line_lenth = len(line)
	# 		if line_lenth == 2:
	# 			print('line_lenth==2')
	# 			cladeTypeIndex = 1
	# 		else:
	# 			cladeTypeIndex = 18
	#
	# 		clade_input = line[cladeTypeIndex]
	#
	# 		if load_train_ID_from_file == 1 and trainLayered == 1 and fatherLayer == 0:
	# 			if idtmp in ID_train:
	# 				if clade_input != 'None':
	# 					idToLineage[line[0]] = clade_input  # 17: nextstrain types   ;   19  gisaid types
	# 					clades.append(clade_input)
	#
	# 		if load_train_ID_from_file == 0 and trainLayered == 0:
	# 			if clade_input != 'None':
	# 				idToLineage[line[0]] = clade_input  # 17: nextstrain types   ;   19  gisaid types
	# 				clades.append(clade_input)
	#
	# 		if load_train_ID_from_file == 1 and trainLayered == 0:
	# 			if idtmp in ID_train:
	# 				if clade_input != 'None':
	# 					idToLineage[line[0]] = clade_input  # 17: nextstrain types   ;   19  gisaid types
	# 					clades.append(clade_input)
	#
	# cladesdict = {}
	# for key in clades:
	# 	cladesdict[key] = cladesdict.get(key, 0) + 1
	# print('\n' + '****************seq clades distribution******************')
	# print(cladesdict)
	# print('clade type number: ', len(cladesdict))
	# print('***********************************************************' + '\n')
	# # print('seq number been fixed: ',nextstrainFailedNum)

	# close the file
	f.close()

	seq_dict = {rec.id: rec.seq for rec in SeqIO.parse(sequence_file, "fasta")}

	print("files read in, now processing")

	dataListSize = 0
	failed_to_find_lineage = 0
	for key in seq_dict.keys():
		if key in idToSelected and key in idToLabel:  # check if key is in idToLineage, time costing too much
			dataList.append(getDataLine(key, seq_dict[key]))
			dataListSize += 1
		# if dataListSize % 50000 == 0:
		# 	print('dataListSize: ', dataListSize)
		else:
			failed_to_find_lineage += 1
	# print("unable to find the lineage classification for: " + key)
	print('dataListSize: ', len(dataList))
	print('failed_to_find_lineage: ', failed_to_find_lineage)
	# add by rrk 20220315 try
	print('del seq_dict')
	del seq_dict
	gc.collect()

	if FormatDataFileFlag == 1:
		rewriteId = 0
		fa_out = open(FormatDataFile, 'w', newline='')
		for i in range(1, dataListSize):  # add by rrk 20220315 jump the first seq
			if rewriteId % 1000 == 0:
				print('rewrite seq num: ', rewriteId)
			fa_out.write(">" + dataList[i][0] + "\n")
			fa_out.write(str(dataList[i][1]) + "\n")
			rewriteId += 1
		print('rewrite seq num: ', rewriteId)


# find columns in the data list which always have the same value
def findColumnsWithoutSNPs():
	# for each index in the length of each sequence
	for index in range(len(dataList[0][1])):
		keep = False

		# loop through all lines
		for line in dataList:

			# if there is a difference somewhere, then we want to keep it
			if dataList[0][1][index] != line[1][index] or index == 0:
				keep = True
				break

		# otherwise, save it
		if keep:
			indiciesToKeep[index] = True


# remove columns from the data list which don't have any SNPs. We do this because
# these columns won't be relevant for a logistic regression which is trying to use
# differences between sequences to assign lineages
def removeOtherIndices(indiciesToKeep):
	# instantiate the final list
	finalList = []

	indicies = list(indiciesToKeep.keys())
	indicies.sort()

	# while the dataList isn't empty
	while len(dataList) > 0:

		# pop the first line
		line = dataList.pop(0)
		seqId = line.pop(0)

		line = line[0]
		# initialize the finalLine
		finalLine = []

		for index in indicies:
			if index == 0:
				# if its the first index, then that's the lineage assignment, so keep it
				finalLine.append(seqId)
			else:
				# otherwise keep everything at the indices in indiciesToKeep
				finalLine.append(line[index])

		# save the finalLine to the finalList
		finalList.append(finalLine)

	# return
	return finalList


def allEqual(list):
	entries = dict()

	for i in list:
		if i not in entries:
			entries[i] = True

	return len(entries) == 1


def removeAmbiguous():
	idsToRemove = set()
	lineMap = dict()
	idMap = dict()

	for line in dataList:
		keyString = ",".join(line[1:])

		if keyString not in lineMap:
			lineMap[keyString] = []
			idMap[keyString] = []

		# print('debug output')
		# print(keyString)
		# print(line[0])
		# print(idToLineage[line[0]])

		lineMap[keyString].append(idToLabel[line[0]])
		idMap[keyString].append(line[0])

	for key in lineMap:
		if not allEqual(lineMap[key]):

			skipRest = False

			# see if any protected lineages are contained in the set, if so keep those ids
			for lineage in lineMap[key]:
				if lineage in mustKeepLineages:
					skipRest = True

					for i in idMap[key]:
						if lineage != idToLabel[i] and i not in mustKeepIds:
							idsToRemove.add(i)

			# none of the lineages are protected, fire at will
			if not skipRest:

				lineageToCounts = dict()

				aLineage = False
				# find most common lineage
				for lineage in lineMap[key]:
					if lineage not in lineageToCounts:
						lineageToCounts[lineage] = 0

					lineageToCounts[lineage] = lineageToCounts[lineage] + 1
					aLineage = lineage

				m = aLineage
				for lineage in lineageToCounts:
					if lineageToCounts[lineage] > lineageToCounts[m]:
						m = lineage

				for i in idMap[key]:
					if m != idToLabel[i]:
						# print(m)
						idsToRemove.add(i)

	newList = []

	print("keeping indicies:")

	for line in dataList:
		seqIdtmp = line[0]  # add by rrk 20220313
		if removeAmbiguousFlag == 1:
			if line[0] not in idsToRemove:
				# print(line[0])
				line[0] = idToLabel[line[0]]
				if 0:  # line[0] == 'A':
					print('*******************Warning: remove the reference type!**********************')
				else:
					newList.append(line)
					newSeqId.append(seqIdtmp)
		else:
			if 1:  # line[0] not in idsToRemove:
				# print(line[0])
				line[0] = idToLabel[line[0]]
				if 0:  # line[0] == 'A':
					print('*******************Warning: remove the reference type!**********************')
				else:
					age = idToAge[seqIdtmp]

					if age_valid:
						line.append(age)
					else:
						line.append(0)

					if gender_valid == True:
						if ambious_age_gender_flag:
							if idToGender[seqIdtmp] == 'Male':
								line.append(-1)
							elif idToGender[seqIdtmp] == 'Female':
								line.append(1)
							else:
								line.append(0)
						else:
							if idToGender[seqIdtmp] == 'Male':
								line.append(1)
							elif idToGender[seqIdtmp] == 'Female':
								line.append(0)
							else:
								print('*******************error**********************')
								exit()
					else:
						if idToGender[seqIdtmp] == 'Male':
							line.append(0)
						elif idToGender[seqIdtmp] == 'Female':
							line.append(0)
						else:
							print('*******************error**********************')
							exit()

					if continent_valid:
						if idToContinent[seqIdtmp] == 'Asia':
							line.append(1)
						elif idToContinent[seqIdtmp] == 'Europe':
							line.append(2)
						elif idToContinent[seqIdtmp] == 'NorthAmerica':
							line.append(3)
						elif idToContinent[seqIdtmp] == 'SouthAmerica':
							line.append(4)
						elif idToContinent[seqIdtmp] == 'Africa':
							line.append(5)
						elif idToContinent[seqIdtmp] == 'Oceania':
							line.append(6)
					else:
						line.append(0)

					newList.append(line)
					newSeqId.append(seqIdtmp)

	return newList


print("reading in data " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

referenceSeq = findReferenceSeq()

readInAndFormatData()


print("processing snps, formatting data " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

if LoadOfflineIndiciesToKeepFile == 0:  
	findColumnsWithoutSNPs()
else:
	with open(OffIndiciesFilePath, "r") as f:
		data1 = f.readlines()

	indiciesToKeep[0] = True  

	for i in data1:
		indiciesToKeep[int(i)] = True


feature_nn_size = len(indiciesToKeep)

print('*************************************read feature size: ', feature_nn_size)

dataList = removeOtherIndices(indiciesToKeep)

print("# sequences before blacklisting")
print(len(dataList))

dataList = removeAmbiguous()

print("# sequences after blacklisting")
print(len(dataList))

processDataLenth = len(dataList)

print("new seq Id lenth: ", len(newSeqId))

if dataModetest == 1:
	seq_num = len(dataList)
	seq_len = len(dataList[0]) - 1  # num of indices

	refseq_original = dataList[0][1:]

	refseq = finalLineRef[0:]

	if refseq != refseq_original:
		print(refseq)
		print(refseq_original)
		for i in range(0, len(refseq)):
			if refseq[i] != refseq_original[i]:
				print(i, refseq[i], '/', refseq_original[i])
	for index in range(1, seq_num):  # other than the reference
		queryline = dataList[index][1:]
		for queryIndex in range(
				len(queryline) - feature_plus):  
			query = queryline[queryIndex]
			ref = refseq[queryIndex]
			if inputModeCharacter: 
				if gene_valid:
					if query == ref or query == '-':
						dataList[index][queryIndex + 1] = ord(query)   # query
					else:
						dataList[index][queryIndex + 1] = ord(query)   # query
				else:
					if query == ref or query == '-':
						dataList[index][queryIndex + 1] = 0
					else:
						dataList[index][queryIndex + 1] = 0
			else:
				if query == ref or query == '-':
					dataList[index][queryIndex + 1] = 0
				else:
					dataList[index][queryIndex + 1] = 1
	# update the reference line
	for i in range(len(queryline)):
		dataList[0][i + 1] = 0

# headers are the original genome locations
headers = list(indiciesToKeep.keys())



print("setting up training " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

pima = pd.DataFrame(dataList, columns=headers)

if dataModetest == 0:
	# nucleotide symbols which can appear
	categories = ['A', 'C', 'G', 'T', '-']

	# one hot encoding of all headers other than the first which is the lineage
	dummyHeaders = headers[1:]

	# add extra rows to ensure all of the categories are represented, as otherwise
	# not enough columns will be created when we call get_dummies
	for i in categories:
		line = [i] * len(dataList[0])
		pima.loc[len(pima)] = line

	# get one-hot encoding
	pima = pd.get_dummies(pima, columns=dummyHeaders)

	# get rid of the fake data we just added
	pima.drop(pima.tail(len(categories)).index, inplace=True)

feature_cols = list(pima)
# print(feature_cols)

# remove the last column from the data frame. This is because we are trying to predict these values.
h = feature_cols.pop(0)
X = pima[feature_cols]
y = pima[h]

X = X.drop(X.index[0])
y = y.drop(y.index[0])

print('del dataList')
del dataList
gc.collect()

if loadTrainedModelFlag == 0:
	
# train and test data report
severe_train = 0
severe_test = 0
mild_train = 0
mild_test = 0
for y_label in y_train:
	if y_label == 'SevereLabel' or y_label == 1:
		severe_train += 1
	elif y_label == 'MildLabel' or y_label == 0:
		mild_train += 1
	else:
		print('Error, invalid label')
for y_label in y_test:
	if y_label == 'SevereLabel' or y_label == 1:
		severe_test += 1
	elif y_label == 'MildLabel' or y_label == 0:
		mild_test += 1
	else:
		print('Error, invalid label')

print('Total train and test:', severe_train + mild_train, severe_test + mild_test)
print('Total severe and mild:', severe_train + severe_test, mild_train + mild_test)
print('Train severe and mild:', severe_train, mild_train)
print('Test severe and mild:', severe_test, mild_test)
total_samples = severe_train + mild_train + severe_test + mild_test



if loadTrainedModelFlag == 0:

	if testRF == 1:
		RF = RandomForestClassifier(n_estimators=500, n_jobs=14, criterion='gini', \
										min_samples_split=4, min_samples_leaf=1,\
										max_depth=None,max_features='auto', random_state = 20)  
		print('RandomForestClassifier model is training .....')
		print("training start " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)
		if k_fold_train == 0:
			start_time = time.time()
			RF.fit(X_train, y_train)
			end_time = time.time()
			print("RandomForestClassifier: {:.3f}".format(end_time - start_time))
			print("training end " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)
			print('RandomForestClassifier model has been trained .....')

			y_predRF = RF.predict(X_test)
			p = precision_score(y_test, y_predRF, average='weighted')  #
			r = recall_score(y_test, y_predRF, average='weighted')
			f1 = f1_score(y_test, y_predRF, average='weighted')
			print('precision: ', p)
			print('recall: ', r)
			print('f1_score: ', f1)
			print("--------------------------------------------")
			print("Confusion Matrix of RandomForest")
			cnf_matrix_rf = metrics.confusion_matrix(y_test, y_predRF)
			print("--------------------------------------------")
			print("Classification report")
			print(metrics.classification_report(y_test, y_predRF, digits=5))
			print(cnf_matrix_rf)
			# for AUC
			y_Pred_n_2 = RF.predict_proba(X_test)  # test
			y_Pred = []
			for i in range(len(y_test)):
				y_Pred.append(y_Pred_n_2[i][0])
			y_True = []
			for y in y_test:
				if y == RF.classes_[0]:
					y_True.append(1)
				else:
					y_True.append(0)
			print('AUC: ', roc_auc_score(y_True, y_Pred))

			print('\n **************Step2 Grid Search CV******************\n')

			criterion = ['gini']  # 'entropy'
			n_estimators = [100,200,300,400,500,1000,1500,2000,2500,3000]
			max_depth = [None]# ,1,2,4,6]
			min_samples_split = [2,4,6] #
			min_samples_leaf = [1,2,4] #
			max_features = ['auto'] #,'sqrt', 'log2'] #
			# n_estimators = [100,1000,1500,2000]
			# min_samples_split = [2,4,6,8]
			# min_samples_leaf = [1,2,4,8]
			#criterion_best:  gini    n_estimators_best: 1000
			#min_samples_split_best:  4    min_samples_leaf_best: 1
			#max_depth_best: None    max_features auto


			parameters = {'n_estimators': n_estimators, 'criterion': criterion, 'min_samples_split': min_samples_split, \
						  'min_samples_leaf': min_samples_leaf,'max_depth':max_depth, 'max_features':max_features}


			grid = GridSearchCV(estimator=RandomForestClassifier(), param_grid=parameters, cv=10, scoring='f1_weighted',n_jobs=-1)
			print('grid fit start...')
			start_time = time.time()
			grid.fit(X_train, y_train)
			end_time = time.time()
			print("GridSearchCV: {:.3f}".format(end_time - start_time))
			print('cv_results_', grid.cv_results_)  
			print('best_score_', grid.best_score_) 
			print('best_params_', grid.best_params_)  
			print('best_estimator_', grid.best_estimator_) 
			criterion_best = grid.best_params_['criterion']
			n_estimators_best = grid.best_params_['n_estimators']
			max_depth_best = grid.best_params_['max_depth']
			min_samples_split_best = grid.best_params_['min_samples_split']
			min_samples_leaf_best = grid.best_params_['min_samples_leaf']
			max_features_best = grid.best_params_['max_features']
			print('criterion_best: ', criterion_best,'   n_estimators_best:',n_estimators_best)
			print('min_samples_split_best: ', min_samples_split_best, '   min_samples_leaf_best:', min_samples_leaf_best)
			print('max_depth_best:', max_depth_best, '   max_features', max_features_best)

			print('\n **************Step3 After Grid Search CV******************\n')
			RF2 = RandomForestClassifier(n_estimators=n_estimators_best, n_jobs=14, criterion=criterion_best, \
										min_samples_split=min_samples_split_best, min_samples_leaf=min_samples_leaf_best,\
										max_depth=max_depth_best,max_features=max_features_best)  # n_estimators default: 100
			start_time = time.time()
			RF2.fit(X_train, y_train)
			end_time = time.time()
			print("RandomForestClassifier2: {:.3f}".format(end_time - start_time))
			print("training end " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)
			print('RandomForestClassifier2 model has been trained .....')

			modelsaved = os.path.join('./RandomForestBestModel.joblib')  
			joblib.dump(RF2, modelsaved, compress=9)
			print('saving model to',modelsaved)

			print('________________out testing_______After grid search cv_____________size', len(X_test))
			y_predRF = RF2.predict(X_test)
			p = precision_score(y_test, y_predRF, average='weighted')  #
			r = recall_score(y_test, y_predRF, average='weighted')
			f1 = f1_score(y_test, y_predRF, average='weighted')
			print('precision2: ', p)
			print('recall2: ', r)
			print('f1_score2: ', f1)
			print("--------------------------------------------")
			print("Confusion Matrix of RandomForest2")
			cnf_matrix_rf = metrics.confusion_matrix(y_test, y_predRF)
			print("--------------------------------------------")
			print("Classification report2")
			print(metrics.classification_report(y_test, y_predRF, digits=5))
			print(cnf_matrix_rf)
			# for AUC
			y_Pred_n_2 = RF2.predict_proba(X_test)  # test
			y_Pred = []
			for i in range(len(y_test)):
				y_Pred.append(y_Pred_n_2[i][0])
			y_True = []
			for y in y_test:
				if y == RF.classes_[0]:
					y_True.append(1)
				else:
					y_True.append(0)
			print('AUC: ', roc_auc_score(y_True, y_Pred))

			if shap_flag:
				explainer = shap.TreeExplainer(RF)
				shap.initjs()  
				shap_values = explainer.shap_values(X_train)
				# shap.force_plot(explainer.expected_value[0], shap_values[0][100,:], X_train.iloc[100,:],matplotlib=True,show=False)
				# shap.force_plot(explainer.expected_value[0], shap_values[0],matplotlib=True)
				shap.summary_plot(shap_values, X_train, plot_type="bar", max_display=75)
				shap.summary_plot(shap_values[1], X_train, plot_type="violin", max_display=75,
								  sort=True)  # plot_type: violin, layered_violin
		
		else:
			print('_____________________train with k fold inside validation____________________')
			# X = X_train
			# y = y_train
			print('para set:', feature_nn_size - 1, max_ter)
			kf = KFold(n_splits=k_set, random_state = 30)
			kfold_id = 0
			acc_ave = 0
			ad_ave = 0
			f1score_ave = 0
			p_ave = 0
			start_time = time.time()
			acc_all = []
			for train_index, test_index in kf.split(X_train):
				kfold_id += 1
				print('kfold_id:', kfold_id)
				# print('train_index', train_index, 'test_index', test_index)
				# train_X, train_y = X[train_index], y[train_index]
				# test_X, test_y = X[test_index], y[test_index]
				train_X = X_train.iloc[train_index, :]
				train_y = y_train.iloc[train_index]
				test_X = X_train.iloc[test_index, :]
				test_y = y_train.iloc[test_index]
				RF.fit(train_X, train_y)
				# acc = RF.score(test_X, test_y)
				# print("RandomForestClassifier k fold validation",kfold_id, acc)
				# acc_ave += acc
				# acc_all.append(acc)
				#
				predict_y = RF.predict(test_X)
				p = precision_score(test_y, predict_y, average='weighted')  #
				r = recall_score(test_y, predict_y, average='weighted')
				f1 = f1_score(test_y, predict_y, average='weighted')
				acc = accuracy_score(test_y, predict_y)
				print('acc: ', acc)
				print('precision: ', p)
				print('recall: ', r)
				print('f1_score: ', f1)
				print("RF k fold validation", kfold_id, acc)
				print(metrics.classification_report(test_y, predict_y, digits=5))
			end_time_i = time.time()
			print("i  {:.3f}".format(end_time_i - start_time_i))


		if test_flag == 1:
			# print("RandomForestClassifier score", RF.score(X_test, y_test))
			print('________________out testing____________________size', len(X_test))
			y_predRF = RF.predict(X_test)
			p = precision_score(y_test, y_predRF, average='weighted')  #
			r = recall_score(y_test, y_predRF, average='weighted')
			f1 = f1_score(y_test, y_predRF, average='weighted')
			print('precision: ', p)
			print('recall: ', r)
			print('f1_score: ', f1)
			print("--------------------------------------------")
			print("Confusion Matrix of RandomForest")
			cnf_matrix_rf = metrics.confusion_matrix(y_test, y_predRF)
			print("--------------------------------------------")
			print("Classification report")
			print(metrics.classification_report(y_test, y_predRF, digits=5))
			print(cnf_matrix_rf)
			# for AUC
			y_Pred_n_2 = RF.predict_proba(X_test)  # test
			y_Pred = []
			for i in range(len(y_test)):
				y_Pred.append(y_Pred_n_2[i][0])
			y_True = []
			for y in y_test:
				if y == RF.classes_[0]:
					y_True.append(1)
				else:
					y_True.append(0)
			print('AUC: ', roc_auc_score(y_True, y_Pred))

	
	if testXGBoost == 1:
		print('XGBoost model is training .....')
		model = XGBClassifier(learning_rate=0.1, n_jobs=14, gamma=0.1, \
										 max_depth=6,
										 min_child_weight=1, \
										 n_estimators=500,
										 reg_alpha=0)
		# (max_depth=depth_max, learning_rate=0.1, n_estimators=iter_max, objective='binary:logistic')
		print('XGBoost model is training step2.....')
		print("training start " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)
		if k_fold_train == 0:

			print('\n **************Step1 Before Grid Search CV******************\n')
			print('train data percent:', 1 - testing_percentage)
			start_time = time.time()
			model.fit(X_train, y_train)
			end_time = time.time()
			print("XGBoost : {:.3f}".format(end_time - start_time))
			print("training end " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)
			print('XGBoost model has been trained .....')

			print('________________out testing_______Default_____________size', len(X_test))
			y_predRF = model.predict(X_test)
			p = precision_score(y_test, y_predRF, average='weighted')  #
			r = recall_score(y_test, y_predRF, average='weighted')
			f1 = f1_score(y_test, y_predRF, average='weighted')
			print('precision: ', p)
			print('recall: ', r)
			print('f1_score: ', f1)
			print("--------------------------------------------")
			print("Confusion Matrix of XGBoost")
			cnf_matrix_rf = metrics.confusion_matrix(y_test, y_predRF)
			print("--------------------------------------------")
			print("Classification report")
			print(metrics.classification_report(y_test, y_predRF, digits=5))
			print(cnf_matrix_rf)
			# for AUC
			y_Pred_n_2 = model.predict_proba(X_test)  # test
			y_Pred = []
			for i in range(len(y_test)):
				y_Pred.append(y_Pred_n_2[i][0])
			y_True = []
			for y in y_test:
				if y == model.classes_[0]:
					y_True.append(1)
				else:
					y_True.append(0)
			print('AUC: ', roc_auc_score(y_True, y_Pred))

			print('\n **************Step2 Grid Search CV******************\n')


			n_estimators = [500]  # [100, 200, 500,1000,1500,2000,2500,3000,5000] best 500
			learning_rate = [0.1] # [0.05, 0.1, 0.2, 0.5, 0.75, 1.0]  # best 0.1
			max_depth = [6] # [0, 1, 2, 4, 6, 8]  # ,1,2,4,6] # best 6
			gamma = [0.1] #[0, 0.1, 0.2, 0.5, 1.0] # 0.1 or None best
			min_child_weight = [1] # [0, 1, 2, 4, 6, 10]  # 1 best
			reg_alpha = 0 # [0, 0.1,0.25, 0.5,0.75,1.0]


			parameters = {'learning_rate': learning_rate, 'gamma': gamma, 'max_depth': max_depth, \
						  'min_child_weight': min_child_weight, 'n_estimators': n_estimators, 'reg_alpha': reg_alpha}


			grid = GridSearchCV(estimator=XGBClassifier(), param_grid=parameters, cv=10, scoring='f1_weighted',
								n_jobs=15)
			print('XGBoost grid fit start...')
			start_time = time.time()
			grid.fit(X_train, y_train)
			end_time = time.time()
			print("GridSearchCV: {:.3f}".format(end_time - start_time))
			print('cv_results_', grid.cv_results_)  
			print('best_score_', grid.best_score_) 
			print('best_params_', grid.best_params_) 
			print('best_estimator_', grid.best_estimator_)  
			learning_rate_best = grid.best_params_['learning_rate']
			gamma_best = grid.best_params_['gamma']
			max_depth_best = grid.best_params_['max_depth']
			min_child_weight_best = grid.best_params_['min_child_weight']
			n_estimators_best = grid.best_params_['n_estimators']
			reg_alpha_best = grid.best_params_['reg_alpha']
			print('learning_rate_best: ', learning_rate_best, '   gamma_best:', gamma_best)
			print('max_depth_best: ', max_depth_best, '   min_child_weight:',min_child_weight_best)
			print('n_estimators_best:', n_estimators_best, '   reg_alpha_best', reg_alpha_best)

			print('\n **************Step3 After Grid Search CV******************\n')
			model2 = XGBClassifier(learning_rate=learning_rate_best, n_jobs=14, gamma=gamma_best, \
										 max_depth=max_depth_best,
										 min_child_weight=min_child_weight_best, \
										 n_estimators=n_estimators_best,
										 reg_alpha=reg_alpha_best)  
			start_time = time.time()
			model2.fit(X_train, y_train)
			end_time = time.time()
			print("XGBoost2 : {:.3f}".format(end_time - start_time))
			print("training end " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)
			print('XGBoost model has been trained .....')

			modelsaved = os.path.join('./XGBoostModel.joblib')  
			joblib.dump(model2, modelsaved, compress=9)
			print('saving model to', modelsaved)

			print('________________out testing_______After grid search cv_____________size', len(X_test))
			y_predRF = model2.predict(X_test)
			p = precision_score(y_test, y_predRF, average='weighted')  #
			r = recall_score(y_test, y_predRF, average='weighted')
			f1 = f1_score(y_test, y_predRF, average='weighted')
			print('precision2: ', p)
			print('recall2: ', r)
			print('f1_score2: ', f1)
			print("--------------------------------------------")
			print("Confusion Matrix of XGBoost2")
			cnf_matrix_rf = metrics.confusion_matrix(y_test, y_predRF)
			print("--------------------------------------------")
			print("Classification report2")
			print(metrics.classification_report(y_test, y_predRF, digits=5))
			print(cnf_matrix_rf)
			# for AUC
			y_Pred_n_2 = model2.predict_proba(X_test)  # test
			y_Pred = []
			for i in range(len(y_test)):
				y_Pred.append(y_Pred_n_2[i][0])
			y_True = []
			for y in y_test:
				if y == model2.classes_[0]:
					y_True.append(1)
				else:
					y_True.append(0)
			print('AUC: ', roc_auc_score(y_True, y_Pred))

		else:
			print('_____________________XGBoost train with k fold inside validation____________________')
			# X = X_train
			# y = y_train
			print('para set:', feature_nn_size - 1)
			kf = KFold(n_splits=k_set, random_state = 30)
			kfold_id = 0
			acc_ave = 0
			ad_ave = 0
			f1score_ave = 0
			p_ave = 0
			start_time = time.time()
			acc_all = []
			for train_index, test_index in kf.split(X_train):
				start_time_i = time.time()
				kfold_id += 1
				print('kfold_id:', kfold_id)
				# print('train_index', train_index, 'test_index', test_index)
				# train_X, train_y = X[train_index], y[train_index]
				# test_X, test_y = X[test_index], y[test_index]
				train_X = X_train.iloc[train_index, :]
				train_y = y_train.iloc[train_index]
				test_X = X_train.iloc[test_index, :]
				test_y = y_train.iloc[test_index]
				model.fit(train_X, train_y)

				predict_y = model.predict(test_X)
				p = precision_score(test_y, predict_y, average='weighted')  #
				r = recall_score(test_y, predict_y, average='weighted')
				f1 = f1_score(test_y, predict_y, average='weighted')
				acc = accuracy_score(test_y, predict_y)
				print('acc: ', acc)
				print('precision: ', p)
				print('recall: ', r)
				print('f1_score: ', f1)
				print("XGBoost k fold validation", kfold_id, acc)
				print(metrics.classification_report(test_y, predict_y, digits=5))

			end_time_i = time.time()
			print("i  {:.3f}".format(end_time_i - start_time_i))

sys.exit()
