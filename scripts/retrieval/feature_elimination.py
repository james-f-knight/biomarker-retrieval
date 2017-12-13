#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 10:52:48 2017

@author: james
"""

import csv
import sys
from collections import defaultdict
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from itertools import groupby
from sklearn.model_selection import cross_val_predict

#Read Normalised gene expression matrix and removes uniformative cols
def read_tidy(input_matrix):

    matrix = defaultdict(list)
    with open(input_matrix, 'r') as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            for (k, v) in list(row.items()):
                matrix[k].append(v)
    #Remove name column
    del matrix["Name"]
    
    return(matrix)


#Select data corresponding to the genes whose locus is contained in
#identifiers list
def select_row(identifiers, dictionary):

    ind = []
    dictionary_trim = {}
    for a in identifiers:
        ind.append(dictionary['Locus'].index(a))
    for a in dictionary:
        dictionary_trim[a] = []
        for z in ind:
            dictionary_trim[a].append(dictionary[a][z])
        dictionary_trim[a].append(dictionary[a][len(dictionary[a])-1])
    return(dictionary_trim)


data = read_tidy("Results/RNAseq_and_tiling_preprocessed")

biomarker_ids = []
with open(str("Results/polices_result/Experimentunion_model.txt"), 'r') as file:
    for line in file:
        line = line.strip()
        biomarker_ids.append(line.strip('"'))

biomarker_data = select_row(biomarker_ids, data)

attributes = {'CONTROL': 0, 'STRESS': 1}
Y = []
for a in biomarker_data:
    if a != 'Locus':
       Y.append(attributes[biomarker_data[a][-1]])

X = []

for a in biomarker_data:
    if a != 'Locus':
        Xi = []
        for z in range(len(biomarker_data[a])-1):
            Xi.append(float(biomarker_data[a][z]))
        X.append(Xi)

pairs = [[a,b] for a in range(len(biomarker_ids)) for b in range(len(biomarker_ids)) if b > a]

classifier = RandomForestClassifier(n_estimators = 3000, class_weight='balanced')

for n in range(len(pairs)):
    reduced_data = [[Xi[pairs[n][0]], Xi[pairs[n][1]]] for Xi in X]
    test = cross_val_predict(classifier, reduced_data, Y, cv=min(10, Y.count(1)))
    tp, tn, fp, fn = 0, 0, 0, 0
    for m in range(len(test)):
        if(test[m] == 1 & Y[m] == 1):
            tp = tp + 1
        elif(Y[m] == 1):
            fn = fn + 1
        elif(test[m] == 1):
            fp = fp + 1
        else:
            tn = tn + 1
    gmean = np.sqrt((tp/(tp + fp)*(tn/(tn + fn))))
    pairs[n].append(gmean)
    print(".")
print("\n")    
    

print("LocusTags            CV gmean")
print("---------------------------------------------------------------------")
for n in range(len(pairs)):
    print("{} {}    {}".format(biomarker_ids[pairs[n][0]], biomarker_ids[pairs[n][1]], pairs[n][2]))

