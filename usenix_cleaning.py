import numpy as np 
import os
import GEOparse
import pandas as pd
import networkx as nx
from sklearn.decomposition import PCA
from scipy import optimize
import csv
#import matplotlib.pyplot as plt

gse = GEOparse.get_GEO(filepath='datasets/GSE68951_family.soft')

df = pd.read_csv('GSM_conversion.txt', header=None, delim_whitespace=True)
df.rename(columns={0: 'sample', 1: 'patient', 2: 'time'}, inplace=True)

def two_timepoints(t1, t2, gse, df):
	filtered = df[(df['time'] == t1) | (df['time'] == t2)]
	e = filtered['patient'].value_counts()
	filtered = filtered[filtered['patient'].isin(e[e>1].index)]
	t1_df = filtered[filtered['time'] == t1]
	t2_df = filtered[filtered['time'] == t2]
	return t1_df, t2_df


def maximum_matching(t1s, t2s, pca):
	matching_matrix = np.zeros((len(t1s), len(t2s)))
	t1_pcas = []
	t2_pcas = []

	for sample1 in list(t1s['sample']):
		t1 = gse.gsms[sample1].table['VALUE'].to_numpy().reshape(1,-1)
		t1_pca = pca.transform(t1)
		t1_pcas.append(t1_pca)
	for sample2 in list(t2s['sample']):
		t2 = gse.gsms[sample2].table['VALUE'].to_numpy().reshape(1,-1)
		t2_pca = pca.transform(t2)
		t2_pcas.append(t2_pca)

	for i in range(matching_matrix.shape[0]):
		for j in range(matching_matrix.shape[1]):
			pca_dist = np.linalg.norm(t1_pcas[i]-t2_pcas[j])
			matching_matrix[i][j] = pca_dist
	rows, cols = optimize.linear_sum_assignment(matching_matrix)
	count = 0
	for r, c in zip(rows, cols):
		if r == c:
			count += 1
	return(count/len(rows))


def apply_pca(n_comp, df):
	for_pca = []

	for sample1 in list(df['sample']):
		t1 = gse.gsms[sample1].table['VALUE'].to_numpy()
		for_pca.append(t1)
	#for sample2 in list(t2s['sample']):
	#	t2 = gse.gsms[sample2].table['VALUE'].to_numpy()
	#	for_pca.append(t2)
	arr = np.stack(for_pca, axis=0)
	print(arr.shape)
	pca = PCA(n_components=n_comp, svd_solver='arpack')
	pca.fit(arr)

	#accuracy = maximum_matching(t1s, t2s, pca)
	#print(n_comp, accuracy)
	return(pca)


t1, t2 = two_timepoints(1,2,gse,df)
n_comps = np.linspace(4, 100, 25)

accs = []
for n_comp in n_comps:
	pca = apply_pca(int(n_comp), df)
	acc = maximum_matching(t1, t2, pca)
	print(n_comp, acc)
	with open("temp/two_timepoints.csv", 'a', newline='\n') as fd:
		fieldnames = ['n_comp', 'acc']
		writer = csv.DictWriter(fd, fieldnames=fieldnames)
		writer.writerow({"n_comp": n_comp, "acc": acc})
#plt.plot(n_comps, accs)







