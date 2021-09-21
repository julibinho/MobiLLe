"""
Evaluating the quality of any clustering tool even on sequences that we dont know the real cluster labeles.
"""
#Have to install optparse and tqdm

import sys
import math
import time
import resource
import Levenshtein
from Levenshtein import distance
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from collections import Counter
from optparse import OptionParser
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation
#####################################################################
def read_file (nomFi):
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	return lines
#####################################################################
def delete_duplicate(nomFi):
	lines = read_file (nomFi)
	uniq_seq_dico = {}
	filtered_lines =[]
	dup_corresp = {}
	for l in range(0,len(lines)):
		seq= lines[l].split("\t")
		dup_id = seq[2].rstrip()+"_"+seq[4].rstrip()+"_"+seq[5].rstrip()
		if dup_id in  dup_corresp.keys() :
			dup_corresp[dup_id].append(seq[0])
		else :
			dup_corresp[dup_id] = [seq[0]]
			filtered_lines.append(lines[l])
	for key in dup_corresp.keys():
		#print dup_corresp[key]
		uniq_seq_dico[dup_corresp[key][0]] = dup_corresp[key][1:]
	return uniq_seq_dico,filtered_lines
#####################################################################
def read_vjunction(lines):
	seq_V_CDR3_J = {}
	len_max_CDR3 = 0
	len_max_J = 0
	for l in lines:
		seq = l.split("\t")
		#seq_V_CDR3_J[seq id] = [Functionality,V_id,J_id,CDR3_seq_aa,J_seq_nt]
		seq_V_CDR3_J[seq[0]] = [seq[1].rstrip(),seq[2].rstrip(),seq[3].rstrip(),seq[4].rstrip(),seq[5].rstrip()]
		if len(seq[3].rstrip()) > len_max_CDR3 :
			len_max_CDR3 = len(seq[3].rstrip())
		if len(seq[4].rstrip()) > len_max_J :
			len_max_J = len(seq[4].rstrip())
	return seq_V_CDR3_J, len_max_CDR3,len_max_J
#####################################################################
# read the clustering result
def readClusteringResults(nomFi):
	lines = read_file (nomFi)
	Clustering_lables = {}
	for l in range(len(lines)):
		Seq_nom = lines[l].split("\t")[1].rstrip().split(" ")
		cluster = lines[l].split("\t")[0].rstrip()
		if cluster in Clustering_lables.keys():
			Clustering_lables[cluster].append(Seq_nom)
		else:
			Clustering_lables[cluster] = Seq_nom
	#print(Clustering_lables)
	return Clustering_lables
#####################################################################
def filter_clustering_output (Clustering_lables,uniq_seq_dico):
	filtered_clustering_label = {}
	for cluster in Clustering_lables.keys():
		filtered_clustering_label[cluster] = []
		for seq in Clustering_lables[cluster]:
			if seq in uniq_seq_dico.keys():
				filtered_clustering_label[cluster].append(seq)
	return filtered_clustering_label
#####################################################################
def hamming_similarity(s1,s2,l):
    # computes the hamming similarity between two sequences
    # s1 and s2 are sequencs and l is the length of the sequences
    delta = 0.0
    for a1, a2 in zip(s1, s2):
        if a1 == a2:
            delta = delta + 1
    return (delta / l) * 100
#####################################################################
#computes the sequence similairy between s1 with length l1 and s2
# with length l2. Hamiing similarity if l1==l2 and edit distance/levenstein 
# distance if l1!=l2
def get_similarity_score(s1,s2, l1, l2):
	l_avg=(l1+l2)/2.0
	if l1==l2:
	   sim=hamming_similarity(s1,s2,l1)
	else:
	    if l1>l2 :
	        long_seq = s1
	        short_seq = s2
	    else : 
	        long_seq = s2
	        short_seq = s1
	    five_prime = long_seq[ : len(short_seq)]
	    three_prime = long_seq[-(len(short_seq)) :]
	    sim_five=hamming_similarity(five_prime,short_seq,len(short_seq))
	    sim_three=hamming_similarity(three_prime,short_seq,len(short_seq))
	    sim_lev=100.0*(1.0-distance(s1,s2)/l_avg)
	    sim = max(round(sim_five),round(sim_three),round(sim_lev))
	return sim
#####################################################################
def CalculateMedoid(dico_vjunc,Dicoresult):
	centroid={}
	#print (Dicoresult)
	for key in Dicoresult.keys() :
		listloc=[]
		for seq in Dicoresult[key]:
			if seq.rstrip() in dico_vjunc.keys():
				listloc.append(dico_vjunc[seq.rstrip()][3])
		if len(listloc) != 0 :
			centroid[key]=Levenshtein.median(listloc)
	#print ("centroid",centroid)
	return centroid
#####################################################################
def Creat_dico_neighbour(Dicocentroid):
	Cluster_list=[]
	Centroid_list=[]
	dicoNeighbour={}
	#print (Dicocentroid)	
	for cluster in Dicocentroid.keys():
		Cluster_list.append(cluster)
		Centroid_list.append(Dicocentroid[cluster])
	list__for_final_elem=[]
	for i in range(len(Centroid_list)-1) :	
		listloc=[]
		for j in range(i+1,len(Centroid_list)):
			listloc.append(Levenshtein.distance(Centroid_list[i],Centroid_list[j]))
			if j == len(Centroid_list)-1:
				list__for_final_elem.append(Levenshtein.distance(Centroid_list[i],Centroid_list[j]))
		argmin=listloc.index(min(listloc))
		if Cluster_list[i] not in dicoNeighbour.keys():
			dicoNeighbour[Cluster_list[i]]=Cluster_list[argmin+i+1]
		if Cluster_list[argmin+i] not in dicoNeighbour.keys():
			dicoNeighbour[Cluster_list[argmin+i+1]]=Cluster_list[i]

	if dicoNeighbour != {} :
		dicoNeighbour[Cluster_list[-1]]=Cluster_list[list__for_final_elem.index(min(list__for_final_elem))]
	else:
		return dicoNeighbour
	return dicoNeighbour
#####################################################################
def run_silhouette(Dicofasta,Dicocentroid,Dicoresult,DicoNeighbour, len_max_CDR3,len_max_J,dico_vjunc):
	i=0
	j=0

	new_dico_res = Dicoresult
	new_dico_neig = DicoNeighbour
	new_dico_cent = Dicocentroid
	while (i <10000 and j<5 and new_dico_neig != {}) :
		before = new_dico_res
		a,b,c = silhouette(Dicofasta,new_dico_cent,new_dico_res,new_dico_neig, len_max_CDR3,len_max_J,dico_vjunc)
		if a == before :
			j +=1
		else :
			j = 0
		new_dico_res = a
		new_dico_neig = b
		new_dico_cent = c
		i +=1

	return new_dico_res

#####################################################################

# SILHOUETTE

def silhouette(Dicofasta,Dicocentroid,Dicoresult,DicoNeighbour, len_max_CDR3,len_max_J,dico_vjunc):
	summe=0
	for cluster in Dicoresult.keys() :
		for seq in Dicoresult[cluster]:
			dist_intra = 0
			dist_neighb = {}
			for seq_same_clust in Dicoresult[cluster] :
				if seq_same_clust != seq :
					V_component = 0
					# same V, without considering allel
					if Dicofasta[seq][0].split("*")[0] != Dicofasta[seq_same_clust][0].split("*")[0] :
						V_component += 1
					# normalize the distance depending on the longest seq in the whole repertoire.
					CDR_component = float(Levenshtein.distance(Dicofasta[seq][2].split(" ")[0],Dicofasta[seq_same_clust][2].split(" ")[0])) / len_max_CDR3
					J_component = (100-get_similarity_score(Dicofasta[seq][3],Dicofasta[seq_same_clust][3], len(Dicofasta[seq][3]), len(Dicofasta[seq_same_clust][3])) )/float(len_max_J)
					#print ('J_component',J_component)
					dist_intra += CDR_component + V_component+J_component/3.0
			if len(Dicoresult[cluster]) != 1 :
				ai = float(dist_intra) /(len(Dicoresult[cluster]) - 1)
			else :
				ai = 0
			for seq_neighb in Dicoresult[DicoNeighbour[cluster]] :
				V_component_b = 0
				CDR_component_b = 0
				J_component_b = 0
				# same V, without considering allel
				if Dicofasta[seq][0].split("*")[0] != Dicofasta[seq_neighb][0].split("*")[0] :
					V_component_b += 1
				CDR_component_b = Levenshtein.distance(Dicofasta[seq][2].split(" ")[0],Dicofasta[seq_neighb][2].split(" ")[0]) / float(len_max_CDR3)
				J_component_b = (100 -get_similarity_score(Dicofasta[seq][3],Dicofasta[seq_neighb][3], len(Dicofasta[seq][3]), len(Dicofasta[seq_neighb][3])) )/float(len_max_J)
				#print('V_component_b',V_component_b,'CDR_component_b',CDR_component_b,'J_component_b',J_component_b)
				dist_neighb[seq_neighb] =CDR_component_b + V_component_b + J_component_b /3.0
			#print ('dist_neigh',dist_neighb,ai)
			bi = min(dist_neighb.values())
			#print ("================")
			#print (seq)
			#print ("ai = ",ai )
			if bi<ai :
				if len(Dicofasta[seq][2].split(" ")[0]) == len(Dicofasta[seq_neighb][2].split(" ")[0]):
					#print (len(Dicofasta[seq][2].split(" ")[0]),len(Dicofasta[seq_neighb][2].split(" ")[0]))
					 #print ("lalalalaala",ai,bi)
					#print("disonnnnn,Dicoresult",)
					
					to_move = (list(dist_neighb.keys())[list(dist_neighb.values()).index(bi)]) 
					#print("tomove",to_move)
					Dicoresult[DicoNeighbour[cluster]].remove(to_move)
					if len(Dicoresult[DicoNeighbour[cluster]]) == 0:
						del Dicoresult[DicoNeighbour[cluster]]
					Dicoresult[cluster].append(to_move)
					#print(Dicoresult)
					Dicocentroid = CalculateMedoid(dico_vjunc,Dicoresult)
					DicoNeighbour = Creat_dico_neighbour(Dicocentroid)
					#silhouette(Dicofasta,Dicocentroid,Dicoresult,DicoNeighbour, len_max_CDR3,len_max_J,dico_vjunc)


			return (Dicoresult,DicoNeighbour,Dicocentroid)
			#print ("sil = ", calculeSil(ai,bi))
			#summe+=calculeSil(ai,bi)
	#return summe/float(len(Dicofasta))

#####################################################################

def calculeSil(ai,bi):
    si=0
    if ai!=bi:
        if ai<bi:
            si=1-(ai/float(bi))
        else:
            si=(float(bi)/ai)-1
    return si

def write_clone_V_cdr3_(Dicoresult,dicoSeq,uniq_seq_dico,output_file):
	file_name = output_file.split(".")[0]+"_clone_V_CDR3_J.txt"
	filetowrite=open(file_name,"w")
	clone_number = 0
	for cluster in Dicoresult.keys():
		for seq in Dicoresult[cluster]:
				sequence = str(clone_number) + "\t" + str(seq) + "\t" +dicoSeq[seq][0] + "\t" + dicoSeq[seq][1] + "\t" +dicoSeq[seq][2].split(" ")[0] +"\t" +dicoSeq[seq][3]+"\t" +dicoSeq[seq][4]+ "\n"
				filetowrite.write(sequence)
				if len(uniq_seq_dico[seq]) != 0 :
					for dup in uniq_seq_dico[seq] :
						sequence = str(clone_number) + "\t" + str(dup) + "\t" +dicoSeq[seq][0] + "\t" + dicoSeq[seq][1] + "\t" +dicoSeq[seq][2].split(" ")[0] +"\t" +dicoSeq[seq][3]+"\t" +dicoSeq[seq][4]+ "\n"
						filetowrite.write(sequence)
		clone_number += 1

	filetowrite.close()
	return 0
####################################################################
def main():
	usage = "usage: refinement.py -f FastaFile -c ClusteringFile"
	parser = OptionParser(usage)
	parser.add_option("-f", "--FastaFile", dest="FastaFile",
	      help="read data from FILENAME")
	parser.add_option("-c", "--ClusteringFile",dest="ClusteringFile",
	      help="read data from ClusteringFile")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 5:
		parser.error("incorrect number of arguments")
	FastaFile = options.FastaFile
	ClusteringFile = options.ClusteringFile
	uniq_seq_dico,filtered_lines = delete_duplicate(FastaFile)

	dico_vjunc, len_max_CDR3,len_max_J =read_vjunction(filtered_lines)
	
	dico_result = readClusteringResults(ClusteringFile)
	filtered_clustering_label = filter_clustering_output (dico_result,uniq_seq_dico)
	#print("filtered_clustering_label",filtered_clustering_label)

	dico_centroid = CalculateMedoid(dico_vjunc,filtered_clustering_label)
	dico_neighbour = Creat_dico_neighbour(dico_centroid)
	#print ("dicoNeighbour",dico_neighbour)
	Dicoresult = run_silhouette(dico_vjunc,dico_centroid,filtered_clustering_label,dico_neighbour,len_max_CDR3,len_max_J,dico_vjunc)
	write_clone_V_cdr3_(Dicoresult,dico_vjunc,uniq_seq_dico,FastaFile)

#####################################################################
if __name__ == "__main__":
	main()
  


