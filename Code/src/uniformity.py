#Have to install optparse and tqdm

import sys
import time
from refinement import *

#####################################################################

def get_clusters_uniformity(Dicofasta, Dicoresult, max_lens, typeDistanceT) :
	"""
	Recover the uniformity for all the clusters
	input Dicofasta: 		Dict()		key=sequence ID, value= Tuple(status, IGHV, IGHJ, CDR3 nuc sequence, CDR3 aa sequence)
	input Dicoresult: 		Dict() 		key=cluster label, value= list of unique sequence IDs within the cluster
	input max_lens:			Tuple(int, int, int) containg max lengths of CDR3, J, and V sequences
	input typeDistance: 	Tuple(int, int, int, int) type distance of each componnent
	"""
	uniformity = {}
	clone_number = 0
	for cluster in Dicoresult.keys() :
		seqs = Dicoresult[cluster]
		uniformity[clone_number] = computeUniformity(Dicofasta, seqs, max_lens, typeDistanceT)
		clone_number += 1
	return uniformity
		
#####################################################################

def save_cluster_uniformity (uniformity, output_file) :
	file_name = output_file+"_clone_uniformity.txt"
	filetowrite=open(file_name,"w")
	for cluster in uniformity.keys():
		filetowrite.write(str(cluster)+"\t"+str(uniformity[cluster])+"\n")
	filetowrite.close()
	return 0

####################################################################

def main():
	start_time = time.time()
	usage = "usage: uniformity.py -f FastaFile -l ClusteringFile -v VDistance -j JDistance -c CDR3Distance -d CombineDistance -m MergeSingleton "
	parser = OptionParser(usage)
	parser.add_option("-f", "--FastaFile", dest="FastaFile", help="read data from FILENAME")
	parser.add_option("-l", "--ClusteringFile",dest="ClusteringFile", help="read data from ClusteringFile")
	parser.add_option("-v", "--VDistance",dest="VDistance", help="method for calculating the V gene distance, 1- binaire, 2-levenstein,3- GIANA, 4-K-mers")
	parser.add_option("-j", "--JDistance",dest="JDistance", help="method for calculating the J gene distance, 1- binaire, 2-levenstein,3- GIANA, 4-K-mers")
	parser.add_option("-c", "--CDR3Distance",dest="CDR3Distance", help="method for calculating the CDR3 region distance, 1- binaire, 2-levenstein,3- GIANA, 4-K-mers")
	parser.add_option("-d", "--CombineDistance",dest="CombineDistance", help="method for combining three distances (1-mean, 2-harmonic mean,  ponderate mean(requires a Tuple with three weights), )")
	parser.add_option("-m", "--MergeSingleton", dest="MergeSingleton", help="Merging the singletons")
	parser.add_option("-o", "--output", dest="outputname", help="output name")


	(options, args) = parser.parse_args()
	if len(sys.argv) != 17:
		parser.error("incorrect number of arguments")

	FastaFile = options.FastaFile
	outputname = options.outputname
	ClusteringFile = options.ClusteringFile
	VDistance = options.VDistance
	JDistance = options.JDistance
	CDR3Distance = options.CDR3Distance
	CombineDistance = options.CombineDistance
	MergeSingleton = options.MergeSingleton
	

	typeDistanceT = (VDistance, CDR3Distance, JDistance, CombineDistance)

	uniq_seq_dico, filtered_lines = delete_duplicate(FastaFile)
	dico_vjunc, max_lens = read_vjunction(filtered_lines)

	(len_max_CDR3, len_max_J, len_max_V)= max_lens
	dico_result = readClusteringResults(ClusteringFile)
	filtered_clustering_label = filter_clustering_output (dico_result, uniq_seq_dico)
	
	######
	#save uniformity for each cluster
	######
	clusters_uniformity = get_clusters_uniformity(dico_vjunc, filtered_clustering_label, max_lens, typeDistanceT)
	save_cluster_uniformity(clusters_uniformity, outputname)
	
	print("The uniformity compute time : %s seconds " % (time.time() - start_time))
	
#####################################################################
if __name__ == "__main__":
	main()
  


