"""
The program gather all sequences with the same V gene and allele, J gene and
Junction AA with an identity percentage value higher than i into one clone.
"""

import sys
from optparse import OptionParser
import operator
from Levenshtein import distance as levenshtein_distance
import skbio 
from skbio import TabularMSA
from ClusterCDR3Only import FaIR_CDR3Only
import time
#from skbio.alignment.parasail import global_pairwise_align_protein
#import parasail
# =============================================================================
#               	Read IMGT high/vquest output 
# =============================================================================	

def read_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines


def filter_seq(lines):
	filtered_lines = []
	seq_unannotated = []
	for l in range(0,len(lines)):
		seq= lines[l].split("\t")
		if seq[2].rstrip() == "_" or seq[3].rstrip() == "_" or seq[4].rstrip() == "_" :
			seq_unannotated.append(seq[0]+ "\t" + seq[1] + "\t" + seq[2] + "\t" + seq[3] + "\t" + seq[4] + "\n")
		else :
			filtered_lines.append(lines[l])
	print ("Number of discarded sequences because of incompelete gene annoation by IMGT: ", len(seq_unannotated))
	print ("Number of sequences to be clustered : ", len(filtered_lines))
	return filtered_lines



def delete_duplicate(lines):
	uniq_seq_dico = {}
	filtered_lines =[]
	dup_corresp = {}

	for l in range(0,len(lines)):

		seq = lines[l].split("\t")

		dup_id = seq[2].rstrip()+"_"+seq[3].rstrip()+"_"+seq[4].rstrip()

		if dup_id in  dup_corresp.keys() :

			dup_corresp[dup_id].append(seq[0])
		else :
			dup_corresp[dup_id] = [seq[0]]
			filtered_lines.append(lines[l])

	for key in dup_corresp.keys():
		uniq_seq_dico[dup_corresp[key][0]] = dup_corresp[key][1:]
	return uniq_seq_dico,filtered_lines


#=============================================================================#

def group_same_VJ(lines):

	dico_same_VJ = {} # same Vgene, same J gene
	dicoSeq = {}

	for l in range(0,len(lines)):
		NumClone= lines[l].split("\t")
		
		dicoSeq[NumClone[0]] = [NumClone[1].rstrip(),NumClone[2].rstrip(),NumClone[3].rstrip(),NumClone[4].rstrip(),NumClone[5].rstrip(),NumClone[6].rstrip(),NumClone[7].rstrip(),NumClone[8].rstrip()]
		Clone_identity = ""
		Clone_identity = str(NumClone[2].split("*")[0]+"_"+NumClone[3].split("*")[0]) # same V gene  + same J gene
		if Clone_identity in dico_same_VJ.keys():
			dico_same_VJ[Clone_identity].append(NumClone[0])
		else:
			dico_same_VJ[Clone_identity] = [NumClone[0]]

	return dico_same_VJ,dicoSeq

# =============================================================================

def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))
# =============================================================================	
def group_clone_VJ_cdr3(dico_same_VJ,directory,dicoSeq,Clone_threshold):
	VJ_ID_diff_CDR3 = {}
	dicoclone_vj_cdr3 = {}
	loclist =[]
	for VJ_id in dico_same_VJ.keys():
		if VJ_id.split("_")[1] == "IGHJ6" : 
			#for J6 we toleart 1 AA of length difference
			loclist =[]
			for i in range(len(dico_same_VJ[VJ_id])):
				VJ_ID_diff_CDR3[VJ_id]={}
				loclist.append({'Id':dico_same_VJ[VJ_id][i] ,'CDR3' : dicoSeq[dico_same_VJ[VJ_id][i]][3] , 'Length' : len(dicoSeq[dico_same_VJ[VJ_id][i]][3]) })
			clusters=FaIR_CDR3Only(loclist,directory, th=Clone_threshold*100,tolerance=1, mth=Clone_threshold*100,split_size=250)
			for c in clusters :
				CDR3_seq=c[0]["CDR3"]
				VJ_ID_diff_CDR3[VJ_id][CDR3_seq] =[]
				for seq in c:
					VJ_ID_diff_CDR3[VJ_id][CDR3_seq].append(seq["Id"])

		else : 
			loclist =[]
			for i in range(len(dico_same_VJ[VJ_id])):
				VJ_ID_diff_CDR3[VJ_id]={}
				loclist.append({'Id':dico_same_VJ[VJ_id][i] ,'CDR3' : dicoSeq[dico_same_VJ[VJ_id][i]][3] , 'Length' : len(dicoSeq[dico_same_VJ[VJ_id][i]][3]) })
			clusters=FaIR_CDR3Only(loclist,directory,th=Clone_threshold*100,tolerance=0, mth=Clone_threshold*100,split_size=250)
			for c in clusters :
				CDR3_seq=c[0]["CDR3"]
				VJ_ID_diff_CDR3[VJ_id][CDR3_seq] =[]
				for seq in c:
					VJ_ID_diff_CDR3[VJ_id][CDR3_seq].append(seq["Id"])
	return VJ_ID_diff_CDR3



def write_clone_VJ_cdr3_all(VJ_ID_diff_CDR3,dicoSeq,output_file,Clone_threshold):
	file_name = output_file+"_sameVJ_noallele_CDR3_"+ Clone_threshold +".txt"
	filetowrite=open(file_name,"w")
	clone_number = 0
	for VJ in VJ_ID_diff_CDR3.keys():
		for cdr3 in VJ_ID_diff_CDR3[VJ]:
			for seq in VJ_ID_diff_CDR3[VJ][cdr3] :
				sequence = str(clone_number) + "\t" + str(seq) + "\t" + dicoSeq[seq.rstrip()][0] + "\t" +dicoSeq[seq.rstrip()][1] +"\t" +dicoSeq[seq.rstrip()][2]+ '\t' + dicoSeq[seq.rstrip()][3]+ '\t' + dicoSeq[seq.rstrip()][4]+ '\t' + dicoSeq[seq.rstrip()][5]+ '\t' + dicoSeq[seq.rstrip()][6]+ '\t' + dicoSeq[seq.rstrip()][7]+"\n"
				filetowrite.write(sequence)
			clone_number += 1
	filetowrite.close()
	return 0



def find_max_list(list):
    list_len = [len(i) for i in list]
    return(max(list_len))
#=============================================================================#

def main():
    start_time = time.time()
    usage = "python  initial_clustering.py  -i <formated IMGT highvquest statistics output> -o <output file name> -s <Clone identity between zero and one>  -a adress for temporary files \n "
    parser = OptionParser(usage)
    parser.add_option("-i", "--hv_stat_output", dest="hv_stat_output",
          help="formated IMGT highvquest statistics output")
    parser.add_option("-o", "--output_file_name", dest="output_file_name",
          help="the name for the file to write")
    parser.add_option("-s", "--Clone_threshold", dest="Clone_threshold",
          help="Clone identity percentage between 0 and 1")
    parser.add_option("-a", "--adress_tempo_file", dest="adress_tempo_file",
          help="adress for temporary files")
    (options, args) = parser.parse_args()
    if len(sys.argv) != 9:
        parser.error("incorrect number of arguments")

    IMGT_file = options.hv_stat_output
    output_file_name = options.output_file_name
    Clone_threshold = options.Clone_threshold
    adress_tempo_file = options.adress_tempo_file


    CloneID = read_file(IMGT_file)
    annotated_CloneID = filter_seq(CloneID)
    uniq_seq_dico,filtered_lines = delete_duplicate(annotated_CloneID)


    dico_same_VJ,dicoSeq = group_same_VJ(annotated_CloneID)

    
    VJ_ID_diff_CDR3 = group_clone_VJ_cdr3(dico_same_VJ,adress_tempo_file,dicoSeq,float(Clone_threshold),)

    write_clone_VJ_cdr3_all(VJ_ID_diff_CDR3,dicoSeq,output_file_name ,Clone_threshold)
    print("The initial clustering step execution time : %s seconds " % (time.time() - start_time))


#=============================================================================#

if __name__ == "__main__":
    main()
