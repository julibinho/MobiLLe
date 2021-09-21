"""
The program gather all sequences with the same V gene and allele, J gene and
Junction AA with an identity percentage value higher than Z into one clone.
"""

import sys
from optparse import OptionParser
import operator
from Levenshtein import distance as levenshtein_distance
import skbio 
# =============================================================================
#               	Read IMGT high/vquest output 
# =============================================================================	

def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines

def filter_seq(lines,output_file):
	file_name = output_file + "_unannotated_seq.txt"
	filetowrite=open(file_name,"w")
	filtered_lines = []
	for l in range(0,len(lines)):
		seq= lines[l].split("\t")
		#print (seq)
		if seq[2].rstrip() == "_" or seq[3].rstrip() == "_" or seq[4].rstrip() == "_" :
			seq_unannotated =  seq[0]+ "\t" + seq[1] + "\t" + seq[2] + "\t" + seq[3] + "\t" + seq[4] + "\n"
			#print seq
			filetowrite.write(seq_unannotated)
		else :
			filtered_lines.append(lines[l])
	#print (filtered_lines)		
	return filtered_lines

def delete_duplicate(lines):
	uniq_seq_dico = {}
	filtered_lines =[]
	dup_corresp = {}
	#print "lines before filter : ",len(lines)
	for l in range(0,len(lines)):
		#print lines[l]
		seq= lines[l].split("\t")
		#print seq
		dup_id = seq[2].rstrip()+"_"+seq[3].rstrip()+"_"+seq[4].rstrip()
		#print dup_id
		if dup_id in  dup_corresp.keys() :
			#print seq[0]
			dup_corresp[dup_id].append(seq[0])
		else :
			dup_corresp[dup_id] = [seq[0]]
			filtered_lines.append(lines[l])
	#print dup_corresp,"dup_corresp"
	#print "lines after filter : ",len(filtered_lines)
	for key in dup_corresp.keys():
		#print dup_corresp[key]
		uniq_seq_dico[dup_corresp[key][0]] = dup_corresp[key][1:]
	#print (uniq_seq_dico)
	return uniq_seq_dico,filtered_lines


#=============================================================================#

def group_same_VJ(lines,formatted_file):
	#print lines
	dico_same_VJ = {} # same Vgene, same J gene
	dicoSeq = {}
	#print lines
	for l in range(0,len(lines)):
		NumClone= lines[l].split("\t")
		dicoSeq[NumClone[0]] = [NumClone[1].rstrip(),NumClone[2].rstrip(),NumClone[3].rstrip(),NumClone[4].rstrip()]
		Clone_identity = ""
		#Clone_identity = str(NumClone[1]+"_"+NumClone[2]) # same V gene and allele + same J gene and allele
		#Clone_identity = str(NumClone[1]+"_"+NumClone[2].split("*")[0]) # same V gene and allele + same J gene
		Clone_identity = str(NumClone[2].split("*")[0]+"_"+NumClone[3].split("*")[0]) # same V gene  + same J gene
		if Clone_identity in dico_same_VJ.keys():
			dico_same_VJ[Clone_identity].append(NumClone[0])
		else:
			dico_same_VJ[Clone_identity] = [NumClone[0]]
	#print "dico_same_VJ",len(dico_same_VJ['IGHV3-74_IGHJ4'])
	#print "dico seq", dicoSeq
	return dico_same_VJ,dicoSeq

# =============================================================================

def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))
# =============================================================================	
def group_clone_VJ_cdr3(dico_same_VJ,dicoSeq,Clone_threshold):
	VJ_ID_diff_CDR3 = {}
	dicoclone_vj_cdr3 = {}

	for VJ_ID in dico_same_VJ.keys():
		sub_sub_group = 0
		#print (VJ_ID,"VJ_ID,VJ_ID")
		VJ_ID_diff_CDR3[VJ_ID] = {}
		for seq in dico_same_VJ[VJ_ID]:
			sub_gourp_dist ={}
			CDR3_seq = dicoSeq[seq.rstrip()][3]
			if len(VJ_ID_diff_CDR3[VJ_ID].keys()) != 0 :
				Sub_gourp = VJ_ID_diff_CDR3[VJ_ID].keys()
				for g in Sub_gourp :

					if len(dicoSeq[seq.rstrip()][3]) == len(g) :
						if 1-(hamming_distance(dicoSeq[seq.rstrip()][3] , g)/float(len(g))) >= Clone_threshold :
							sub_gourp_dist[g] =('+',1-(hamming_distance(dicoSeq[seq.rstrip()][3] , g)/float(len(g))))

					elif dicoSeq[seq.rstrip()][2].split("*")[0][-1] == "6" :
						length = max(len(dicoSeq[seq.rstrip()][3]),len(g))
						if 1-(levenshtein_distance(dicoSeq[seq.rstrip()][3] , g)/float(length)) >= Clone_threshold :
							sub_gourp_dist[g] =('+',1-(levenshtein_distance(dicoSeq[seq.rstrip()][3] , g)/float(length)))
				
				if sub_gourp_dist == {}:
					VJ_ID_diff_CDR3[VJ_ID][CDR3_seq] = [seq]
				else :
					dist_loc ={}
					for key in sub_gourp_dist.keys():
						seqs = [skbio.Protein(CDR3_seq, metadata={'id': "CDR3_seq"}),skbio.Protein(key, metadata={'id': "key"})]
						#print (seqs[0],seqs[1])
						msa = skbio.alignment.global_pairwise_align_protein(seqs[0],seqs[1],25)
						dist_loc[key]= float(msa[1])
					#print (dist_loc,"dist_loc")
					best_coressp = (max(dist_loc.items(), key=operator.itemgetter(1))[0])
					VJ_ID_diff_CDR3[VJ_ID][best_coressp].append(seq)
			else :
				VJ_ID_diff_CDR3[VJ_ID][CDR3_seq]= [seq]
	#print (VJ_ID_diff_CDR3)
	return VJ_ID_diff_CDR3

def write_clone_VJ_cdr3(VJ_ID_diff_CDR3,dicoSeq,uniq_seq_dico,output_file,Clone_threshold):
	file_name = output_file+"_sameVJ_noallele_CDR3_"+ Clone_threshold +".txt"
	filetowrite=open(file_name,"w")
	clone_number = 0
	for VJ in VJ_ID_diff_CDR3.keys():
		#print (VJ,"VJ")
		for cdr3 in VJ_ID_diff_CDR3[VJ]:
			#print (cdr3,"cdr3")
			for seq in VJ_ID_diff_CDR3[VJ][cdr3] :
				#print (dicoSeq[seq.rstrip()],"dicoSeq")
				sequence = str(clone_number) + "\t" + str(seq) + "\t" + dicoSeq[seq.rstrip()][0] + "\t" +dicoSeq[seq.rstrip()][1] +"\t" +dicoSeq[seq.rstrip()][2]+ '\t' + dicoSeq[seq.rstrip()][3]+"\n"
				filetowrite.write(sequence)
				if len(uniq_seq_dico[seq]) != 0 :
					for dup in uniq_seq_dico[seq] :
						sequence = str(clone_number) + "\t" + str(dup) + "\t" + dicoSeq[seq.rstrip()][0] + "\t" +dicoSeq[seq.rstrip()][1] +"\t" +dicoSeq[seq.rstrip()][2]+ '\t' +dicoSeq[seq.rstrip()][3]+"\n"
						filetowrite.write(sequence)
			clone_number += 1
	filetowrite.close()
	return 0




def find_max_list(list):
    list_len = [len(i) for i in list]
    return(max(list_len))
#=============================================================================#

def main():
    usage = "python  initial_clustering.py -i <formated IMGT highvquest statistics output> -o <output file name> -s <Clone identity between zero and one> \n "
    parser = OptionParser(usage)
    parser.add_option("-i", "--hv_stat_output", dest="hv_stat_output",
          help="formated IMGT highvquest statistics output")
    parser.add_option("-o", "--output_file_name", dest="output_file_name",
          help="the name for the file to write")
    parser.add_option("-s", "--Clone_threshold", dest="Clone_threshold",
          help="Clone identity percentage between 0 and 1")
    
    (options, args) = parser.parse_args()
    if len(sys.argv) != 7:
        parser.error("incorrect number of arguments")
    IMGT_file = options.hv_stat_output
    output_file_name = options.output_file_name
    Clone_threshold = options.Clone_threshold
    CloneID = read_output_file(IMGT_file)
    print ("all seq : " , len(CloneID))
    annotated_CloneID = filter_seq(CloneID,output_file_name)
    print ("annotated seq : ", len(annotated_CloneID))
    uniq_seq_dico,filtered_lines = delete_duplicate(annotated_CloneID)
    print ("uniq seq : ", len(filtered_lines))
    dico_same_VJ,dicoSeq = group_same_VJ(filtered_lines,output_file_name)
    
    #write_clone_VJ(dico_same_VJ,dicoSeq,output_file_name)
    
    VJ_ID_diff_CDR3 = group_clone_VJ_cdr3(dico_same_VJ,dicoSeq,float(Clone_threshold))
    write_clone_VJ_cdr3(VJ_ID_diff_CDR3,dicoSeq,uniq_seq_dico,output_file_name ,Clone_threshold)

    print ("Done!")

#=============================================================================#

if __name__ == "__main__":
    main()
