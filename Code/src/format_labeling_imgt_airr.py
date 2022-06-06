
import sys
from collections import Counter
from optparse import OptionParser
import time

#####################################################################
def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines


#####################################################################
def filter_dico(Dico,remain_seq) :
	list_all_seq = Dico.keys()
	to_delet = set(list_all_seq) - set(remain_seq)
	for seq in to_delet:
		del Dico[seq]
	return Dico, len(to_delet)

#####################################################################
def dico_V_J_CDR3_format(AIRR,abundance_threshold,quality_filter):

	lines = read_output_file(AIRR)
	low_quality_seq, remain_seq =[], []
	Dico={}
	Dico_uniq_adundance = {}
	for l in range(1,len(lines)):
		split=lines[l].split("\t")
		#split[13] is the sequence alignment, is better than the sequence because it starts at the beginig of the V and ends at the end of J 
		
		if quality_filter == "True" or quality_filter == "TRUE" :
			#sequences contaning N are discarded 
			if 'n' in split[13].lower(): 
				low_quality_seq.append(l)
			else :
				functionality,V,J,CDR3,Jseq,Vseq,Jseq_aa,Vseq_aa = "_","_","_","_","_","_","_","_"
				sequence_id = split[0]
				if split[4] != "":
					functionality = split[4]
				if split[9] != "":
					if len(split[9].split(" "))>1 :
						V = split[9].split(" ")[1]
					else: 
						V = split[9]
				
				if split[11] != "":
					if len(split[11].split(" "))>1 :
						J = split[11].split(" ")[1]
					else :
						J = split[11]
				
				if split[28] != "":
					CDR3 = split[28].replace("#", ".")
				if split[89] != "":
					Jseq = split[89].replace('.', '')
				if split[85] != "":
					Vseq = split[85].replace('.', '')
				if split[90] != "":
					Jseq_aa = split[90].replace('.', '')
				if split[86] != "":
					Vseq_aa = split[86].replace('.', '')
				Dico[sequence_id] = [functionality,V,J,CDR3,Jseq,Vseq,Jseq_aa,Vseq_aa ]
				# keep the sequences based on the read abundace
				if split[13] in Dico_uniq_adundance.keys() :
					Dico_uniq_adundance[split[13]].append(sequence_id)
				else : 
					Dico_uniq_adundance[split[13]] = [sequence_id]

		else:
			functionality,V,J,CDR3,Jseq,Vseq,Jseq_aa,Vseq_aa = "_","_","_","_","_","_","_","_"
			sequence_id = split[0]
			if split[4] != "":
				functionality = split[4]
			if split[9] != "":
				if len(split[9].split(" "))>1 :
					V = split[9].split(" ")[1]
				else: 
					V = split[9]
			
			if split[11] != "":
				if len(split[11].split(" "))>1 :
					J = split[11].split(" ")[1]
				else :
					J = split[11]
			
			if split[28] != "":
				CDR3 = split[28].replace("#", ".")
			if split[89] != "":
				Jseq = split[89].replace('.', '')
			if split[85] != "":
				Vseq = split[85].replace('.', '')
			if split[90] != "":
				Jseq_aa = split[90].replace('.', '')
			if split[86] != "":
				Vseq_aa = split[86].replace('.', '')
			Dico[sequence_id] = [functionality,V,J,CDR3,Jseq,Vseq,Jseq_aa,Vseq_aa ]
			# keep the sequences based on the read abundace
			if split[13] in Dico_uniq_adundance.keys() :
				Dico_uniq_adundance[split[13]].append(sequence_id)
			else : 
				Dico_uniq_adundance[split[13]] = [sequence_id]

	for key in Dico_uniq_adundance.keys():
		if len(Dico_uniq_adundance[key]) > abundance_threshold : 
			remain_seq.append(Dico_uniq_adundance[key])
	flat_list_remain_seq = [item for sublist in remain_seq for item in sublist]
	#delete low abundant unique sequences 
	final_dico, deleted = filter_dico(Dico,flat_list_remain_seq)
	number_of_analyzed_seq = len(Dico.keys()) - deleted
	print("Total sequence count : " , len(lines)-1)
	print("Number of low quality sequences that have been eliminated from the analysis (N in the sequence): ",  len(low_quality_seq))
	print("Number of low abundant sequences that have been eliminated from the analysis, using the threshold of "+ str(abundance_threshold)+": ",  deleted)
	print("Number of sequences to be analysed : ",  (len(lines)-1)-len(low_quality_seq)-deleted)
	return final_dico

#####################################################################				
def write_file(Dico_VJCDR3,output_file):
	outputname = output_file+"_V_CDR3_Jseq.txt"

	#the output file contains the following information for each sequence 
	#sequence_id	productive	v_call	j_call	cdr3_aa	j_sequence_alignment	v_sequence_alignment	j_sequence_alignment_aa	v_sequence_alignment_aa
	#the gaps in all sequences have been removed 

	f = open(outputname,"w")
	for key in Dico_VJCDR3.keys():
		line = key + "\t" + Dico_VJCDR3[key][0] + "\t" + Dico_VJCDR3[key][1] + "\t" + Dico_VJCDR3[key][2] + "\t" + Dico_VJCDR3[key][3] + "\t" + Dico_VJCDR3[key][4]+"\t" + Dico_VJCDR3[key][5]+"\t"+ Dico_VJCDR3[key][6]+"\t" + Dico_VJCDR3[key][7]+"\n"
		f.write(line)
	f.close()
	return 0


####################################################################
def main():
	start_time = time.time()
	usage = "usage: format_labeling_imgt_airr.py -a AIRR -o output_file -t abundance_threshold -q quality_filter"
	parser = OptionParser(usage)
	parser.add_option("-a", "--AIRR", dest="AIRR",
	      help="read data from AIRR")
	parser.add_option("-o", "--output_file",dest="output_file",
	      help="write data to output_file")
	parser.add_option("-t", "--abundance_threshold",dest="abundance_threshold",
	      help="the abundance threshold of unique sequences to be analyzed, default = 3")
	parser.add_option("-q", "--quality_filter",dest="quality_filter",
	      help="If True, sequences contaning N are discarded from the analysis")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 9:
		parser.error("incorrect number of arguments")
	
	AIRR_file = options.AIRR
	output_file = options.output_file
	abundance_threshold = options.abundance_threshold
	quality_filter = options.quality_filter 

	time_start = time.perf_counter()
	Dico_VJCDR3 = dico_V_J_CDR3_format(AIRR_file,int(abundance_threshold),quality_filter)
	write_file(Dico_VJCDR3,output_file)
	print("The AIRR file reading step execution time : %s seconds " % (time.time() - start_time))


#####################################################################
if __name__ == "__main__":
	main()
