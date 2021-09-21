
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
def dico_V_CDR3_format(Summary_file):
	lines = read_output_file(Summary_file)
	Dico={}
	for l in range(1,len(lines)):
		functionality,V,J,CDR3 = "_","_","_","_"
		split=lines[l].split("\t")
		sequence_id = split[1]
		if len(split)>5:
			if split[2] != "":
				functionality = split[2].split(" ")[0]
			if split[3] != "":
				V = split[3].split(" ")[1]
				#print (V)
			if split[9] != "":
				J = split[9].split(" ")[1]
				#print (J)
			if split[20] != "":
				CDR3 = split[20].split(" ")[0].replace("#", ".")
				#print (CDR3)
		Dico[sequence_id] = [functionality,V,J,CDR3]
	return Dico
#####################################################################
def dico_J_seq_format(Dico_from_summary, gapped_nt_file ):
	Dico_VJCDR3 = Dico_from_summary
	lines = read_output_file(gapped_nt_file)
	#print(lines[0].split("\t")[16])
	for l in range(1,len(lines)):
		split=lines[l].split("\t")
		sequence_id = split[1]
		if len(split)>5:
			#print (split[17])
			Dico_VJCDR3[sequence_id].append(split[16])
		else :
			Dico_VJCDR3[sequence_id].append("_")
	#print (Dico_VJCDR3)
	return Dico_VJCDR3
#####################################################################				
def write_file(Dico_VJCDR3,output_file):
	outputname = output_file.split(".")[0]+"_V_CDR3_Jseq.txt"
	f = open(outputname,"w")
	for key in Dico_VJCDR3.keys():
		#print (Dico_VJCDR3[key])
		line = key + "\t" + Dico_VJCDR3[key][0] + "\t" + Dico_VJCDR3[key][1] + "\t" + Dico_VJCDR3[key][2] + "\t" + Dico_VJCDR3[key][3] + "\t" + Dico_VJCDR3[key][4] + "\n"
		f.write(line)
	f.close()
	return 0


####################################################################
def main():
	usage = "usage: format_labeling_imgt.py -s Summary_file -g gapped_nt_file -o output_file"
	parser = OptionParser(usage)
	parser.add_option("-s", "--Summary_file", dest="Summary_file",
	      help="read data from 1_Summary.txt")
	parser.add_option("-g", "--gapped_nt_file",dest="gapped_nt_file",
	      help="read data from 2_IMGT-gapped-nt-sequences.txt")
	parser.add_option("-o", "--output_file",dest="output_file",
	      help="write data to output_file")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 7:
		parser.error("incorrect number of arguments")
	
	Summary_file = options.Summary_file
	gapped_nt_file = options.gapped_nt_file
	output_file = options.output_file
	time_start = time.perf_counter()

	Dico_from_summary = dico_V_CDR3_format(Summary_file)
	Dico_VJCDR3 = dico_J_seq_format(Dico_from_summary, gapped_nt_file )
	write_file(Dico_VJCDR3,output_file)


#####################################################################
if __name__ == "__main__":
	main()
