#!/bin/bash


mkdir $2
str=$1
IFS='/' # delimiter
read -ra ADDR <<< "$str" # str is read into an array as tokens separated by IFS
name=${ADDR[-1]}
IFS=''

export T=$3

if [[ $T ]]; then     
	printf 'You have set T = %s\n' "$T"; 
else
	export T=0.7
fi


mkdir $2/${name}

python format_labeling_imgt.py -s $1/1_Summary.txt -g $1/2_IMGT-gapped-nt-sequences.txt -o ${name}_seq_Fo.txt
python initial_clustering.py -i ${name}_seq_Fo_V_CDR3_Jseq.txt -o ${name} -s $T
python format_clustering_output.py -i ${name}_sameVJ_noallele_CDR3_0.7.txt -o ${name}_initial_clusters_Fo.txt
python refinement.py -f ${name}_seq_Fo_V_CDR3_Jseq.txt -c ${name}_initial_clusters_Fo.txt
python format_clustering_output.py -i ${name}_seq_Fo_V_CDR3_Jseq_clone_V_CDR3_J.txt -o ${name}_final_clusters_Fo.txt
python two_level_clonal_info.py -f ${name}_seq_Fo_V_CDR3_Jseq.txt -c ${name}_final_clusters_Fo.txt -n ${name}
python Cluster_distribution_plot.py -n ${name} -c ${name}_final_clusters_Fo.txt



#mv $1_repertoire_two_levels_info.txt Output/$1


mv ${name}_final_clusters_seq_info.txt $2/${name}
mv ${name}_unannotated_seq.txt $2/${name}
mv ${name}_initial_clusters_Fo.txt $2/${name}
mv ${name}_final_clusters_Fo.txt $2/${name}
mv ${name}_cluster_distribution.txt $2/${name}
mv ${name}_repertoire.png $2/${name}

rm ${name}_seq_Fo_V_CDR3_Jseq.txt
rm ${name}_sameVJ_noallele_CDR3_0.7.txt
rm ${name}_seq_Fo_V_CDR3_Jseq_clone_V_CDR3_J.txt
rm ${name}_repertoire_two_levels_info.txt

