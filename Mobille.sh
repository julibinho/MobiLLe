#!/usr/bin/env bash

############################
start() {
    echo "##########################################################"
    echo "#                    Mobille                             #"
    echo "##########################################################"

}

stop() {
    echo "Stoping $0"
    exit 1; 
}

help() {
    echo "Usage: $0 -p Path to input -i Analysis name -o Path to output -t 0 -s 0.7 -q 0 -m 0 -r 1 " >&2
        echo " A multi-objective clonal grouping methodto to explore the adaptive immune receptor repertoires clonality"
    echo " 
where :
    -p [string] Path to the file contaning IMGT/High-VQUEST AIRR  file (vquest_airr.tsv)

    -i [string] Analysis name 

    -s [float between 0 and 1] Identity threshold for the initial clustering step 

    -t [int] Abundance filter, the minimum count of sequence to be considered in the analysis, if -t 0, all the sequences will be analysed

    -q [0 (no) or 1 (yes)] Quality filter, if -q 1, sequences contaning N will be discarded from the analysis

    -r [0 (no) or 1 (yes)] Apply refining step  {if -r 0, there will be no need to provide v,j, and c information}

    -v [1 or 2 or 3 or 4] V-distance (1- binaire, 2-levenstein,3- GIANA, 4-K-mers) by default -v 1

    -j [1 or 2 or 3 or 4] J-distance (1- binaire, 2-levenstein,3- GIANA, 4-K-mers) by default -j 2

    -c [1 or 2 or 3 or 4] CDR3-distance (1- binaire, 2-levenstein,3- GIANA, 4-K-mers) by default -c 2

    -d [1 or 2 or 3 or 4] Combine distance (1-mean, 2-  weighted mean(requires a Tuple with weights), 3-harmonic mean)

     example for weighted mean  -d (0.5,0.25,0.25) [first value is for V gene, the second one for the CDR3 and the thired one for the J gene]

    -m [0 (no) or 1 (yes)] Merging the singletons
"

exit 1; 
}



# PARSIN ARGUMENTS 
while getopts "a:b:c:d:e:f:g:hi:j:k:l:m:n:o:p:q:r:s:t:u:v:y:z:w:" OPTION
do
    case $OPTION in
        h)
            help
            ;;
        i)
            AnalysisName=$OPTARG
            ;;
        s)
            SeqIdThre=$OPTARG
            ;;
        t)
            SeqNumThre=$OPTARG
            ;;
        q)
            QualityFilter=$OPTARG
            ;;
        r)
            Refinement=$OPTARG
            ;;
        m)
            MergeSingleton=$OPTARG
            ;;
        d)
            CombineDistance=$OPTARG
            ;;

        v)
            VDistance=$OPTARG
            ;;

        j)
            JDistance=$OPTARG
            ;;

        c)
            CDR3Distance=$OPTARG
            ;;
        p)
            InputPath=$OPTARG
            ;;
        o)
            OutPath=$OPTARG
            ;;
        f)
            ParameterFile=$OPTARG
            ;;
        ?)
            help
            ;;
    esac
done


# Catching no arguments 
if [ $OPTIND -eq 1 ]; then help; fi




# Set the delimiter
IFS='='




if [ -n $Parameters ]; #-n STRING the length of STRING is nonzero
    then
    while read line  || [ -n "$line" ];
        do
            read -a data <<< "$line"
            if [[ ${data[0]} == inputDir* ]]; # * is used for pattern matching
            then
                InputPath=${data[1]};
            elif [[ ${data[0]} == outputDir* ]];
                then
                OutputPath=${data[1]};
            elif [[ ${data[0]} == q* ]];
                then
                QualityFilter=${data[1]};
            elif [[ ${data[0]} == i* ]];
                then
                AnalysisName=${data[1]};

            elif [[ ${data[0]} == s* ]];
                then
                SeqIdThre=${data[1]};

            elif [[ ${data[0]} == t* ]];
                then
                SeqNumThre=${data[1]};
            elif [[ ${data[0]} == d* ]];
                then
                CombineDistance=${data[1]};

            elif [[ ${data[0]} == m* ]];
                then
                MergeSingleton=${data[1]};

            elif [[ ${data[0]} == r* ]];
                then
                Refinement=${data[1]};

            elif [[ ${data[0]} == v* ]];
                then
                VDistance=${data[1]};

            elif [[ ${data[0]} == j* ]];
                then
                JDistance=${data[1]};
            elif [[ ${data[0]} == c* ]];
                then
                CDR3Distance=${data[1]};
            fi
        done < $ParameterFile

fi



# Catching no input path 
if [ -z $InputPath ];
    then 
    	echo ""
        echo "Please provide the path to the the input folder "
        echo "run -h for more information"
        echo ""
        exit 1;
    fi


# if missing parameters, default values are fixed here 

if [ -z $SeqNumThre ];
    then 
        SeqNumThre=0 ;
    fi

if [ -z $QualityFilter ];
    then 
        QualityFilter=1 ;
    fi

if [ -z $SeqIdThre ];
    then 
        SeqIdThre=0.7 ;
    fi

if [ -z $MergeSingleton ];
    then 
        MergeSingleton=0 ;
    fi


if [ -z $AnalysisName ];
    then 
        AnalysisName=MyAnalysis ;
    fi

if [ -z $Refinement ];
    then 
        Refinement=1 ;
    fi

if [ -z $CombineDistance ];
    then 
        CombineDistance=1 ;
    fi


#if the refinement is selected, by default values for the distances are fixed here 

if [ -n $Refinement ]; #-n STRING the length of STRING is nonzero
    then
        echo ""
         if [ -z $VDistance ];
            then 
                VDistance=1 ;
         fi

         if [ -z $JDistance ];
            then 
                JDistance=2 ;
         fi

         if [ -z $CDR3Distance ];
            then 
                CDR3Distance=2 ;
         fi
fi






export BASE_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export SCRIPTS_DIR=$BASE_DIR/Code/src/
export RESULT_DIR=$OutPath

echo ">>> $OutPath"


mkdir $RESULT_DIR
#mkdir ${RESULT_DIR}${AnalysisName}/


#-------------------------------------------------------------

# Create log file
LogFile=$RESULT_DIR/${AnalysisName}_log.txt
# if exist remove
    if [ -e ${LogFile} ]
    then
      rm $LogFile
    fi
# Create log file
touch $LogFile
#fi 

#-------------------------------------------------------------



start
echo "Clustering begins..."
echo $RESULT_DIR/${AnalysisName}_seq_Fo.txt
python ${SCRIPTS_DIR}format_labeling_imgt_airr.py -a ${InputPath}/vquest_airr.tsv -o $RESULT_DIR/${AnalysisName}_seq_Fo -t $SeqNumThre -q $QualityFilter >> $LogFile
python ${SCRIPTS_DIR}initial_clustering_f.py -i $RESULT_DIR/${AnalysisName}_seq_Fo_V_CDR3_Jseq.txt -o $RESULT_DIR/${AnalysisName} -s $SeqIdThre -a $RESULT_DIR/ >> $LogFile
python ${SCRIPTS_DIR}format_clustering_output.py -i $RESULT_DIR/${AnalysisName}_sameVJ_noallele_CDR3_${SeqIdThre}.txt -o $RESULT_DIR/${AnalysisName}_initial_clusters_Fo.txt >> $LogFile

echo "end of clustering"
echo ""


if [[ $Refinement == 1 ]]
	then 
	echo "Refinement begins..."
	echo ""

	python ${SCRIPTS_DIR}refinement.py -f $RESULT_DIR/${AnalysisName}_seq_Fo_V_CDR3_Jseq.txt -l $RESULT_DIR/${AnalysisName}_initial_clusters_Fo.txt -v $VDistance -j $JDistance -c $CDR3Distance -d $CombineDistance -m $MergeSingleton -o $RESULT_DIR/${AnalysisName} >> $LogFile

	echo "end of refinement"
	echo ""
    echo "Formating the output files..."
    echo ""

    python ${SCRIPTS_DIR}format_clustering_output.py -i $RESULT_DIR/${AnalysisName}_clone_V_CDR3_J.txt -o $RESULT_DIR/${AnalysisName}_final_clusters_Fo.txt >> $LogFile
    python ${SCRIPTS_DIR}two_level_clonal_info.py -f $RESULT_DIR/${AnalysisName}_seq_Fo_V_CDR3_Jseq.txt -c $RESULT_DIR/${AnalysisName}_final_clusters_Fo.txt -n $RESULT_DIR/${AnalysisName} >> $LogFile
    python ${SCRIPTS_DIR}uniformity.py -f $RESULT_DIR/${AnalysisName}_seq_Fo_V_CDR3_Jseq.txt -l $RESULT_DIR/${AnalysisName}_final_clusters_Fo.txt -v $VDistance -j $JDistance -c $CDR3Distance -d $CombineDistance -m $MergeSingleton -o $RESULT_DIR/${AnalysisName} >> $LogFile
    python ${SCRIPTS_DIR}Cluster_distribution_plot.py -n ${AnalysisName} -c $RESULT_DIR/${AnalysisName}_final_clusters_Fo.txt -a $RESULT_DIR/ -u $RESULT_DIR/${AnalysisName}_clone_uniformity.txt
else
    echo "Formating the output files..."
    echo ""

    python ${SCRIPTS_DIR}two_level_clonal_info.py -f $RESULT_DIR/${AnalysisName}_seq_Fo_V_CDR3_Jseq.txt -c $RESULT_DIR/${AnalysisName}_initial_clusters_Fo.txt -n $RESULT_DIR/${AnalysisName} >> $LogFile
    python ${SCRIPTS_DIR}uniformity.py -f $RESULT_DIR/${AnalysisName}_seq_Fo_V_CDR3_Jseq.txt -l $RESULT_DIR/${AnalysisName}_initial_clusters_Fo.txt -v 1 -j 2 -c 2 -d 1 -m $MergeSingleton -o $RESULT_DIR/${AnalysisName} >> $LogFile
    python ${SCRIPTS_DIR}Cluster_distribution_plot.py -n ${AnalysisName} -c $RESULT_DIR/${AnalysisName}_initial_clusters_Fo.txt -a $RESULT_DIR/ -u $RESULT_DIR/${AnalysisName}_clone_uniformity.txt

fi




rm $RESULT_DIR/Sorted_clusters.p

echo "The results are ready to explore here : "$RESULT_DIR/
echo "##########################################################"


