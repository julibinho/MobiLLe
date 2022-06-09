# MobiLLe

**A multi-objective based clustering for inferring BCR clonal lineages from high-throughput B-cell repertoire data**

MobiLLe is a new method based on multi-objective clustering to detect clonally-related  sequences in BCR repertoires. It requires V(D)J annotations to obtain the initial clones and iteratively applies two objective functions that optimize cohesion and separation within clones simultaneously. 
MobiLLe can accurately identify clone members, has fewer parameter settings and presents low running time and minimal memory requirements. 
All these features constitute an attractive option for repertoire analysis, particularly in the clinical context for diagnosing and monitoring B cell malignancies.

**REFERENCE**  
Nika Abdollahi, Lucile Jeusset, Anne de Septenville, Hugues Ripoche,  Frederic Davi and Juliana Silva Bernardes. A multi-objective based clustering for inferring BCR clonal lineages from high-throughput B cell repertoire data. 2022 (under revision)

**CONTACT**  
  E-mail: 
  juliana.silva_bernardes@sorbonne-universite.fr 
  nikaabdollahi@gmail.com 
  
## Inputs
 
  * The IMGT/HighV-QUEST's AIRR file must be provided:
    * vquest_airr.tsv
  * See [example input file](https://github.com/julibinho/MobiLLe/blob/main/Input/toy_dataset/vquest_airr.tsv)
  * You can use any V(D)J annotation software, but input should be formatted as above.

## Outputs

  * MobiLLe returns:

    - 9 tab delimited files:

      * [repertoire_name]\_repertoire_two_levels_info.txt : each lines contains informations of a clone [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_repertoire_two_levels_info.txt)

      The columns are:
      ```
      Clone number   abundance   number of reads   Clonotype abundance, functionality   IGHV_and_allele   IGHJ_and_allele   Clonotypes CDR3, sequence_id
      ```
      * [repertoire_name]\_cluster_distribution.txt : clusters and their abundance sorted from highest to lowest [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_cluster_distribution.txt)

      The columns are:
      ```
      cluster_Id   abundance
      ```

      * [repertoire_name]\_initial_clusters_Fo.txt : initial clustering output. Sequences with the same IGHV and IGHJ genes, same CDR3 sequence length, and CDR3 identity higher than 70% are grouped together [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_initial_clusters_Fo.txt)

      Each line contains the cluster id and sequence ids of its members.
      ```
      cluster_Id   seqid1 seqid2 ...
      ```
      * [repertoire_name]\_final_clusters_Fo.txt : final clustering output, after minimizing intraclonal distances and maximizing interclonal distances [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_final_clusters_Fo.txt)
      ```
      cluster_Id   seqid1 seqid2 ...
      ```
      * [repertoire_name]\_clusters_seq_info.txt : each line contains the following information for each sequence [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_clusters_seq_info.txt):
      ```
      Cluster_id__clonotype_id   seq Id  functionality  IGHV_and_allele IGHJ_and_allele CDR3 Junction
      ```
      * [repertoire_name]\_clone_uniformity.txt : uniformity of each clusters [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_clone_uniformity.txt)
      ```
      cluster_Id   uniformity
      ```
      * [repertoire_name]\_clone_V_CDR3_J.txt [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_clone_V_CDR3_J.txt)
      * [repertoire_name]\_sameVJ_noallele_CDR3_0.7.txt [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_sameVJ_noallele_CDR3_0.7.txt)
      * [repertoire_name]\_seq_Fo_V_CDR3_Jseq.txt [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_seq_Fo_V_CDR3_Jseq.txt) 
      These 3 files contain informations of the sequences of each clusters.
      

    - A png file containing [example](https://github.com/julibinho/MobiLLe/blob/main/Output/toy_dataset/toy_dataset_repertoire.png)

      A) Circle representation of the clones' uniformity and abundances. Each circle symbolizes a clone and for the 20 most abundant clones, its size represents the clone's abundance.

      B) Number of sequences in each clone, all clones are represented, vertical axe is in log scale.

      C) Lorenz curve and Gini index. A Lorenz curve shows the graphical representation of clonal inequality. On the horizontal axe, it plots the cumulative fraction of total clones when ordered from the less to the most abundant; On the vertical axe, it shows the cumulative fraction of sequences.

      D) Size distribution (percentage) of the 100 most abundant clones.
       
## Requirements 

  * We strongly recommend [anaconda](https://docs.anaconda.com/anaconda/install/) environment. 
  
  * Python version 3 or later

  * numpy :
      ```
      conda install numpy
      ```
      or 
      ```
      pip install numpy
      ```

  * matplotlib
    ```
      conda install -c conda-forge matplotlib
     ```
     or
  
      ```
      pip install matplotlib
      ```
      
  * Palettable :
      ```
      conda install -c conda-forge palettable
      ```
      or
      ```
      pip install palettable
      ```

  * skbio
      ```
      conda install -c https://conda.anaconda.org/biocore scikit-bio
      ```
      or 
      ```
      pip install scikit-bio
      ```
  * Levenshtein
      ```
      conda install -c conda-forge python-levenshtein 
      ```
      or
      ```
      pip install python-Levenshtein
      ```


## Using MobiLLe 
  In the MobiLLe file run the following command:
  ```
  $ bash Mobille.sh -p [input_repertoire_name] -o [output_repertoire_name] -i [analysis_name] [options]
  ```
  or you can pass a file that contains your parameters ([example](https://github.com/julibinho/MobiLLe/blob/main/parameter.txt)) :
  ```
  $ bash Mobille.sh -f [parameters_file] 
  ```

  ### required arguments 
  * [input_repertoire_name] is the path directory where are input file, for instance: the IMGT/highVquest's output folder path.
  * [output_newick_file] is the output directory path
   Output files will be placed as such:
  ```
  ~[output_repertoire_name]/[analysis_name]_cluster_distribution.txt
                            [analysis_name ]_final_clusters_Fo.txt
                            [analysis_name]_final_clusters_seq_info.txt
                            [analysis_name ]_initial_clusters_Fo.txt
                            [analysis_name]_unannotated_seq.txt
                            [analysis_name]_repertoire.png
 ```

  ### optional arguments [...options]

*  s : CDR3 amino acid identity threshold (by default 0.7) for the initial clustering step (between 0 and 1)
  *  t : Abundance filter, the minimum count of sequence to be considered in the analysis, if -t 0, all the sequences will be analysed
  *  q : Quality filter, if -q 1, sequences contaning N will be discarded from the analysis (0 : no, 1 :yes)
  *  r : Apply refining step (0 : no , 1 :yes). If -r 0, there will be no need to provide v, j, c, and m parameters.
  *  v : V-distance (1- binaire, 2-levenstein,3- GIANA, 4-K-mers) by default -v 1
  *  j : J-distance (1- binaire, 2-levenstein,3- GIANA, 4-K-mers) by default -j 2
  *  c : CDR3-distance (1- binaire, 2-levenstein,3- GIANA, 4-K-mers) by default -c 2
  *  d : Combine distance (1-mean, Tuple with three weights for  (IGHV, CDR3, IGHJ)). 
  *  m : Merging singletons (0 : no, 1 :yes)
                      
 Exemple : 
  ```
  $ bash Mobille.sh -p Input/toy_dataset/ -o Output/toy_dataset/ -i toy_dataset -t 0 -s 0.8 -q 1 -m 0 -r 1 -v 2 -j 2 -c 2 -d 123
  ```


## License, Patches, and Ongoing Developements

  * The program is distributed under the CeCILL licence.  
  * [Feature requests and open issues](https://github.com/julibinho/MobiLLe/issues).
