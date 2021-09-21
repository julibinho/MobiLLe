# MobiLLe

**A multi-objective based clustering for inferring BCR clones from high-throughput B cell repertoire data**

MobiLLe is a new method based on multi-objective clustering to detect clonally-related  sequences in BCR repertoires. It requires V(D)J annotations to obtain the initial clones and iteratively applies two objective functions that optimize cohesion and separation within clones simultaneously. 
MobiLLe can accurately identify clone members, has fewer parameter settings and presents low running time and minimal memory requirements. 
All these features constitute an attractive option for repertoire analysis, particularly in the clinical context for diagnosing and monitoring B cell malignancies.

**REFERENCE**  
Nika Abdollahi, Anne de Septenville, Hugues Ripoche,  Frederic Davi and Juliana S. Bernardes. A multi-objective based clustering for inferring BCR clones from high-throughput B cell repertoire data. 2021 To be submit

**CONTACT**  
  E-mail: 
  juliana.silva_bernardes@sorbonne-universite.fr 
  nikaabdollahi@gmail.com 
  
## Inputs
 
  * The IMGT/HighV-QUEST's output in one unzipped folder.
    The following files must be provided:
    * 1_Summary
    * 2_IMGT-gapped-nt-sequences.txt
  * See [example input files](https://github.com/julibinho/MobiLLe/tree/main/Data/Real_datasets/IMGT_highvquest_output/toy_dataset)
  * You can use any V(D)J annotation software, but inputs shoude be formatted as above.

## Outputs

  * MobiLLe returns:

    - 5 tab delimited files:

      * [repertoire_name]\_unannotated_seq.txt : any sequence that could not be annotated fully [example](https://github.com/julibinho/MobiLLe/tree/main/Data/Real_datasets/MobiLLe_output/toy_dataset/toy_dataset_unannotated_seq.txt)

      The columns are:
      ```
      seq Id   functionality  IGHV_and_allele IGHJ_and_allele CDR3
      ```
      * [repertoire_name]\_cluster_distribution.txt : clusters and their abundance sorted from highest to lowest [example](https://github.com/julibinho/MobiLLe/tree/main/Data/Real_datasets/MobiLLe_output/toy_dataset/toy_dataset_cluster_distribution.txt)

      The columns are:
      ```
      cluster_Id   abundance
      ```

      * [repertoire_name]\_initial_clusters_Fo.txt : initial clustering output. Sequences with the same IGHV and IGHJ genes, same CDR3 sequence length, and CDR3 identity higher than 70% are grouped together [example](https://github.com/julibinho/MobiLLe/tree/main/Data/Real_datasets/MobiLLe_output/toy_dataset/toy_dataset_initial_clusters_Fo.txt)

      Each line contains the cluster id and sequence ids of its members.
      ```
      cluster_Id   seqid1 seqid2 ...
      ```
      * [repertoire_name]\_final_clusters_Fo.txt : final clustering output, after minimizing intraclonal distances and maximizing interclonal distances [example](https://github.com/julibinho/MobiLLe/tree/main/Data/Real_datasets/MobiLLe_output/toy_dataset/toy_dataset_final_clusters_Fo.txt)
      ```
      cluster_Id   seqid1 seqid2 ...
      ```
      * [repertoire_name]\_final_clusters_seq_info.txt : each line contains the following information for each sequence ([example](https://github.com/julibinho/MobiLLe/tree/main/Data/Real_datasets/MobiLLe_output/toy_dataset/toy_dataset_final_clusters_seq_info.txt)):
      ```
      Cluster_id__clonotype_id   seq Id  functionality  IGHV_and_allele IGHJ_and_allele CDR3 Junction
      ```
      

    - A png file containing [example](https://github.com/julibinho/MobiLLe/tree/main/Data/Real_datasets/MobiLLe_output/toy_dataset/toy_dataset_repertoire.png)

      A) Circle representation of the clones' abundances. Each circle symbolizes a clone, and its size represents the clone's abundance.

      B) Number of sequences in each clone, all clones are represented, vertical axe is in log scale.

      C) Lorenz curve and Gini index. A Lorenz curve shows the graphical representation of clonal inequality. On the horizontal axe, it plots the cumulative fraction of total clones when ordered from the less to the most abundant; On the vertical axe, it shows the cumulative fraction of sequences.

      D) Size distribution (percentage) of the 100 most abundant clones.
       
## Requirements 

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
  In the src/ run the following command:
  ```
  $ bash run_MobiLLe.sh [input_repertoire_name] [output_repertoire_name] [options]
  ```
  
  ### required arguments 
  * [input_repertoire_name] is the path directory where are input files, for instance: the IMGT/highVquest's output folder path.
  * [output_newick_file] is the output directory path
   Output files will be placed as such:
  ```
  ~[output_repertoire_name]/[input_repertoire_name]_cluster_distribution.txt
                            [input_repertoire_name]_final_clusters_Fo.txt
                            [input_repertoire_name]_final_clusters_seq_info.txt
                            [input_repertoire_name]_initial_clusters_Fo.txt
                            [input_repertoire_name]_unannotated_seq.txt
                            [input_repertoire_name]_repertoire.png
 ```

  ### optional arguments [...options]

  *  CDR3 amino acid identity threshold (by default 0.7)
                      
 For instance the following command can be run in the src/ folder:
  ```
  $ bash run_MobiLLe.sh  ../Data/Real_datasets/IMGT_highvquest_output/toy_dataset ../Data/Real_datasets/MobiLLe_output/toy_dataset 0.7
  ```


## License, Patches, and Ongoing Developements

  * The program is distributed under the CeCILL licence.  
  * [Feature requests and open issues](https://github.com/julibinho/MobiLLe/issues).

