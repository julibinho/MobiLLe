# MobiLLe

**A multi-objective based clustering for inferring BCR clones from high-throughput B cell repertoire data **

MobiLLe is a new method based on multi-objective clustering to detect clonally-related  sequences in BCR repertoires. It Our requires V(D)J annotations to obtain the initial clones and iteratively applies two objective functions that optimize cohesion and separation within clones simultaneously. 
MobiLLe can accurately identify clone members, has fewer parameter settings and presents low running time and  minimal memory requirements. 
All these features constitute an attractive option for repertoire analysis, particularly in the clinical context for diagnosing and monitoring B cell malignancies.

**REFERENCE**  


**CONTACT**  
  E-mail: 
  juliana.silva_bernardes@sorbonne-universite.fr 
  nikaabdollahi@gmail.com 
  
## Inputs
 
  * The IMGT/HighV-QUEST's output in one unzipped folder.
    The following files must be provided:
    * 1_Summary
    * 2_IMGT-gapped
  * See [example input files](https://github.com/NikaAb/BCR_GTG/tree/master/Data/Real_datasets/IMGT_highvquest_output/toy_dataset)

## Outputs

  * MobiLLe returns:

    - 5 tab delimited file:

      * [repertoire_name]\_unannotated_seq.txt : any sequences that could not be annotated fully [example](https://github.com/NikaAb/BCR_GTG/blob/master/Data/GTM_output/I1_IMGT/I1_IMGT_unannotated_seq.txt)

      The columns are:
      ```
      seq Id   functionality  IGHV_and_allele IGHJ_and_allele CDR3
      ```
      * [repertoire_name]\_cluster_distribution.txt : clusters and their abundance sorted from highest to lowest [example](https://github.com/NikaAb/BCR_GTG/blob/master/Data/GTM_output/I1_IMGT/I1_IMGT_cluster_distribution.txt)

      The columns are:
      ```
      cluster_Id   abundance
      ```

      * [repertoire_name]\_initial_clusters_Fo.txt : initial clustering output. Sequences with the same IGHV and IGHJ genes, same CDR3 sequence length, and CDR3 identity higher than 70% are grouped together [example](https://github.com/NikaAb/BCR_GTG/blob/master/Data/GTM_output/I1_IMGT/I1_IMGT_initial_clusters_Fo.txt)

      Each line contains the one cluster id and all the sequence ids of it's members.
      ```
      cluster_Id   seqid1 seqid2 ...
      ```
      * [repertoire_name]\_final_clusters_Fo.txt : final clustering output, after minimizing intraclonal distances and maximizing interclonal distances [example](https://github.com/NikaAb/BCR_GTG/blob/master/Data/GTM_output/I1_IMGT/I1_IMGT_final_clusters_Fo.txt)
      ```
      cluster_Id   seqid1 seqid2 ...
      ```
      * [repertoire_name]\_final_clusters_seq_info.txt : each line contains the following information for each sequence ([example](https://github.com/NikaAb/BCR_GTG/blob/master/Data/GTM_output/I1_IMGT/I1_IMGT_final_clusters_seq_info.txt)):
      ```
      Cluster_id__clonotype_id   seq Id  functionality  IGHV_and_allele IGHJ_and_allele CDR3 Junction
      ```
      

    - A png file containing:

      ![alt text](https://github.com/NikaAb/BCR_GTG/blob/master/Data/GTM_output/I1_IMGT/I1_IMGT_repertoire.png )

      A) Circle representation of the clone abundance. Each  circle  symbolizes  a  clone, and the clone’s abundance is represented by its size.

      B) Number of sequences in each clone, all clones are represented, vertical axe is in log scale.

      C) Lorenz curve and Gini index. A Lorenz curve shows the graphical represen-tation of clonal inequality. On the horizontal axe, it plots the cumulative fraction of total clones when ordered from the less to the most abundant; On the vertical axe, it shows the cumulative fraction of sequences.

      D) Percentage of the 100 most abundant clones.
       
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
  $ bash run_MobiLLe.sh [input_repertoire_name] [output_repertoire_name]
  ```
                      
  Output files will be placed as such:
  ```
  ~[output_repertoire_name]/[input_repertoire_name]_cluster_distribution.txt
                            [input_repertoire_name]_final_clusters_Fo.txt
                            [input_repertoire_name]_final_clusters_seq_info.txt
                            [input_repertoire_name]_initial_clusters_Fo.txt
                            [input_repertoire_name]_unannotated_seq.txt
                            [input_repertoire_name]_repertoire.png
 ```
 [input_repertoire_name] is the IMGT/highVquast's output folder name.
## License, Patches, and Ongoing Developements

  * The program is distributed under the .  
  * [Feature requests and open issues](https://github.com/NikaAb/BCR_GTG/issues).

