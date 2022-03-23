# breakfast - FAST outBREAK detection and sequence clustering

`breakfast` is a simple and fast script developed for clustering SARS-CoV-2 genomes using precalculated sequence features (e.g. nucleotide substitutions) from [covSonar](https://gitlab.com/s.fuchs/covsonar). 

**This project is under development and in experimental stage**

## Installation

### System Dependencies
Breakfast runs under Python 3.9 and later. The base requirements are networkx, pandas, numpy, scikit-learn, and scipy. 

### Install using conda
We recommend using conda for installing all necessary dependencies:
```
conda env create -n breakfast -f breakfast/envs/sc2-breakfast.yml
```
## Example Commandline Usage
Clustering with a maximum SNP-distance of 1 and excluding clusters below a size of 5 sequences.

```
./src/breakfast.py --input-file INPUT.tsv.gz --max-dist 1 --min-cluster-size 5
```

## Parameter description

#### --input-file
A TSV file with mutation profiles to be clustered. The file requires a column specifiying the sequence ID and a column specifying the DNA or AA mutation profile.
We recommend using [covSonar](https://gitlab.com/s.fuchs/covsonar) to generate mutation profiles from a fasta file. 

#### --id-col
The name of the sequence identifier column of the input file. If none is given, it will be 'accession' by default.

#### --clust-col
The name of the mutation profile column of the input file. If none is given, it will be 'dna_profile' by default.

#### --var-type
Specify if DNA or AA substitutions are used for the mutation profiles. By default DNA.

#### --sep, --sep2
Specify the clustering column seperator (between each mutation) and the input file separator. Will be tab-separated by default.

#### --outdir
Specify the output directory for the output file. The output file will be a tsv file including a sequence identifier column and a column for specifiying the cluster ID of each sequence. If no output directory is provided, a folder named 'output/' will be created in the current directory.
  
#### --max-dist   
Two sequences will be grouped together, if their pairwise edit distance does not exceed this threshold. If none is given, a pairwise edit distance of 1 will be used for clustering sequences.
  
#### --min-cluster-size 
Specify the minimum number of sequences a cluster needs to have to be included in the result file. Clusters below this threshold will not be shown in the result file. If none is given, single-sequence clusters will be excluded by default.

#### --trim-start 
Specify the bases to trim from the beginning. (default: 264)
  
#### --trim-end    
Specify the bases to trim from the end (default: 228)
  
#### --reference-length 
Specify the length of reference genome (defaults to NC_045512.2 length = 29903)
  
#### --skip-del, --no-skip-del
By default deletions will be skipped for calculating the pairwise distance of your input sequences.

#### --skip-ins, --no-skip-ins
By default insertions will be skipped for calculating the pairwise distance of your input sequences.

#### --input-cache
Specify path to import results from previous run as pickle file. Only new sequences will be used for distance matrix computation to reduce the runtime. If none is given, the complete dataset will be used for computing the pairwise distances via sparse matrix computation.

#### --output-cache
Specify path to export results as pickle file. Results can be saved and used in the next run to reduce the runtime. 
