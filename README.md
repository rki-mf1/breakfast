# BREAKFAST - FAST outBREAK detection and sequence clustering

`BREAKFAST` is a simple and fast script developed for clustering SARS-CoV-2 genomes using precalculated sequence features (e.g. nucleotide substitutions) from [covSonar](https://gitlab.com/s.fuchs/covsonar). 

**This project is under development and in experimental stage**

<img src="/img/breakfast_logo_2.png" width="300">

## Installation

### System Dependencies
BREAKFAST runs under Python 3.9 and later. The base requirements are networkx, pandas, numpy, scikit-learn, and scipy. 

### Install using conda
We recommend using conda for installing all necessary dependencies:
```
conda env create -n sonar -f covsonar/sonar.env.yml
conda env create -n breakfast -f breakfast/envs/sc2-breakfast.yml
```

## Example Commandline Usage
Sequence processing with [covSonar](https://gitlab.com/s.fuchs/covsonar)

```
conda activate sonar
covsonar/sonar.py add -f genomes.fasta --db mydb --cpus 8
covsonar/sonar.py match --tsv --db mydb > genomic_profiles.tsv
```

Clustering with a maximum SNP-distance of 1 and excluding clusters below a size of 5 sequences

```
conda activate breakfast
breakfast/src/breakfast.py --input-file genomic_profiles.tsv --max-dist 1 --min-cluster-size 5
```

## Parameter description

| Parameter              | Type    	| Required | Default 	| Description                                |           
|----------------------- |---------	|----------|----------|------------------------------------------  |
| --input-file           | String     	|TRUE	     | '../input/covsonar/rki-2021-05-19-minimal.tsv.gz'    	| Path of the input file (in tsv format)     |
| --max-dist              | Integer  	|FALSE	     | 1     	| Two sequences will be grouped together, if their pairwise edit distance does not exceed this threshold |
| --min-cluster-size  | Integer  	|FALSE      | 2     	| Minimum number of sequences a cluster needs to include to be defined in the result file      |
| --id-col    | String 	|FALSE     | 'accession'      	| Name of the sequence identifier column of the input file          |
| --clust-col              | String 	|FALSE     | 'dna_profile'      | Name of the mutation profile column of the input file         |
| --var-type              | String 	|FALSE     | 'DNA'       | Specify if DNA or AA substitutions are used for the mutation profiles         |
| --sep              | String 	|FALSE     | '\t'      | Input file separator       |
|  --sep2              | String 	|FALSE     | '  '      | Secondary clustering column separator (between each mutation)        |
| --outdir              | String 	|FALSE     | 'output/'       | Path of output directory        |
| --trim-start               | Integer 	|FALSE     |264       | Bases to trim from the beginning         |
| --trim-end               | Integer 	|FALSE     | 228       | Bases to trim from the end         |
| --reference-length              | Integer 	|FALSE     | 29903      | Length of reference genome (defaults to NC_045512.2)        |
| --skip-del               | Bool 	|FALSE     | TRUE       | Deletions will be skipped for calculating the pairwise distance of your input sequences.|
| --skip-ins               | Bool 	|FALSE     | TRUE       | Insertions will be skipped for calculating the pairwise distance of your input sequences.         |
| --input-cache           | Integer 	|FALSE     | None   | Path to import results from previous run as pickle file. Only new sequences will be used for distance matrix computation to reduce the runtime. If none is given, the complete dataset will be used for computing the pairwise distances via sparse matrix computation        |
| --output-cache              | String 	|FALSE     | None       | Path to export results as pickle file. Results can be saved and used in the next run to reduce the runtime.       |
| --help                   | N/A     	|FALSE	   | N/A     	| Show this help message and exit            |
| --version                | N/A     	|FALSE	   | N/A     	| Show version and exit            |


