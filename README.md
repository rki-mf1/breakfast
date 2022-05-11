# breakfast - FAST outBREAK detection and sequence clustering

[![Tests](https://github.com/rki-mf1/breakfast/workflows/Tests/badge.svg)](https://github.com/rki-mf1/breakfast/actions?workflow=Tests)

`breakfast` is a simple and fast script developed for clustering SARS-CoV-2 genomes using precalculated sequence features (e.g. nucleotide substitutions) from [covSonar](https://gitlab.com/s.fuchs/covsonar) or [Nextclade](https://clades.nextstrain.org/).

**This project is under development and in experimental stage**

<img src="/img/breakfast_logo_2.png" width="300">

## Installation

### Installation using pip

```
$ pip install breakfast
```

### System Dependencies

`breakfast` runs under Python 3.10 and later. The base requirements are networkx, pandas, numpy, scikit-learn, click, and scipy.

### Install using conda

We recommend using conda for installing all necessary dependencies:

```
conda env create -n sonar -f covsonar/sonar.env.yml
conda env create -n breakfast -f breakfast/envs/sc2-breakfast.yml
```

## Example Command Line Usage

### Simple test run
```
conda activate breakfast
breakfast/src/breakfast.py \
   --input-file breakfast/test/testfile.tsv  \
   --max-dist 1 \
   --outdir test-run/
```
You will find your results in `test-run/cluster.tsv`, which should be identical to `breakfast/test/expected_clusters_dist1.tsv`


### 1) covSonar + breakfast
Sequence processing with [covSonar](https://gitlab.com/s.fuchs/covsonar)
```
conda activate sonar
covsonar/sonar.py add -f genomes.fasta --db mydb --cpus 8
covsonar/sonar.py match --tsv --db mydb > genomic_profiles.tsv
```

Clustering with a maximum SNP-distance of 1 and excluding clusters below a size of 5 sequences
```
conda activate breakfast
breakfast/src/breakfast.py \
   --input-file genomic_profiles.tsv \
   --max-dist 1 \
   --min-cluster-size 5 \
   --outdir covsonar-breakfast-results/
```

### 2) Nextclade + breakfast

Sequence processing with [Nextclade CLI](https://clades.nextstrain.org/).

```
conda install -c bioconda nextclade
nextclade dataset get --name 'sars-cov-2' --output-dir 'data/sars-cov-2'
nextclade \
   --in-order \
   --input-fasta genomes.fasta \
   --input-dataset data/sars-cov-2 \
   --output-tsv output/nextclade.tsv \
   --output-tree output/nextclade.auspice.json \
   --output-dir output/ \
   --output-basename nextclade
```

Alternatively, you can also use [Nextclade Web](https://clades.nextstrain.org/) to process your fasta and export the genomic profile as "nextclade.tsv".

Clustering with a maximum SNP-distance of 1 and excluding clusters below a size of 5 sequences. Since the input tsv of Nextclade looks a little different from the covSonar tsv, you need to specify the additional parameters `--id-col`, `--clust-col` and `--sep2` for identifying the correct columns.

```
conda activate breakfast
breakfast/src/breakfast.py \
   --input-file output/nextclade.tsv \
   --max-dist 1 \
   --min-cluster-size 5 \
   --id-col "seqName" \
   --clust-col "substitutions" \
   --sep2 "," \
   --outdir nextclade-breakfast-results/
```

## Parameter description

| Parameter              | Type    	| Required | Default 	| Description                                |
|----------------------- |---------	|----------|----------|------------------------------------------  |
| --input-file           | String     	|✅	     | 'genomic_profiles.tsv.gz'    	| Path of the input file (in tsv format)     |
| --max-dist              | Integer  	|	     | 1     	| Two sequences will be grouped together, if their pairwise edit distance does not exceed this threshold |
| --min-cluster-size  | Integer  	|      | 2     	| Minimum number of sequences a cluster needs to include to be defined in the result file      |
| --id-col    | String 	|     | 'accession'      	| Name of the sequence identifier column of the input file          |
| --clust-col              | String 	|     | 'dna_profile'      | Name of the mutation profile column of the input file         |
| --var-type              | String 	|     | 'dna'       | Specify if DNA or AA substitutions are used for the mutation profiles         |
| --sep              | String 	|     | '\t'      | Input file separator       |
|  --sep2              | String 	|     | '  '      | Secondary clustering column separator (between each mutation)        |
| --outdir              | String 	|     | 'output/'       | Path of output directory        |
| --trim-start               | Integer 	|     |264       | Bases to trim from the beginning         |
| --trim-end               | Integer 	|     | 228       | Bases to trim from the end         |
| --reference-length              | Integer 	|     | 29903      | Length of reference genome (defaults to NC_045512.2)        |
| --skip-del               | Bool 	|     | TRUE       | Deletions will be skipped for calculating the pairwise distance of your input sequences.|
| --skip-ins               | Bool 	|     | TRUE       | Insertions will be skipped for calculating the pairwise distance of your input sequences.         |
| --input-cache           | Integer 	|     | None   | Path to import results from previous run |
| --output-cache              | String 	|     | None       | Path to export results which can be used in the next run to decrease runtime.  |
| --help                   | N/A     	|	   | N/A     	| Show this help message and exit            |
| --version                | N/A     	|	   | N/A     	| Show version and exit            |
