# breakfast - FAST outBREAK detection and sequence clustering

[![Tests](https://github.com/rki-mf1/breakfast/workflows/Tests/badge.svg)](https://github.com/rki-mf1/breakfast/actions?workflow=Tests)

`breakfast` is a simple and fast script developed for clustering SARS-CoV-2 genomes using precalculated sequence features (e.g. nucleotide substitutions) from [covSonar](https://github.com/rki-mf1/covsonar) or [Nextclade](https://clades.nextstrain.org/).

**This project is under development and in experimental stage**

<img src="/img/breakfast_logo_2.png" width="300">

## Installation

### Installation using conda/mamba

Breakfast is available in [bioconda](http://bioconda.github.io/recipes/breakfast/README.html). You can install it using either the conda command, or if you've installed [mamba](https://github.com/mamba-org/mamba) you can use that:

```
$ conda install breakfast
# or
$ mamba install breakfast
```

### Installation using pip

Conda is available from [PyPI](https://pypi.org/project/breakfast/) and can be installed using pip:

```
$ pip install breakfast
```

## Example Command Line Usage

### Simple test run

```
breakfast \
   --input-file breakfast/test/testfile.tsv  \
   --max-dist 1 \
   --outdir test-run/
```

You will find your results in `test-run/cluster.tsv`, which should be identical to `breakfast/test/expected_clusters_dist1.tsv`

### Using breakfast with input from covsonar

breakfast uses pre-calculated sequence features (= mutations) as input rather than raw sequences. These features can be calculated with several different programs, but the one we mainly use is [covSonar](https://github.com/rki-mf1/covsonar). It can be used to maintain a database of mutations for a large number of sequences, which can then be easily queried.

```
conda activate sonar
covsonar/sonar.py add -f genomes.fasta --db mydb --cpus 8
covsonar/sonar.py match --tsv --db mydb > genomic_profiles.tsv
```

Clustering with a maximum SNP-distance of 1 and excluding clusters below a size of 5 sequences:

```
breakfast \
   --input-file genomic_profiles.tsv \
   --max-dist 1 \
   --min-cluster-size 5 \
   --outdir covsonar-breakfast-results/
```

### Using breakfast with input from Nextclade

An alternative to covsonar that is commonly used is [Nextclade CLI](https://clades.nextstrain.org/).

```
conda install -c bioconda nextclade  # If nextclade isn't already installed
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
breakfast \
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

## Dependencies

`breakfast` runs under Python 3.10 and later. We rely heavily on some excellent open source python libraries: networkx, pandas, numpy, scikit-learn, click, and scipy.
