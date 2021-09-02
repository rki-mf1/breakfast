# breakfast - Fast outbreak detection and sequence clustering

`breakfast` is a simple and fast clustering script developed for SARS-CoV-2
genomes and using precalculated sequence features (e.g. nucleotide
substitutions). 

## Usage

```
usage: breakfast.py [-h] [--input-file INPUT_FILE] [--id-col ID_COL] [--clust-col CLUST_COL] [--var-type VAR_TYPE]
                    [--sep2 SEP2] [--sep SEP] [--outdir OUTDIR] [--max-dist MAX_DIST] [--min-cluster-size MIN_CLUSTER_SIZE]
		    [--trim-start TRIM_START] [--trim-end TRIM_END] [--reference-length REFERENCE_LENGTH]
		    [--skip-del | --no-skip-del] [--skip-ins | --no-skip-ins]

optional arguments:
  -h, --help            show this help message and exit
  --input-file INPUT_FILE
                        Input file (default: ../input/covsonar/rki-2021-05-19-minimal.tsv.gz)
  --id-col ID_COL       Column with the sequence identifier (default: accession)
  --clust-col CLUST_COL
                        Metadata column to cluster (default = 'dna_profile') (default: dna_profile)
  --var-type VAR_TYPE   Type of variants (dna or aa, default = 'dna') (default: dna)
  --sep2 SEP2           Secondary clustering column separator (between each mutation) (default: )
  --sep SEP             Input file separator (default: )
  --outdir OUTDIR       Output directory for all output files (default: output)
  --max-dist MAX_DIST   Maximum parwise distance (default: 1)
  --min-cluster-size MIN_CLUSTER_SIZE
                        Minimum cluster size (default: 2)
  --trim-start TRIM_START
                        Bases to trim from the beginning (0 = disable) (default: 264)
  --trim-end TRIM_END   Bases to trim from the end (0 = disable) (default: 228)
  --reference-length REFERENCE_LENGTH
                        Length of reference genome (defaults to NC_045512.2 length = 29903) (default: 29903)
  --skip-del, --no-skip-del
                        Skip deletions (default: True) (default: True)
  --skip-ins, --no-skip-ins
                        Skip insertions (default: True) (default: True)
```

## Logo

The breakfast logo is CC licensed:

Breakfast by Tippawan Sookruay from the Noun Project
