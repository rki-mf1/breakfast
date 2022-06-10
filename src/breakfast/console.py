import os
import pathlib

import click

from . import breakfast, __version__


# fmt: off
@click.command(context_settings={'show_default': True})
@click.option("--input-file", type=click.Path(exists=True, path_type=pathlib.Path), help="Input file", required=True)
@click.option("--sep", default="\t", help="Input file separator")
@click.option("--outdir", type=click.Path(path_type=pathlib.Path), default="output", help="Output directory for all output files")
@click.option("--max-dist", type=click.IntRange(0), default=1, help="Maximum parwise distance")
@click.option("--min-cluster-size", type=click.IntRange(1), default=2, help="Minimum cluster size")
@click.option("--input-cache", type=click.Path(exists=True, path_type=pathlib.Path), help="Input cached pickle file from previous run")
@click.option("--output-cache", type=click.Path(path_type=pathlib.Path), help="Path to Output cached pickle file")
@click.option("--id-col", default="accession", help="Column with the sequence identifier")
@click.option("--clust-col", default="dna_profile", help="Metadata column to cluster")
@click.option("--var-type", type=click.Choice(['covsonar_dna', 'covsonar_aa', 'nextclade_dna', 'nextclade_aa', 'raw']), default="covsonar_dna", help="Type of variants")
@click.option("--sep2", default=" ", help="Secondary clustering column separator (between each mutation)")
@click.option("--trim-start", type=click.IntRange(0), default=264, help="Bases to trim from the beginning (0 = disable)")
@click.option("--trim-end", type=click.IntRange(0), default=228, help="Bases to trim from the end (0 = disable)")
@click.option("--reference-length", type=click.IntRange(0), default=29903, help="Length of reference genome (defaults to NC_045512.2 length)")
@click.option("--skip-del/--no-skip-del", default=True, help="Skip deletions")
@click.option("--skip-ins/--no-skip-ins", default=True, help="Skip insertions")
@click.option("--jobs", type=click.IntRange(1), default=1, envvar="OMP_NUM_THREADS", help="Number of jobs (threads) to run simultaneously")
@click.version_option(version=__version__)
# fmt: on
def main(
    input_file,
    outdir,
    input_cache,
    output_cache,
    id_col,
    clust_col,
    var_type,
    sep,
    sep2,
    max_dist,
    min_cluster_size,
    trim_start,
    trim_end,
    reference_length,
    skip_del,
    skip_ins,
    jobs,
):

    # Trimming and skipping indels is not currently implemented for amino acid features
    if var_type not in ["covsonar_dna", "nextclade_dna"]:
        if trim_start != 0 or trim_end != 0:
            raise click.BadParameter("Can not trim non-DNA features")
        if skip_ins or skip_del:
            raise click.BadParameter("Can not skip indels in non-DNA features")

    # Can't trim more than the reference length
    if trim_start > reference_length or trim_end > reference_length:
        raise click.BadParameter("Can not trim more than the reference length")

    print("Clustering sequences")
    print(f"  Input file = {input_file}")
    print(f"  Input file separator = '{sep}'")
    print(f"  ID column = {id_col}")
    print(f"  clustering feature type = {var_type}")
    print(f"  clustering feature column = {clust_col}")
    print(f"  clustering feature column separator = '{sep2}'")
    print(f"  max dist = {max_dist}")
    print(f"  minimum cluster size = {min_cluster_size}")
    print(f"  trim start (bp) = {trim_start}")
    print(f"  trim end (bp) = {trim_end}")
    print(f"  reference length (bp) = {reference_length}")
    print(f"  skip deletions = {skip_del}")
    print(f"  skip insertions = {skip_ins}")
    print(f"  Input cache file = {input_cache}")
    print(f"  Output cache file = {output_cache}")

    # Set the number of parallel threads to use for scikit learn distance
    # calculation
    os.environ["OMP_NUM_THREADS"] = str(jobs)

    meta = breakfast.read_input(input_file, sep, id_col, clust_col)

    meta["feature"] = breakfast.filter_features(
        meta["feature"],
        sep2,
        var_type,
        skip_ins,
        skip_del,
        trim_start,
        trim_end,
        reference_length,
    )
    meta_nodups = breakfast.collapse_duplicates(meta)

    meta_clustered = breakfast.cluster(
        meta_nodups, sep2, max_dist, min_cluster_size, input_cache, output_cache
    )

    breakfast.write_output(meta_clustered, meta, outdir)
