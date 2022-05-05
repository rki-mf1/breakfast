import os

import click
import pandas as pd

from . import breakfast, __version__


def aa_no_trim_skip(ctx, param, value):
    if value != "dna" and (ctx.params['trim_start'] != 0 or ctx.params['trim_end'] != 0):
        raise click.BadParameter("Can not trim non-DNA features")

    if value != "dna" and (ctx.params['skip_ins'] or ctx.params['skip_del']):
        raise click.BadParameter("Can not skip insertions and deletions in non-DNA features")


@click.command(context_settings={'show_default': True})
@click.option("--input-file", type=click.Path(exists=True), help="Input file", required=True)
@click.option("--sep", default="\t", help="Input file separator")
@click.option("--outdir", type=click.Path(), default="output", help="Output directory for all output files")
@click.option("--max-dist", type=int, default=1, help="Maximum parwise distance")
@click.option("--min-cluster-size", type=int, default=2, help="Minimum cluster size")
@click.option("--input-cache", type=click.Path(), help="Input cached pickle file from previous run")
@click.option("--output-cache", type=click.Path(), help="Path to Output cached pickle file")
@click.option("--id-col", default="accession", help="Column with the sequence identifier")
@click.option("--clust-col", default="dna_profile", help="Metadata column to cluster")
@click.option("--var-type", type=click.Choice(['dna', 'aa']), default="dna", help="Type of variants", callback=aa_no_trim_skip)
@click.option("--sep2", default=" ", help="Secondary clustering column separator (between each mutation)")
@click.option("--trim-start", type=int, default=264, help="Bases to trim from the beginning (0 = disable)")
@click.option("--trim-end", type=int, default=228, help="Bases to trim from the end (0 = disable)")
@click.option("--reference-length", type=int, default=29903, help="Length of reference genome (defaults to NC_045512.2 length)")
@click.option("--skip-del", is_flag=True, default=True, help="Skip deletions")
@click.option("--skip-ins", is_flag=True, default=True, help="Skip insertions")
@click.version_option(version=__version__)
def cli(input_file,
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
        ):

    print("Clustering sequences")
    print(f"  Input file = {input_file}")
    print(f"  Input file separator = '{sep}'")
    print(f"  ID column = {id_col}")
    print(f"  clustering feature type ('dna' or 'aa') = {var_type}")
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

    meta = pd.read_table(
        input_file,
        usecols=[id_col, clust_col],
        dtype={id_col: str, clust_col: str},
        sep=sep,
    ).rename(columns={id_col: "id", clust_col: "feature"})
    print(f"Number of sequences: {meta.shape[0]}")

    # Remove Indels from mutation profiles before grouping sequences together
    if skip_del or skip_ins or (trim_start > 0) or (trim_end > 0):
        meta["feature"] = breakfast.remove_indels(meta["feature"],
                                                  sep2,
                                                  var_type,
                                                  skip_ins,
                                                  skip_del,
                                                  trim_start,
                                                  trim_end,
                                                  reference_length)

    print(f"Number of duplicates: {meta['feature'].duplicated().sum()}")

    # Group IDs with identical sequences together
    meta_nodups = meta.groupby("feature", as_index=False, sort=False).agg(
        {"id": lambda x: tuple(x), "feature": "first"}
    )
    print(f"Number of unique sequences: {meta_nodups.shape[0]}")

    if max_dist == 0:
        meta_withoutDUPS = breakfast.calc_without_sparse_matrix(meta, min_cluster_size)
    else:
        meta_withoutDUPS = breakfast.calc_sparse_matrix(meta,
                                                        sep2,
                                                        max_dist,
                                                        min_cluster_size,
                                                        input_cache,
                                                        output_cache)

    # Assign correct ID
    meta_clusterid = []
    meta_accession = []
    accession_list = meta_withoutDUPS["id"].tolist()
    cluster_ids = meta_withoutDUPS["cluster_id"].tolist()
    for accession, clust_id in zip(accession_list, cluster_ids):
        for seq in accession:
            meta_accession.append(seq)
            meta_clusterid.append(clust_id)
    meta_out = pd.DataFrame()
    meta_out["id"] = meta_accession
    meta_out["cluster_id"] = meta_clusterid

    # Sort according to input file
    meta_out = meta_out.set_index("id")
    meta_out = meta_out.reindex(index=meta["id"])
    meta_out = meta_out.reset_index()

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    meta_out[["id", "cluster_id"]].to_csv(
        os.path.join(outdir, "clusters.tsv"), sep="\t", index=False
    )
