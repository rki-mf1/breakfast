import argparse
import os
import sys

import pandas as pd

from . import breakfast, __version__


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input-file", help="Input file", required=True)
    parser.add_argument(
        "--input-cache", help="Input cached pickle file from previous run"
    )
    parser.add_argument("--output-cache", help="Path to Output cached pickle file")
    parser.add_argument("--id-col", help="Column with the sequence identifier")
    parser.add_argument(
        "--clust-col", help="Metadata column to cluster (default = 'dna_profile')"
    )
    parser.add_argument(
        "--var-type",
        default="dna",
        help="Type of variants (dna or aa, default = 'dna')",
    )
    parser.add_argument(
        "--sep2",
        default=" ",
        help="Secondary clustering column separator (between each mutation)",
    )
    parser.add_argument("--sep", default="\t", help="Input file separator")
    parser.add_argument(
        "--outdir", default="output", help="Output directory for all output files"
    )
    parser.add_argument(
        "--max-dist", type=int, default=1, help="Maximum parwise distance"
    )
    parser.add_argument(
        "--min-cluster-size", type=int, default=2, help="Minimum cluster size"
    )
    parser.add_argument(
        "--trim-start",
        type=int,
        default=264,
        help="Bases to trim from the beginning (0 = disable)",
    )
    parser.add_argument(
        "--trim-end",
        type=int,
        default=228,
        help="Bases to trim from the end (0 = disable)",
    )
    parser.add_argument(
        "--reference-length",
        type=int,
        default=29903,
        help="Length of reference genome (defaults to NC_045512.2 length = 29903)",
    )
    parser.add_argument(
        "--skip-del",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Skip deletions",
    )
    parser.add_argument(
        "--skip-ins",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Skip insertions",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    parser.set_defaults(
        input_file="genomic_profiles.tsv.gz",
        input_cache=None,
        output_cache=None,
        id_col="accession",
        clust_col="dna_profile",
    )

    args = parser.parse_args()

    if args.var_type != "dna" and (args.trim_start != 0 or args.trim_end != 0):
        print(
            "Sequence trimming on non-dna variants is not currently implemented. Aborting."
        )
        sys.exit(1)

    if args.var_type != "dna" and (args.skip_ins or args.skip_del):
        print(
            "Skipping insertions and deletions for non-dna variants is not currently implemented. Aborting."
        )
        sys.exit(1)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("Clustering sequences")
    print(f"  Input file = {args.input_file}")
    print(f"  Input file separator = '{args.sep}'")
    print(f"  ID column = {args.id_col}")
    print(f"  clustering feature type ('dna' or 'aa') = {args.var_type}")
    print(f"  clustering feature column = {args.clust_col}")
    print(f"  clustering feature column separator = '{args.sep2}'")
    print(f"  max dist = {args.max_dist}")
    print(f"  minimum cluster size = {args.min_cluster_size}")
    print(f"  trim start (bp) = {args.trim_start}")
    print(f"  trim end (bp) = {args.trim_end}")
    print(f"  reference length (bp) = {args.reference_length}")
    print(f"  skip deletions = {args.skip_del}")
    print(f"  skip insertions = {args.skip_ins}")
    print(f"  Input cache file = {args.input_cache}")
    print(f"  Output cache file = {args.output_cache}")

    if os.path.isfile(args.input_file):
        meta = pd.read_table(
            args.input_file,
            usecols=[args.id_col, args.clust_col],
            dtype={args.id_col: str, args.clust_col: str},
            sep=args.sep,
        ).rename(columns={args.id_col: "id", args.clust_col: "feature"})
    else:
        print(f"The input file {args.input_file} cannot be found!")
        exit()

    print(f"Number of sequences: {meta.shape[0]}")

    # Remove Indels from mutation profiles before grouping sequences together
    if args.skip_del or args.skip_ins:
        meta = breakfast.remove_indels(meta, args)

    print(f"Number of duplicates: {meta['feature'].duplicated().sum()}")

    # Group IDs with identical sequences together
    meta_withoutDUPS = meta.groupby("feature", as_index=False, sort=False).agg(
        {"id": lambda x: tuple(x), "feature": "first"}
    )
    print(f"Number of unique sequences: {meta_withoutDUPS.shape[0]}")

    if args.max_dist == 0:
        meta_withoutDUPS = breakfast.calc_without_sparse_matrix(meta_withoutDUPS, args)
    else:
        meta_withoutDUPS = breakfast.calc_sparse_matrix(meta_withoutDUPS, args)

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

    meta_out[["id", "cluster_id"]].to_csv(
        os.path.join(args.outdir, "clusters.tsv"), sep="\t", index=False
    )
