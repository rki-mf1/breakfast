#!/usr/bin/env python
VERSION="0.3.0"

import argparse
import collections
import sys
from itertools import chain
from scipy.sparse import csr_matrix

from sklearn.metrics import pairwise_distances_chunked
from sklearn.metrics import pairwise_distances
import sklearn.datasets
import pandas as pd
import numpy as np
import math
import os
import datetime
import re

import networkx
from networkx.algorithms.components.connected import connected_components


# Merge connected components
def _to_graph(l):
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(_to_edges(part))
    return G


def _to_edges(l):
    """
    treat `l` as a Graph and returns it's edges
    to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current


def remove_indels(meta, args):
        subs = meta[args.clust_col]
        new_sub = []
        insertion = re.compile(".*[A-Z][A-Z]$")
        for subt in subs:
            if isinstance(subt, float):
                d = []
            else:
                if subt.find(args.sep2) != -1:
                    d = subt.split(args.sep2)
                else:
                    d = [subt]
            new_d = []
            for term in d:
                if args.var_type == "dna":
                    if term.startswith("del"):
                        if args.skip_del:
                            continue
                        pos = int(term.split(":")[1])
                        if ((args.trim_start is not None) and (pos <= args.trim_start)) or (
                            (args.trim_end is not None)
                            and (pos >= (args.reference_length - args.trim_end))
                        ):
                            continue
                    elif args.skip_ins and insertion.match(term) is not None:
                        continue
                    else:
                        # Blindly remove reference and alt NT, leaving the position. Then
                        # check if it is in the regions we want to trim away
                        pos = int(term.translate(str.maketrans("", "", "ACGTN")))
                        if ((args.trim_start is not None) and (pos <= args.trim_start)) or (
                            (args.trim_end is not None)
                            and (pos >= (args.reference_length - args.trim_end))
                        ):
                            continue
                new_d.append(term)
            new_sub.append(' '.join(new_d))
        meta[args.clust_col] = new_sub
        return meta
    

def calc_sparse_matrix(meta_withoutDUPS, args):
    print("Convert list of substitutions into a sparse matrix")
    insertion = re.compile(".*[A-Z][A-Z]$")
    subs = meta_withoutDUPS[args.clust_col]
    indptr = [0]
    indices = []
    data = []
    vocabulary = {}
    for subt in subs:
        if isinstance(subt, float):
            d = []
        else:
            if subt.find(args.sep2) != -1:
                d = subt.split(args.sep2)
            else:
                d = [subt]
        for term in d:
            if args.var_type == "dna":
                if term.startswith("del"):
                    if args.skip_del:
                        continue
                    pos = int(term.split(":")[1])
                    if ((args.trim_start is not None) and (pos <= args.trim_start)) or (
                        (args.trim_end is not None)
                        and (pos >= (args.reference_length - args.trim_end))
                    ):
                        continue
                elif args.skip_ins and insertion.match(term) is not None:
                    continue
                else:
                    # Blindly remove reference and alt NT, leaving the position. Then
                    # check if it is in the regions we want to trim away
                    pos = int(term.translate(str.maketrans("", "", "ACGTN")))
                    if ((args.trim_start is not None) and (pos <= args.trim_start)) or (
                        (args.trim_end is not None)
                        and (pos >= (args.reference_length - args.trim_end))
                    ):
                        continue
            index = vocabulary.setdefault(term, len(vocabulary))
            indices.append(index)
            data.append(1)
        indptr.append(len(indices))
    sub_mat = csr_matrix((data, indices, indptr), dtype=int)
    num_nz = sub_mat.getnnz()
    print(
        f"Number of non-zero matrix entries: {num_nz} ({100 * (num_nz/(sub_mat.shape[0] * sub_mat.shape[1])):.2f}%)"
    )

    print("Use sparse matrix to calculate pairwise distances, bounded by max_dist")


    def _reduce_func(D_chunk, start):
        neigh = [np.flatnonzero(d <= args.max_dist) for d in D_chunk]
        return neigh

    gen = pairwise_distances_chunked(
        sub_mat, reduce_func=_reduce_func, metric="manhattan", n_jobs=1
    )

    neigh = list(chain.from_iterable(gen))


    print("Create graph and recover connected components")
    G = _to_graph(neigh)
    clusters = connected_components(G)


    print("Save clusters")
    meta_withoutDUPS['cluster_id'] = pd.NA
    cluster_id = 0
    accession_list = meta_withoutDUPS[args.id_col].tolist()
    for clust in clusters:
        clust_len = 0
        for set_clust in clust:
          clust_len += len(accession_list[set_clust])
        if clust_len >= args.min_cluster_size:
          cluster_id += 1
          meta_withoutDUPS.iloc[list(clust), meta_withoutDUPS.columns.get_loc("cluster_id")] = cluster_id
    print(f"Number of clusters found: {cluster_id}")
    return meta_withoutDUPS

def calc_without_sparse_matrix(meta_withoutDUPS, args):
    print("Skip sparse matrix calculation since max-dist = 0")
    clusters = list(range(0,len(meta_withoutDUPS)))
    accession_list = meta_withoutDUPS[args.id_col].tolist()
    meta_withoutDUPS['cluster_id'] = pd.NA
    cluster_id = 0
    for clust in clusters:
        clust_len = len(accession_list[clust])
        if clust_len >= args.min_cluster_size:
            cluster_id += 1
            meta_withoutDUPS.iloc[clust, meta_withoutDUPS.columns.get_loc("cluster_id")] = cluster_id
    print(f"Number of clusters found: {cluster_id}")
    return meta_withoutDUPS


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input-file", help="Input file")
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

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    parser.set_defaults(
        input_file="../input/covsonar/rki-2021-05-19-minimal.tsv.gz",
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
    
    if os.path.isfile(args.input_file):
      meta = pd.read_table(
          args.input_file,
          usecols=[args.id_col, args.clust_col],
          dtype={args.id_col: str, args.clust_col: str},
          sep=args.sep,
       )
    else:
      print(f"The input file {args.input_file} cannot be found!")
      exit()

    print(f"Number of sequences: {meta.shape[0]}")

    # Remove Indels from mutation profiles before grouping sequences together
    if args.skip_del or args.skip_ins:
        meta = remove_indels(meta, args)

    print(f"Number of duplicates: {meta[args.clust_col].duplicated().sum()}")

    # Group genomes with identical sequences together
    meta_withoutDUPS = meta.groupby(args.clust_col,as_index=False, sort=False).agg({args.id_col:lambda x : list(x),args.clust_col:'first'})

    print(f"Number of sequences without duplicates: {meta_withoutDUPS.shape[0]}")

    if args.max_dist == 0:
        meta_withoutDUPS = calc_without_sparse_matrix(meta_withoutDUPS, args)       
    else:
        meta_withoutDUPS = calc_sparse_matrix(meta_withoutDUPS, args)


    # Assign correct ID
    meta_clusterid = []
    meta_accession = []
    accession_list = meta_withoutDUPS[args.id_col].tolist()
    cluster_ids = meta_withoutDUPS['cluster_id'].tolist()
    for accession, clust_id in zip(accession_list, cluster_ids):
      for seq in accession:
        meta_accession.append(seq)
        meta_clusterid.append(clust_id)
    meta_out = pd.DataFrame()
    meta_out[args.id_col] = meta_accession
    meta_out['cluster_id'] = meta_clusterid

    # Sort according to input file
    meta_out = meta_out.set_index(args.id_col)
    meta_out = meta_out.reindex(index=meta[args.id_col])
    meta_out = meta_out.reset_index()

    meta_out[[args.id_col, "cluster_id"]].to_csv(
        os.path.join(args.outdir, "clusters.tsv"), sep="\t", index=False
    )

# Main body
if __name__ == "__main__":
    main()
