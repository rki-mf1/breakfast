#!/usr/bin/env python

import argparse
import collections
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


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input-file", help="Input file")
    parser.add_argument("--id-col", help="Column with the sequence identifier")
    parser.add_argument("--clust-col", help="Metadata column to cluster (default = 'dna_profile')")
    parser.add_argument("--sep", help="Clustering column separator")
    parser.add_argument("--outdir", default="output", help="Output directory for all output files")
    parser.add_argument("--max-dist", type=int, default=1, help="Maximum parwise distance")
    parser.add_argument("--min-cluster-size", type=int, default=2, help="Minimum cluster size")
    parser.add_argument("--trim-start", type=int, default=100, help="Bases to trim from the beginning")
    parser.add_argument("--trim-end", type=int, default=50, help="Bases to trim from the end")
    parser.add_argument("--reference-length", type=int, default=29903, help="Length of reference genome (defaults to NC_045512.2 length = 29903)")

    parser.set_defaults(
        input_file="../input/covsonar/rki-2021-05-19-minimal.tsv.gz",
        id_col="accession",
        clust_col="dna_profile",
        sep=" ",
    )

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("Clustering sequences")
    print(f"  Input file = {args.input_file}")
    print(f"  ID column = {args.id_col}")
    print(f"  clustering feature column = {args.clust_col}")
    print(f"  clustering feature column separator = {args.sep}")
    print(f"  max dist = {args.max_dist}")
    print(f"  minimum cluster size = {args.min_cluster_size}")

    meta = pd.read_table(args.input_file, usecols=[args.id_col, args.clust_col], dtype={args.id_col: str, args.clust_col: str})
    print(f"Number of sequences: {meta.shape[0]}")

    insertion = re.compile(".*[A-Z][A-Z]$")
    print("Convert list of substitutions into a sparse matrix")
    subs = meta[args.clust_col]
    indptr = [0]
    indices = []
    data = []
    vocabulary = {}
    for subt in subs:
        if isinstance(subt, float):
            d = []
        else:
            if subt.find(args.sep) != -1:
                d = subt.split(args.sep)
            else:
                d = [subt]

        for term in d:
            # Skip all deletions by default
            if term.startswith("del"):
                #print(f"skipping: {term}");
                continue

            if insertion.match(term) is not None:
                print(f"skipping: {term}");
                continue

            # Blindly remove reference and alt NT, leaving the position. Then
            # check if it is in the regions we want to trim away
            pos = int(term.translate(str.maketrans('', '', 'ACGTN')))
            #pos = int(term.replace("ACGTN", ""))
            if ((args.trim_start is not None) and (pos <= args.trim_start)) or ((args.trim_end is not None) and (pos >= (args.reference_length - args.trim_end))):
                #print(f"skipping: {term}");
                continue

            index = vocabulary.setdefault(term, len(vocabulary))
            indices.append(index)
            data.append(1)
        indptr.append(len(indices))

    sub_mat = csr_matrix((data, indices, indptr), dtype=int)
    num_nz = sub_mat.getnnz()
    print(f"Number of non-zero matrix entries: {num_nz} ({100 * num_nz/(sub_mat.shape[0] * sub_mat.shape[1]):.2f}%)")

    print("Use sparse matrix to calculate pairwise distances, bounded by max_dist")
    def reduce_func(D_chunk, start):
        neigh = [np.flatnonzero(d <= args.max_dist) for d in D_chunk]
        return neigh

    gen = pairwise_distances_chunked(
        sub_mat, reduce_func=reduce_func, metric="manhattan"
    )

    neigh = list(chain.from_iterable(gen))

    # Merge connected components
    def to_graph(l):
        G = networkx.Graph()
        for part in l:
            # each sublist is a bunch of nodes
            G.add_nodes_from(part)
            # it also imlies a number of edges:
            G.add_edges_from(to_edges(part))
        return G

    def to_edges(l):
        """
        treat `l` as a Graph and returns it's edges
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
        """
        it = iter(l)
        last = next(it)

        for current in it:
            yield last, current
            last = current

    print("Create graph and recover connected components")
    G = to_graph(neigh)
    clusters = connected_components(G)

    print("Save clusters")
    meta["cluster_id"] = pd.NA
    cluster_id = 0
    for clust in clusters:
        if len(clust) >= args.min_cluster_size:
            cluster_id += 1
            meta.iloc[list(clust), meta.columns.get_loc("cluster_id")] = cluster_id

    meta[[args.id_col, "cluster_id"]].to_csv(
        os.path.join(args.outdir, "clusters.tsv"), sep="\t", index=False
    )

    print(f"Number of clusters found: {cluster_id}")


# Main body
if __name__ == "__main__":
    main()
