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

import networkx
from networkx.algorithms.components.connected import connected_components


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--nextclade", help="Nextclade metadata (tsv format)")
    parser.add_argument("--clust-col", help="Nextclade metadata column to cluster (default = 'substitutions')")
    parser.add_argument("--outdir", help="Output directory for all output files")
    parser.add_argument("--max-dist", type=int, help="Maximum parwise distance")
    parser.add_argument("--min-cluster-size", type=int, help="Minimum cluster size")
    parser.add_argument("--jobs", type=int, help="Number of jobs to run in parallel")

    parser.set_defaults(
        nextclade="/scratch_slow/sc2/global/2021-05-02/clades/clades-nextstrain.tsv",
        clust_col="substitutions",
        outdir="../output/latest",
        max_dist=1,
        min_cluster_size=4,
        jobs=20,
    )

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("Clustering sequences")
    print(f"  max dist = {args.max_dist}")
    print(f"  minimum cluster size = {args.min_cluster_size}")
    print(f"  nextclade metadata file = {args.nextclade}")
    print(f"  nextclade metadata column = {args.clust_col}")

    meta = pd.read_table(args.nextclade, dtype={args.clust_col: str})
    print(f"Number of sequences: {meta.shape[0]}")

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
            if subt.find(",") != -1:
                d = subt.split(",")
            else:
                d = [subt]

        for term in d:
            index = vocabulary.setdefault(term, len(vocabulary))
            indices.append(index)
            data.append(1)
        indptr.append(len(indices))

    sub_mat = csr_matrix((data, indices, indptr), dtype=int)
    num_nz = sub_mat.getnnz()
    print(f"Number of non-zero matrix entries: {num_nz} ({100 * num_nz/(sub_mat.shape[0] * sub_mat.shape[1])}%)")

    print("Use sparse matrix to calculate pairwise distances, bounded by max_dist")
    def reduce_func(D_chunk, start):
        neigh = [np.flatnonzero(d <= args.max_dist) for d in D_chunk]
        return neigh

    gen = pairwise_distances_chunked(
        sub_mat, reduce_func=reduce_func, metric="manhattan"
    )
        #sub_mat, reduce_func=reduce_func, metric="manhattan", n_jobs=args.jobs, working_memory=2048

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
    meta["cluster_id"] = np.nan
    cluster_id = 0
    for clust in clusters:
        if len(clust) >= args.min_cluster_size:
            cluster_id += 1
            meta.iloc[list(clust), meta.columns.get_loc("cluster_id")] = cluster_id

    meta.to_csv(
        os.path.join(args.outdir, "meta-with-clusters.tsv"), sep="\t", index=False
    )

    meta[["seqName", "clade", "totalMutations", "substitutions", "aaSubstitutions", "cluster_id"]].to_csv(
        os.path.join(args.outdir, "meta-with-clusters-simpler.tsv"), sep="\t", index=False
    )

    meta[["seqName", "cluster_id"]].to_csv(
        os.path.join(args.outdir, "clusters.tsv"), sep="\t", index=False
    )

    print(f"Number of clusters found: {cluster_id}")


# Main body
if __name__ == "__main__":
    main()
