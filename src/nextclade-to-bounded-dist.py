#!/usr/bin/env python

import argparse
import collections
import hdbscan
from itertools import chain
from scipy.sparse import csr_matrix

from sklearn.metrics import pairwise_distances_chunked
from sklearn.metrics import pairwise_distances
import sklearn.datasets
import pandas as pd
import numpy as np
import math
import os

import networkx
from networkx.algorithms.components.connected import connected_components


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--nextclade", help="Nextclade metadata (tsv format)")
    parser.add_argument("--outdir", help="Output directory for all output files")
    parser.add_argument("--max-dist", type=int, help="Maximum parwise distance")
    parser.add_argument("--min-cluster-size", type=int, help="Minimum cluster size")

    parser.set_defaults(
        nextclade="../input/global-clades-nextstrain.tsv",
        outdir="../output",
        max_dist=1,
        min_cluster_size=5,
    )

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    meta = pd.read_table(args.nextclade, dtype={"substitutions": str})

    print("Convert list of substitutions into a sparse matrix")
    subs = meta["substitutions"]
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
    meta["cluster_id"] = np.nan
    cluster_id = 1
    for clust in clusters:
        if len(clust) >= args.min_cluster_size:
            print("Found a cluster:" + str(cluster_id))
            print(clust)
            meta.iloc[list(clust), meta.columns.get_loc("cluster_id")] = cluster_id
            cluster_id += 1

    meta.to_csv(
        os.path.join(args.outdir, "meta-with-clusters.tsv"), sep="\t", index=False
    )

    print(meta.shape)


# Main body
if __name__ == "__main__":
    main()
