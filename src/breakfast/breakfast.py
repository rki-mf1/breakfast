#!/usr/bin/env python

import gzip
import re
from itertools import chain

import _pickle as cPickle
import networkx
import numpy as np
import pandas as pd
from networkx.algorithms.components.connected import connected_components
from scipy.sparse import csr_matrix
from sklearn.metrics import pairwise_distances_chunked

from . import cache as ca, __version__


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
    subs = meta["feature"]
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
        new_sub.append(" ".join(new_d))
    meta["feature"] = new_sub
    return meta


def construct_sub_mat(meta, args):
    print("Convert list of substitutions into a sparse matrix")
    insertion = re.compile(".*[A-Z][A-Z]$")
    subs = meta["feature"]
    indptr = [0]
    indices = []
    data = []
    vocabulary = {}
    for subt in subs:
        if isinstance(subt, float):
            d = []
        else:
            if subt.find(" ") != -1:
                d = subt.split(" ")
            else:
                d = [subt]
        for term in d:
            if args.var_type == "dna":
                if len(term) == 0:
                    continue
                elif term.startswith("del"):
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
    return sub_mat


def calc_sparse_matrix(meta, args):
    # IMPORT RESULTS FROM PREVIOUS RUN
    try:
        with gzip.open(args.input_cache, "rb") as f:
            print("Import from pickle file")
            cache = cPickle.load(f)

        ca.validate(cache, args, __version__)

        feature_map = ca.map_features(cache["meta"]["feature"], meta["feature"])
        neigh_cache_updated = ca.update_neighbours(cache["neigh"], feature_map)

        # construct sub_mat of complete dataset and sub_mat of only new sequences compared to cached meta
        idx_only_new = ca.find_new(feature_map)
        select_ind = np.array(idx_only_new)
        sub_mat = construct_sub_mat(meta, args)
        sub_mat_only_new_seqs = sub_mat[select_ind, :]

        print("Use sparse matrix to calculate pairwise distances, bounded by max_dist")

        def _reduce_func(D_chunk, start):
            neigh = [np.flatnonzero(d <= args.max_dist) for d in D_chunk]
            return neigh

        gen = pairwise_distances_chunked(
            X=sub_mat_only_new_seqs,
            Y=sub_mat,
            reduce_func=_reduce_func,
            metric="manhattan",
            n_jobs=1,
        )

        neigh_new = list(chain.from_iterable(gen))
        neigh = neigh_cache_updated + neigh_new

    except (UnboundLocalError, TypeError) as e:
        print(
            "Imported cached results are not available. Distance matrix of complete dataset will be calculated."
        )

        sub_mat = construct_sub_mat(meta, args)

        print("Use sparse matrix to calculate pairwise distances, bounded by max_dist")

        def _reduce_func(D_chunk, start):
            neigh = [np.flatnonzero(d <= args.max_dist) for d in D_chunk]
            return neigh

        gen = pairwise_distances_chunked(
            sub_mat,
            reduce_func=_reduce_func,
            metric="manhattan",
            n_jobs=1,
        )
        neigh = list(chain.from_iterable(gen))

    # EXPORT RESULTS FOR CACHING
    try:
        print("Export results as pickle")
        d = {
            "max_dist": args.max_dist,
            "version": __version__,
            "neigh": neigh,
            "meta": meta[["id", "feature"]],
        }
        with gzip.open(args.output_cache, "wb") as f:
            cPickle.dump(d, f, 2)  # protocol 2, python > 2.3
    except TypeError:
        print("Export of pickle was not succesfull")

    print("Create graph and recover connected components")
    G = _to_graph(neigh)
    clusters = connected_components(G)

    print("Save clusters")
    meta["cluster_id"] = pd.NA
    cluster_id = 0
    accession_list = meta["id"].tolist()
    for clust in clusters:
        clust_len = 0
        for set_clust in clust:
            clust_len += len(accession_list[set_clust])
        if clust_len >= args.min_cluster_size:
            cluster_id += 1
            meta.iloc[list(clust), meta.columns.get_loc("cluster_id")] = cluster_id
    print(f"Number of clusters found: {cluster_id}")
    return meta


# TODO: Caching results for max-dist 0
def calc_without_sparse_matrix(meta, args):
    print("Skip sparse matrix calculation since max-dist = 0")
    clusters = list(range(0, len(meta)))
    accession_list = meta["id"].tolist()
    meta["cluster_id"] = pd.NA
    cluster_id = 0
    for clust in clusters:
        clust_len = len(accession_list[clust])
        if clust_len >= args.min_cluster_size:
            cluster_id += 1
            meta.iloc[clust, meta.columns.get_loc("cluster_id")] = cluster_id
    print(f"Number of clusters found: {cluster_id}")
    return meta
