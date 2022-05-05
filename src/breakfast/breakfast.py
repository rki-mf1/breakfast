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


def remove_indels(features, feature_sep, feature_type, skip_ins, skip_del, trim_start, trim_end, reference_length):
    filtered_features = []
    insertion = re.compile(".*[A-Z][A-Z]$")
    for feature in features:
        if isinstance(feature, float):
            d = []
        else:
            if feature.find(feature_sep) != -1:
                d = feature.split(feature_sep)
            else:
                d = [feature]
        new_d = []
        for term in d:
            if feature_type == "dna":
                if term.startswith("del"):
                    if skip_del:
                        continue
                    pos = int(term.split(":")[1])
                    if ((trim_start is not None) and (pos <= trim_start)) or (
                        (trim_end is not None) and (pos >= (reference_length - trim_end))
                    ):
                        continue
                elif skip_ins and insertion.match(term) is not None:
                    continue
                else:
                    # Blindly remove reference and alt NT, leaving the position. Then
                    # check if it is in the regions we want to trim away
                    pos = int(term.translate(str.maketrans("", "", "ACGTN")))
                    if ((trim_start is not None) and (pos <= trim_start)) or (
                        (trim_end is not None)
                        and (pos >= (reference_length - trim_end))
                    ):
                        continue
            new_d.append(term)
        filtered_features.append(" ".join(new_d))
    return filtered_features


def construct_sub_mat(features, feature_sep):
    print("Convert list of substitutions into a sparse matrix")
    insertion = re.compile(".*[A-Z][A-Z]$")
    indptr = [0]
    indices = []
    data = []
    vocabulary = {}
    for subt in features:
        if isinstance(subt, float):
            d = []
        else:
            if subt.find(feature_sep) != -1:
                d = subt.split(feature_sep)
            else:
                d = [subt]
        for term in d:
            index = vocabulary.setdefault(term, len(vocabulary))
            indices.append(index)
            data.append(1)
        indptr.append(len(indices))
    sub_mat = csr_matrix((data, indices, indptr), dtype=int)
    return sub_mat


def calc_sparse_matrix(meta,
                       feature_sep,
                       max_dist,
                       min_cluster_size,
                       input_cache,
                       output_cache):
    # IMPORT RESULTS FROM PREVIOUS RUN
    try:
        with gzip.open(input_cache, "rb") as f:
            print("Import from pickle file")
            cache = cPickle.load(f)

        ca.validate(cache, max_dist, __version__)

        feature_map = ca.map_features(cache["meta"]["feature"], meta["feature"])
        neigh_cache_updated = ca.update_neighbours(cache["neigh"], feature_map)

        # construct sub_mat of complete dataset and sub_mat of only new sequences compared to cached meta
        idx_only_new = ca.find_new(feature_map)
        select_ind = np.array(idx_only_new)
        sub_mat = construct_sub_mat(meta["feature"], feature_sep)
        sub_mat_only_new_seqs = sub_mat[select_ind, :]

        print("Use sparse matrix to calculate pairwise distances, bounded by max_dist")

        def _reduce_func(D_chunk, start):
            neigh = [np.flatnonzero(d <= max_dist) for d in D_chunk]
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

        sub_mat = construct_sub_mat(meta["feature"], feature_sep)

        print("Use sparse matrix to calculate pairwise distances, bounded by max_dist")

        def _reduce_func(D_chunk, start):
            neigh = [np.flatnonzero(d <= max_dist) for d in D_chunk]
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
            "max_dist": max_dist,
            "version": __version__,
            "neigh": neigh,
            "meta": meta[["id", "feature"]],
        }
        with gzip.open(output_cache, "wb") as f:
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
        if clust_len >= min_cluster_size:
            cluster_id += 1
            meta.iloc[list(clust), meta.columns.get_loc("cluster_id")] = cluster_id
    print(f"Number of clusters found: {cluster_id}")
    return meta


# TODO: Caching results for max-dist 0
def calc_without_sparse_matrix(meta, min_cluster_size):
    print("Skip sparse matrix calculation since max-dist = 0")
    clusters = list(range(0, len(meta)))
    accession_list = meta["id"].tolist()
    meta["cluster_id"] = pd.NA
    cluster_id = 0
    for clust in clusters:
        clust_len = len(accession_list[clust])
        if clust_len >= min_cluster_size:
            cluster_id += 1
            meta.iloc[clust, meta.columns.get_loc("cluster_id")] = cluster_id
    print(f"Number of clusters found: {cluster_id}")
    return meta
