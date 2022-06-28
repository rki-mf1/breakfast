#!/usr/bin/env python

import re
import sys
from itertools import chain

import networkx
import numpy as np
import pandas as pd
from networkx.algorithms.components.connected import connected_components
from scipy.sparse import csr_matrix

from . import cache as ca


def read_input(input_file, sep, id_col, feature_col):
    meta = pd.read_table(
        input_file,
        usecols=[id_col, feature_col],
        dtype={id_col: str, feature_col: str},
        sep=sep,
    ).rename(columns={id_col: "id", feature_col: "feature"})
    assert not any(meta.id.duplicated())
    meta.fillna(value={"feature": ""}, inplace=True)
    print(f"Number of sequences: {meta.shape[0]}")
    return meta


def write_output(meta_nodups, meta_original, outdir):
    # Assign correct ID
    meta_clusterid = []
    meta_accession = []
    accession_list = meta_nodups["id"].tolist()
    cluster_ids = meta_nodups["cluster_id"].tolist()
    for accession, clust_id in zip(accession_list, cluster_ids):
        for seq_id in accession:
            meta_accession.append(seq_id)
            meta_clusterid.append(clust_id)
    meta_out = pd.DataFrame()
    meta_out["id"] = meta_accession
    meta_out["cluster_id"] = meta_clusterid

    # Sort according to input file
    meta_out = meta_out.set_index("id")
    meta_out = meta_out.reindex(index=meta_original["id"])
    meta_out = meta_out.reset_index()

    # Assign new cluster IDs according to order of input file
    keys = list(dict.fromkeys(meta_out["cluster_id"].tolist()))
    keys = [x for x in keys if not pd.isna(x)]
    dict_id = {}
    new_cluster_id = 0
    for i in range(len(keys)):
        new_cluster_id += 1
        dict_id[keys[i]] = new_cluster_id
    replacer = dict_id.get
    meta_out["cluster_id"] = [replacer(n, n) for n in meta_out["cluster_id"].tolist()]

    assert meta_out.shape[0] == meta_original.shape[0]

    outdir.mkdir(parents=True, exist_ok=True)

    meta_out[["id", "cluster_id"]].to_csv(
        outdir / "clusters.tsv", sep="\t", index=False
    )


def collapse_duplicates(meta):
    print(f"Number of duplicates: {meta['feature'].duplicated().sum()}")
    # Group IDs with identical sequences together
    meta_nodups = meta.groupby("feature", as_index=False, sort=False).agg(
        {"id": lambda x: tuple(x), "feature": "first"}
    )
    print(f"Number of unique sequences: {meta_nodups.shape[0]}")
    return meta_nodups


def cluster(meta_nodups, sep2, max_dist, min_cluster_size, input_cache, output_cache):
    if max_dist == 0:
        meta_nodups = cluster_identical_features(meta_nodups, min_cluster_size)
    else:
        meta_nodups = cluster_features(
            meta_nodups, sep2, max_dist, min_cluster_size, input_cache, output_cache
        )
    return meta_nodups


# Merge connected components
def _to_graph(neighbour_list):
    G = networkx.Graph()
    for part in neighbour_list:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(_to_edges(part))
    return G


def _to_edges(clustered_ids):
    """Treat `clustered_ids` as a Graph and returns it's edges

    to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(clustered_ids)
    last = next(it)

    for current in it:
        yield last, current
        last = current


def filter_features(  # noqa: C901
    features,
    feature_sep,
    feature_type,
    skip_ins,
    skip_del,
    trim_start,
    trim_end,
    reference_length,
):
    # If there's no filtering to be done, we immediately return the original
    # feature list
    if not (skip_del or skip_ins or (trim_start > 0) or (trim_end > 0)):
        return features

    filtered_features = []

    # Set up the regexes we need based on which type of features the input is
    # using
    if feature_type == "covsonar_dna":
        substitution = re.compile(r"^[A-Z](\d+)[A-Z]$")
        insertion = re.compile(r"^.*[A-Z][A-Z]$")
        deletion = re.compile(r"^del:\d+:\d+$")
    elif feature_type == "covsonar_aa":
        substitution = re.compile(r"^[a-zA-Z0-9]+:[A-Z]\d+[A-Z]$")
        insertion = re.compile(r"^[a-zA-Z0-9]+:[A-Z]\d+[A-Z][A-Z]+$")
        deletion = re.compile(r"^[a-zA-Z0-9]+:del:\d+:\d+$")
    elif feature_type == "nextclade_dna":
        substitution = re.compile(r"^[A-Z](\d+)[A-Z]$")
        insertion = re.compile(r"^\d+:[A-Z]+$")
        deletion = re.compile(r"^\d+(-\d+)?$")
    elif feature_type == "nextclade_aa":
        substitution = re.compile(r"^[a-zA-Z0-9]+:[A-Z]\d+[A-Z*]$")
        insertion = re.compile(r"^$")
        deletion = re.compile(r"^[a-zA-Z0-9]+:[A-Z]\d+-$")
    elif feature_type == "raw":
        # Dummy regexes that will be skipped anyway
        substitution = re.compile(r"(?!x)x")
        insertion = re.compile(r"(?!x)x")
        deletion = re.compile(r"(?!x)x")
    else:
        print(
            f"The feature type (--var-type) you chose is not supported: '{feature_type}'"
        )
        sys.exit(1)
    for feature in features:
        d = feature.split(feature_sep)
        new_d = []
        for term in d:
            if feature_type != "raw":
                if submatch := substitution.match(term):
                    # Blindly remove reference and alt NT, leaving the position. Then
                    # check if it is in the regions we want to trim away
                    pos = int(submatch.group(1))
                    if (pos <= trim_start) or (pos >= (reference_length - trim_end)):
                        continue
                elif insertion.match(term):
                    if skip_ins:
                        continue
                elif deletion.match(term):
                    if skip_del:
                        continue
                else:
                    print(f"Skipping invalid feature: '{term}'")
                    continue
            # Skip features that are empty strings ("")
            if not term:
                continue
            new_d.append(term)
        filtered_features.append(feature_sep.join(new_d))
    return filtered_features


def sparse_feature_matrix(features, feature_sep):
    """Convert list of features (mutations) into a sparse matrix"""
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


def filter_where(arr, k_min, k_max):
    arr = arr[np.where(arr >= k_min) and np.where(arr <= k_max)]
    return arr


def get_neighbours_batch(
    fmatrix, n_features_all, n_features_query, max_dist, select_ind=None
):
    def _reduce_func(D_chunk, start):
        neigh = [np.flatnonzero(d <= max_dist) for d in D_chunk]
        return neigh

    # Because some values are possibly cached, we might want to select only
    # new/updated features for clustering. The next two if statements handle
    # this:

    # if select_ind is an empty list, so we return an empty neighbour list
    # immediately. All values are already in the cache or the input is empty.
    if select_ind is not None and select_ind.size == 0:
        return []

    # select_ind = None means we select everything (= no values are cached),
    # otherwise we only select a subset (= some values are cached)
    fmatrix_batch = fmatrix
    n_features_batch = n_features_all
    if select_ind is not None:
        fmatrix_batch = fmatrix[select_ind, :]
        n_features_batch = n_features_all[select_ind]

    # Optimization to not compute distances when the length of the feature
    # vectors indicate that they are guaranteed to be > max_dist from each
    # other.
    close_enough = np.isclose(n_features_all, n_features_query, atol=max_dist)
    close_enough_batch = np.isclose(n_features_batch, n_features_query, atol=max_dist)

    fmatrix = fmatrix[close_enough, :]
    fmatrix_batch = fmatrix_batch[close_enough_batch, :]

    # We "lazy import" this function because the number of threads
    # (=OMP_NUM_THREADS env variable) needs to be set before the import
    # happens, and we can't do that if the import is at the top of the file
    from sklearn.metrics import pairwise_distances_chunked

    gen = pairwise_distances_chunked(
        X=fmatrix_batch,
        Y=fmatrix,
        reduce_func=_reduce_func,
        metric="manhattan",
        n_jobs=1,
    )

    neigh = list(chain.from_iterable(gen))

    # Need to fix the index numbers in the neighbour array, because we have
    # probably subset the input data and we need the indexes to point to the
    # correct element of the original input array
    fix_idx = np.where(close_enough)[0]
    neigh = [fix_idx[x] for x in neigh]
    return neigh


def cluster_features(
    meta, feature_sep, max_dist, min_cluster_size, input_cache, output_cache
):

    # Create sparse matrix of all features once, which we will then subset
    # based on caching and optimization details
    feat_matrix = sparse_feature_matrix(meta["feature"], feature_sep)
    # Count the number of features using a row sum
    meta["n_features"] = feat_matrix.sum(axis=1)

    select_ind = None  # This actually means we select all seqs
    neigh_list = []
    n_features_unique = meta["n_features"].drop_duplicates().tolist()

    # Import results from previous run
    try:
        cache = ca.load(input_cache, max_dist)
        feature_map = ca.map_features(cache["meta"]["feature"], meta["feature"])
        neigh_cache_updated = ca.update_neighbours(cache["neigh"], feature_map)

        # Identify new/updated sequences and only select them for clustering
        select_ind = np.array(ca.find_new(feature_map)).astype(int)
        n_features_unique = meta["n_features"][select_ind].drop_duplicates().tolist()

        # Add the cached neighbours to our neighbour list
        neigh_list += neigh_cache_updated

    except (UnboundLocalError, TypeError):
        print(
            (
                "Imported cached results are not available. "
                "Distance matrix of complete dataset will be calculated."
            )
        )

    for n_features_query in n_features_unique:
        neigh = get_neighbours_batch(
            feat_matrix, meta["n_features"], n_features_query, max_dist, select_ind
        )
        neigh_list = neigh_list + neigh
    neigh = neigh_list

    if output_cache:
        ca.save(output_cache, neigh, meta, max_dist)

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


def cluster_identical_features(meta, min_cluster_size):
    """Identify identical feature sets quickly

    This is a special case for the situation where max_dist is 0, which means
    the feature sets need to be identical. In this case we can skip the
    creation of a sparse matrix and everything else, and just do a quick scan
    of the feature vectors. This is much faster and more memory efficient.

    Note that caching is not supported in this case.
    """
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
