#!/usr/bin/env python
VERSION="0.3.0"

import argparse
import collections
import sys
from itertools import chain
from scipy.sparse import csr_matrix

from sklearn.metrics import pairwise_distances_chunked
from sklearn.metrics import pairwise_distances
import pandas as pd
import numpy as np
import os
import re


import networkx
from networkx.algorithms.components.connected import connected_components


from sklearn.metrics.pairwise import check_pairwise_arrays
from scipy.sparse import issparse
from sklearn.metrics._pairwise_fast import _sparse_manhattan
from scipy.spatial import distance
from sklearn.metrics.pairwise import *
from sklearn.metrics.pairwise import _check_chunk_size, _precompute_metric_params, _return_float_dtype
import sklearn.datasets
from sklearn.utils.validation import _num_samples
import Cython
from setuptools import setup
from Cython.Build import cythonize

import pyximport
pyximport.install()
import _pairwise_fast_breakfast
#from _pairwise_fast_breakfast import _sparse_manhattan_breakfast, _sparse_manhattan


def manhattan_distances_breakfast(X, max_dist, mutation_length_list, Y=None, *, sum_over_features=True):
    X, Y = check_pairwise_arrays(X, Y)
    if issparse(X) or issparse(Y):
        if not sum_over_features:
            raise TypeError(
                "sum_over_features=%r not supported for sparse matrices"
                % sum_over_features
            )

        X = csr_matrix(X, copy=False)
        Y = csr_matrix(Y, copy=False)
        X.sum_duplicates()  # this also sorts indices in-place
        Y.sum_duplicates()
        D = np.zeros((X.shape[0], Y.shape[0]))
        # TODO: Set args.max_dist and list of mutation lengths as parameters
        mutation_length_list = [float(i) for i in mutation_length_list]
        mutation_length_array = np.array(mutation_length_list)
        max_dist = float(max_dist)
        print(mutation_length_array)
        print(max_dist)
        #print(X.indptr)
        #print(len(X.indptr))
        #print(len(mutation_length_array))
        #print(X.shape[0], Y.shape[0])
        print(D.shape[0], D.shape[1])
        #print(X.shape[1], Y.shape[1])
        _sparse_manhattan(X.data, X.indices, X.indptr, Y.data, Y.indices, Y.indptr, D)
        #_sparse_manhattan_breakfast(X.data, X.indices, X.indptr, Y.data, Y.indices, Y.indptr, D, max_dist, mutation_length_array)
        return D

    if sum_over_features:
        return distance.cdist(X, Y, "cityblock")

    D = X[:, np.newaxis, :] - Y[np.newaxis, :, :]
    D = np.abs(D, D)
    return D.reshape((-1, X.shape[1]))

PAIRWISE_DISTANCE_FUNCTIONS = {
    # If updating this dictionary, update the doc in both distance_metrics()
    # and also in pairwise_distances()!
    "cityblock": manhattan_distances,
    "cosine": cosine_distances,
    "euclidean": euclidean_distances,
    "haversine": haversine_distances,
    "l2": euclidean_distances,
    "l1": manhattan_distances,
    "manhattan": manhattan_distances,
    "manhattan_breakfast": manhattan_distances_breakfast,
    "precomputed": None,  # HACK: precomputed is always allowed, never called
    "nan_euclidean": nan_euclidean_distances,
}

_VALID_METRICS = [
    "manhattan_breakfast",
    "euclidean",
    "l2",
    "l1",
    "manhattan",
    "cityblock",
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "cosine",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis",
    "matching",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule",
    "wminkowski",
    "nan_euclidean",
    "haversine",
]

def pairwise_distances_chunked(
    X,
    Y=None,
    *,
    reduce_func=None,
    metric="euclidean",
    n_jobs=None,
    working_memory=None,
    max_dist=None,
    mutation_length_list=[],
    **kwds,
):
    n_samples_X = _num_samples(X)
    if metric == "precomputed":
        slices = (slice(0, n_samples_X),)
    else:
        if Y is None:
            Y = X
        # We get as many rows as possible within our working_memory budget to
        # store len(Y) distances in each row of output.
        #
        # Note:
        #  - this will get at least 1 row, even if 1 row of distances will
        #    exceed working_memory.
        #  - this does not account for any temporary memory usage while
        #    calculating distances (e.g. difference of vectors in manhattan
        #    distance.
        chunk_n_rows = get_chunk_n_rows(
            row_bytes=8 * _num_samples(Y),
            max_n_rows=n_samples_X,
            working_memory=working_memory,
        )
        slices = gen_batches(n_samples_X, chunk_n_rows)

    # precompute data-derived metric params
    params = _precompute_metric_params(X, Y, metric=metric, **kwds)
    kwds.update(**params)

    #TODO: Slice mutation_length_list ....
    for sl in slices:
        print(sl)
        if sl.start == 0 and sl.stop == n_samples_X:
            X_chunk = X  # enable optimised paths for X is Y
            #TODO:
            mutation_length_list_chunk = mutation_length_list
        else:
            X_chunk = X[sl]
            #TODO: 
            mutation_length_list_chunk = mutation_length_list[sl]
        #print(X_chunk, mutation_length_list_chunk)
        D_chunk = pairwise_distances(X_chunk, max_dist, mutation_length_list, Y, metric=metric, n_jobs=n_jobs, **kwds)
        if (X is Y or Y is None) and PAIRWISE_DISTANCE_FUNCTIONS.get(
            metric, None
        ) is euclidean_distances:
            # zeroing diagonal, taking care of aliases of "euclidean",
            # i.e. "l2"
            D_chunk.flat[sl.start :: _num_samples(X) + 1] = 0
        if reduce_func is not None:
            chunk_size = D_chunk.shape[0]
            D_chunk = reduce_func(D_chunk, sl.start)
            _check_chunk_size(D_chunk, chunk_size)
        yield D_chunk

def pairwise_distances(
    X, max_dist, mutation_length_list, Y=None, metric="euclidean", *, n_jobs=None, force_all_finite=True, **kwds
):
    if (
        metric not in _VALID_METRICS
        and not callable(metric)
        and metric != "precomputed"
    ):
        raise ValueError(
            "Unknown metric %s. Valid metrics are %s, or 'precomputed', or a callable"
            % (metric, _VALID_METRICS)
        )

    if metric == "precomputed":
        X, _ = check_pairwise_arrays(
            X, Y, precomputed=True, force_all_finite=force_all_finite
        )

        whom = (
            "`pairwise_distances`. Precomputed distance "
            " need to have non-negative values."
        )
        check_non_negative(X, whom=whom)
        return X
    elif metric in PAIRWISE_DISTANCE_FUNCTIONS:
        func = PAIRWISE_DISTANCE_FUNCTIONS[metric]
    elif callable(metric):
        func = partial(
            _pairwise_callable, metric=metric, force_all_finite=force_all_finite, **kwds
        )
    else:
        if issparse(X) or issparse(Y):
            raise TypeError("scipy distance metrics do not support sparse matrices.")

        dtype = bool if metric in PAIRWISE_BOOLEAN_FUNCTIONS else None

        if dtype == bool and (X.dtype != bool or (Y is not None and Y.dtype != bool)):
            msg = "Data was converted to boolean for metric %s" % metric
            warnings.warn(msg, DataConversionWarning)

        X, Y = check_pairwise_arrays(
            X, Y, dtype=dtype, force_all_finite=force_all_finite
        )

        # precompute data-derived metric params
        params = _precompute_metric_params(X, Y, metric=metric, **kwds)
        kwds.update(**params)

        if effective_n_jobs(n_jobs) == 1 and X is Y:
            return distance.squareform(distance.pdist(X, metric=metric, **kwds))
        func = partial(distance.cdist, metric=metric, **kwds)

    return _parallel_pairwise(X, Y, func, n_jobs, max_dist,
    mutation_length_list, **kwds)

def _parallel_pairwise(X, Y, func, n_jobs, max_dist,
    mutation_length_list, **kwds):
    """Break the pairwise matrix in n_jobs even slices
    and compute them in parallel."""

    if Y is None:
        Y = X
    X, Y, dtype = _return_float_dtype(X, Y)

    if effective_n_jobs(n_jobs) == 1:
        #TODO: Add parameters here 
        return func(X, max_dist, mutation_length_list, Y, **kwds)
        #return func(X, Y, **kwds)

    # enforce a threading backend to prevent data communication overhead
    fd = delayed(_dist_wrapper)
    ret = np.empty((X.shape[0], Y.shape[0]), dtype=dtype, order="F")
    Parallel(backend="threading", n_jobs=n_jobs)(
        fd(func, ret, s, X, Y[s], **kwds)
        for s in gen_even_slices(_num_samples(Y), effective_n_jobs(n_jobs))
    )

    if (X is Y or Y is None) and func is euclidean_distances:
        # zeroing diagonal for euclidean norm.
        # TODO: do it also for other norms.
        np.fill_diagonal(ret, 0)

    return ret

########################################



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
    mutation_len = []
    vocabulary = {}
    for subt in subs:
        mutation_counter = 0
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
            mutation_counter += 1
            index = vocabulary.setdefault(term, len(vocabulary))
            indices.append(index)
            data.append(1)
        indptr.append(len(indices))
        mutation_len.append(mutation_counter)
    sub_mat = csr_matrix((data, indices, indptr), dtype=int)
    num_nz = sub_mat.getnnz()
    print(
        f"Number of non-zero matrix entries: {num_nz} ({100 * (num_nz/(sub_mat.shape[0] * sub_mat.shape[1])):.2f}%)"
    )

    print("Use sparse matrix to calculate pairwise distances, bounded by max_dist")


    # TODO: Save results from previous run? As JSON?

    # TODO: Skip distance matrix calculation if len(sequence1) - len (sequence2) > max.dist
    # need to change the sub_mat or even before that since D_chunk goes through all columns?
    # group into smaller matrices which share similar amount of sequences??

    # I would need to change start or remove sequences?
    def _reduce_func(D_chunk, start):
        neigh = [np.flatnonzero(d <= args.max_dist) for d in D_chunk]
        return neigh

    # TODO SOLUTIONS: CHANGE METRIC FUNCTION? 

    gen = pairwise_distances_chunked(
        sub_mat, 
        reduce_func=_reduce_func, 
        metric="manhattan_breakfast", 
        n_jobs=1,
        max_dist=args.max_dist,
        mutation_length_list=mutation_len
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
