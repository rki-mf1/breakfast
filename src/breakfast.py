#!/usr/bin/env python

VERSION = "0.3.0"

import argparse
import gzip
import os
import re
import sys
from itertools import chain

import _pickle as cPickle
import networkx
import numpy as np
import pandas as pd
from networkx.algorithms.components.connected import connected_components
from scipy.sparse import csr_matrix
from sklearn.metrics import pairwise_distances_chunked

import cache as ca


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
    return sub_mat


def calc_sparse_matrix(meta, args):
    # IMPORT RESULTS FROM PREVIOUS RUN
    try:
        with gzip.open(args.input_cache, "rb") as f:
            print("Import from pickle file")
            cache = cPickle.load(f)

        ca.validate(cache, args, VERSION)

        del_seqs, new_seqs, common_seqs = ca.find_deleted(cache, meta)
        mod_seqs = ca.find_modified(cache, meta)

        # Update the deleted seqs list with the sequences that were
        # modified, and also add them to the "new_seqs" list so that they
        # are properly consdiered "new"
        del_seqs += mod_seqs
        new_seqs += mod_seqs

        # What now can happen is the following
        # cached IDs: [(ID1, ID2), (ID3)]
        # new IDs: [(ID4), (ID5, ID6)]
        # but f.e ID4 will get grouped with ID1, ID2 because they all have identical sequences
        # current IDs: [(ID1, ID2, ID4), (ID3), (ID5, ID6)]
        # so the following loop goes trough every set to check if the sequences are completly new
        # or already known and gives back only new IDs
        # in the example above that would be only (ID5, ID6) because the cluster_id of ID4 is already known
        # same as (ID2, ID1)
        new_seqs_grouped = []
        new_seqs_grouped_idx = []
        del_seqs_grouped_idx = []
        current_id = meta["id"].tolist()
        # go trough all current ID sets
        for grouped_seqs_idx, grouped_seqs in enumerate(current_id):
            counter_new_seq = 0
            # are all set members new?
            for seq in grouped_seqs:
                if seq in new_seqs:
                    counter_new_seq += 1
            # only append to list if all elements are new
            if len(grouped_seqs) == counter_new_seq:
                new_seqs_grouped.append(grouped_seqs)
                new_seqs_grouped_idx.append(grouped_seqs_idx)

        meta_cached = cache["meta"]
        cached_id = meta_cached["id"]
        if (len(del_seqs)) > 0:
            print(f"{len(del_seqs)} deleted sequence(s)")
            print(f"The following sequences got deleted {del_seqs}")
            # Check if a unqiue sequence got deleted
            # this would change the index for neigh
            # Example 1
            # cached IDs: [(ID1, ID2,), (ID3), ...]
            # deleted sequence: ID2
            # This would not alter the cluster indices because we still have ID1 at the same position
            # Example 2
            # same as above but deleted sequence: (ID3)
            # This would affect the indices and results which come after ID3
            for grouped_seqs_idx, grouped_seqs in enumerate(cached_id):
                counter_del_seq = 0
                for seq in grouped_seqs:
                    if seq in del_seqs:
                        counter_del_seq += 1
                if len(grouped_seqs) == counter_del_seq:
                    del_seqs_grouped_idx.append(grouped_seqs_idx)
            # remove deleted sequence from cached meta
            meta_cached_wihoutdeleted_seqs = meta_cached.drop(del_seqs_grouped_idx)
            # build dict to change indices for neigh list
            list_of_old_indices = meta_cached_wihoutdeleted_seqs.index.tolist()
            list_of_new_indices = []
            for idx, i in enumerate(list_of_old_indices):
                list_of_new_indices.append(idx)
            zip_iterator = zip(list_of_old_indices, list_of_new_indices)
            indices_dict = dict(zip_iterator)

        meta_only_new_seqs = meta.iloc[new_seqs_grouped_idx]

        print(
            f"Number of new unique sequences compared to cached results: {len(meta_only_new_seqs)}"
        )

        # construct sub_mat of complete dataset and sub_mat of only new sequences compared to cached meta
        sub_mat = construct_sub_mat(meta, args)
        select_ind = np.array(new_seqs_grouped_idx)
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
        # TODO: If sequences got deleted we would need to change neigh_new and neigh_cached
        # Lets say we have as similar example as above
        # cached IDs: [(ID1, ID2), (ID3), (ID3_dup)]
        # new IDs: [(ID4), (ID5, ID6), ... ]
        # deleted IDs: [(ID3_dup)]
        # current IDs: [(ID1, ID2, ID4), (ID3), (ID5, ID6), ...]
        # cached neigh: [(1), (2), (3)]
        # new neigh: [(3), ...]
        # because ID3 got deleted we need to adjust the index of cached neigh
        # cached neigh [(1), (2)]
        # new neigh [(3), ...]
        neigh_new = list(chain.from_iterable(gen))
        neigh_cached = cache["neigh"]

        if (len(del_seqs)) > 0:
            neigh_cached_updated = []
            # go trough each cached neight array
            for neigh_set in neigh_cached:
                neigh_set = list(neigh_set)
                # remove deleted indices
                neigh_set = [e for e in neigh_set if e not in set(del_seqs_grouped_idx)]
                # change remaining indices
                neigh_set_updatedidx = [indices_dict[ind] for ind in neigh_set]
                if neigh_set_updatedidx:
                    neigh_cached_updated.append(np.array(neigh_set_updatedidx))
            neigh = neigh_cached_updated + neigh_new
        else:
            neigh = neigh_cached + neigh_new

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
            "version": VERSION,
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

    parser.add_argument("--version", action="version", version="%(prog)s " + VERSION)

    parser.set_defaults(
        input_file="../input/covsonar/rki-2021-05-19-minimal.tsv.gz",
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
        meta = remove_indels(meta, args)

    print(f"Number of duplicates: {meta['feature'].duplicated().sum()}")

    # Group IDs with identical sequences together
    meta_withoutDUPS = meta.groupby("feature", as_index=False, sort=False).agg(
        {"id": lambda x: tuple(x), "feature": "first"}
    )
    print(f"Number of unique sequences: {meta_withoutDUPS.shape[0]}")

    if args.max_dist == 0:
        meta_withoutDUPS = calc_without_sparse_matrix(meta_withoutDUPS, args)
    else:
        meta_withoutDUPS = calc_sparse_matrix(meta_withoutDUPS, args)

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


# Main body
if __name__ == "__main__":
    main()
