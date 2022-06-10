import gzip
import _pickle as cPickle

import pandas as pd

from . import __version__


def load(input_file, max_dist):
    """Load and validate cached data from a previous run"""
    with gzip.open(input_file, "rb") as f:
        print("Import from pickle file")
        cache = cPickle.load(f)
    validate(cache, max_dist, __version__)
    return cache


def save(output_file, neigh, meta, max_dist):
    """Save cached results to a file, along with some metadata"""
    try:
        print("Export results as pickle")
        d = {
            "max_dist": max_dist,
            "version": __version__,
            "neigh": neigh,
            "meta": meta[["id", "feature"]],
        }
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(output_file, "wb") as f:
            cPickle.dump(d, f, 2)  # protocol 2, python > 2.3
    except TypeError:
        print("Export of pickle was not succesfull")


def validate(cache, max_dist, version):
    max_dist_cached = cache["max_dist"]

    if max_dist != max_dist_cached:
        print("WARNING: Cached results were created using a differnt max-dist paramter")
        print(f"Current max-dist parameter: {max_dist}")
        print(f"Cached max-dist parameter: {max_dist_cached}")
        raise UnboundLocalError()

    version_cached = cache["version"]
    if version_cached != version:
        print(
            f"WARNING: Cached results were created using breakfast version {version_cached}"
        )


def update_neighbours(neigh, fmap):
    """Update cached neighbours list to index into new data properly

    This includes deleting elements that have been deleted from the new data
    set, and also updating the indexes in the neigh list to properly point to
    the right data in the new data set.
    """
    # Create a list that can do lookups in O(1)
    max_idx = int(max(fmap["idx_cache"].dropna()))
    c2n = [None] * (max_idx + 1)
    for idx, row in fmap.dropna().iterrows():
        c2n[int(row["idx_cache"])] = int(row["idx_new"])

    neigh_updated = []
    for i, nlist in enumerate(neigh):
        nlist_updated = [c2n[item] for item in nlist]
        nlist_updated = [x for x in nlist_updated if x is not None]
        if len(nlist_updated) > 0:
            neigh_updated.append(nlist_updated)

    return neigh_updated


def find_deleted(fmap):
    """Find elements in the cache that are not present in the new data

    Returns the index of the deleted elements in the cache, so that they can be
    removed from the neighbours (neigh) object
    """
    idx_cache_only = fmap[fmap["idx_new"].isna().tolist()]["idx_cache"]
    return idx_cache_only


def find_new(fmap):
    """Find elements in the cache that are not present in the new data

    Returns the index of the deleted elements in the cache, so that they can be
    removed from the neighbours (neigh) object
    """
    idx_new_only = fmap[fmap["idx_cache"].isna().tolist()]["idx_new"]
    return idx_new_only


def map_features(cached_feats, new_feats):
    """Map between the indexes of the cached and new data, based on the features"""

    cached_df = pd.DataFrame(
        {
            "idx": range(cached_feats.size),
            "feature": cached_feats,
        }
    ).set_index("feature")
    new_df = pd.DataFrame(
        {
            "idx": range(new_feats.size),
            "feature": new_feats,
        }
    ).set_index("feature")

    feat_map = cached_df.join(new_df, how="outer", lsuffix="_cache", rsuffix="_new")

    return feat_map
