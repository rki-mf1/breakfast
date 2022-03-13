import numpy as np
import pandas as pd

def validate(cache, args, version):
    max_dist_cached = cache["max_dist"]

    if args.max_dist != max_dist_cached:
        print(
            f"WARNING: Cached results were created using a differnt max-dist paramter"
        )
        print(
            f"Current max-dist parameter: {args.max_dist} \n Cached max-dist parameter: {max_dist_cached}"
        )
        raise UnboundLocalError()

    version_cached = cache["version"]
    if version_cached != version:
        print(
            f"WARNING: Cached results were created using breakfast version {version_cached}"
        )

def update_neighbours(neigh, fmap):
    """Update cached neighrours list to index into new data properly

    This includes deleting elements that have been deleted from the new data
    set, and also updating the indexes in the neigh list to properly point to
    the right data in the new data set.
    """
    c2n = fmap.set_index("idx_cache")
    neigh_updated = []
    for nlist in neigh:
        nlist_updated = c2n.loc[nlist, "idx_new"].dropna().astype(int).tolist()
        if len(nlist_updated) > 0:
            neigh_updated.append(nlist_updated)

    return neigh_updated


def find_deleted(feature_map):
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
    """Map between the indexes of the cached and new data, based on the features
    """

    cached_df = pd.DataFrame({
                            "idx": range(cached_feats.size),
                            "feature": cached_feats,
                            }).set_index("feature")
    new_df = pd.DataFrame({
                            "idx": range(new_feats.size),
                            "feature": new_feats,
                            }).set_index("feature")

    feat_map = cached_df.join(new_df, how="outer", lsuffix="_cache", rsuffix="_new")

    return feat_map

