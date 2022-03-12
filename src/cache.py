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

def find_deleted(cache, meta):
    # compare cached IDs and current IDs to check if sequences have been deleted
    meta_cached = cache["meta"]
    cached_seqs = meta_cached["feature"]
    cached_id = meta_cached["id"]
    current_id = meta["id"].tolist()

    # flatten list of IDs to compare cached and current IDs
    # since we summarized IDs with identical sequences to sets
    # e.g. [(ID1, ID2), (ID3)] -> {ID1, ID2, ID3}
    flat_cached_id = set(item for sublist in cached_id for item in sublist)
    flat_current_id = set(item for sublist in current_id for item in sublist)

    # get all IDs which are part of cached ID but not anymore of current ID
    del_seqs = list(flat_cached_id.difference(flat_current_id))
    new_seqs = list(flat_current_id.difference(flat_cached_id))
    common_seqs = list(flat_cached_id.intersection(flat_current_id))

    return del_seqs, new_seqs, common_seqs

def find_modified(cache, meta):
    # Undo the grouping/deduplication of sequences so we can work with
    # individual sequence ids
    def ungroup_df(df):
        df_ungrouped = pd.DataFrame(
            {
                "id": np.concatenate(df["id"].values),
                "feature": np.repeat(df["feature"].values, df["id"].str.len()),
            }
        )
        return df_ungrouped

    meta_cached = cache["meta"]
    meta_cached_ungrouped = ungroup_df(meta_cached).set_index("id")
    meta_new_ungrouped = ungroup_df(meta).set_index("id")

    meta_merged = meta_cached_ungrouped.join(
        meta_new_ungrouped, how="inner", lsuffix="_cached", rsuffix="_new"
    )
    meta_merged["modified"] = (
        meta_merged["feature_cached"] != meta_merged["feature_new"]
    )
    mod_seqs = meta_merged[meta_merged["modified"]].index.tolist()

    print(f"{len(mod_seqs)} modified sequence(s)")
    print(f"The following sequences were modified {mod_seqs}")

    return mod_seqs
