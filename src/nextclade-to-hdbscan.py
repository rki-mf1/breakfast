#!/usr/bin/env python

import hdbscan
from scipy.sparse import csr_matrix
import pandas as pd
import math

meta = pd.read_table("../input/global-clades-nextstrain.tsv", dtype={'substitutions': str})
meta.info()

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

clusterer = hdbscan.HDBSCAN(min_cluster_size=5, min_samples=1, cluster_selection_epsilon=3)
clusts = clusterer.fit(sub_mat)

df_clust = pd.DataFrame({'ID':meta['seqName'], 'CLUSTER':clusterer.labels_})

df_clust.to_csv("clust.tsv", sep = "\t", index = False)
