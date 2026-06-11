from breakfast import breakfast


def test_filter_correct():
    features = ["C241T"]
    expected = ["C241T"]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_dna", False, False, 0, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_trim():
    features = ["C241T"]
    expected = [""]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_dna", False, False, 250, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_empty():
    features = [""]
    expected = [""]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_dna", False, False, 250, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_whitespace():
    features = ["  "]
    expected = [""]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_dna", False, False, 250, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_skipdel():
    features = ["C241T del:10:1 G5343TT"]
    expected = ["C241T G5343TT"]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_dna", False, True, 0, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_skipins():
    features = ["C241T del:10:1 G5343TT"]
    expected = ["C241T del:10:1"]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_dna", True, False, 0, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_skipindel():
    features = ["C241T del:10:1 G5343TT"]
    expected = ["C241T"]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_dna", True, True, 0, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_allfilters():
    features = ["G24C C241T del:10:1 G533TT A990T"]
    expected = ["C241T"]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_dna", True, True, 100, 100, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_covsonar_aa_skipindel():
    features = ["S:N501Y ORF1:del:12:7 N:A34AK"]
    expected = ["S:N501Y"]
    filtered = breakfast.filter_features(
        features, " ", "covsonar_aa", True, True, 0, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_filter_nextclade_aa_skipdel():
    features = ["S:N501Y S:V70-"]
    expected = ["S:N501Y"]
    filtered = breakfast.filter_features(
        features, " ", "nextclade_aa", False, True, 0, 0, 1000
    )
    assert filtered[0] == expected[0]


def test_sparse_feature_matrix_ignores_empty_features():
    matrix = breakfast.sparse_feature_matrix(["", "C241T"], " ")

    assert matrix.shape == (2, 1)
    assert matrix[0].count_nonzero() == 0
    assert matrix[1].count_nonzero() == 1
