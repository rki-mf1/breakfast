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


# def test_filter_invalid_dna_sub():
#     features = ["C241t"]
#     expected = ["C241T"]
#     filtered = breakfast.filter_features(features, ",", "covsonar_dna", True, False, 0, 0, 1000)
#     assert filtered[0] == expected[0]
