from pathlib import Path

import click.testing
import pandas as pd
import pytest

from breakfast import console


@pytest.fixture
def runner():
    return click.testing.CliRunner()


@pytest.fixture
def init_cache(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    cache_path = tmp_path / "cache"
    result = runner.invoke(console.main, [
                           "--input-file", input_file,
                           "--outdir", str(tmp_path),
                           "--output-cache", str(cache_path),
                           "--max-dist", "1"])
    assert result.exit_code == 0
    assert cache_path.exists
    return cache_path


def do_cache_test(input_file, cache, tmp_path, expected_result, monkeypatch, runner):
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(console.main, [
                           "--input-file", input_file,
                           "--outdir", str(tmp_path),
                           "--input-cache", str(cache),
                           "--max-dist", "1"])
    assert result.exit_code == 0
    output = pd.read_table(tmp_path / "clusters.tsv", sep="\t")
    assert expected_result.equals(output)


def test_cache_added_seqs(monkeypatch, runner, tmp_path, init_cache):
    input_file = "testfile_caching01_AddedSeqs.tsv"
    expected = pd.read_table("expected_clusters_caching01_dist1.tsv", sep="\t")
    do_cache_test(input_file,
                  init_cache,
                  tmp_path,
                  expected,
                  monkeypatch,
                  runner)


def test_cache_disordered_only(monkeypatch, runner, tmp_path, init_cache):
    input_file = "testfile_caching02_DisorderedOnly.tsv"
    expected = pd.read_table("expected_clusters_caching02_dist1.tsv", sep="\t")
    do_cache_test(input_file,
                  init_cache,
                  tmp_path,
                  expected,
                  monkeypatch,
                  runner)


def test_cache_added_and_disordered_seqs(monkeypatch, runner, tmp_path, init_cache):
    input_file = "testfile_caching03_AddedAndDisorderedSeqs.tsv"
    expected = pd.read_table("expected_clusters_caching03_dist1.tsv", sep="\t")
    do_cache_test(input_file,
                  init_cache,
                  tmp_path,
                  expected,
                  monkeypatch,
                  runner)

def test_cache_deleted_single_seq(monkeypatch, runner, tmp_path, init_cache):
    input_file = "testfile_caching04_DeletedSingleSeq.tsv"
    expected = pd.read_table("expected_clusters_caching04_dist1.tsv", sep="\t")
    do_cache_test(input_file,
                  init_cache,
                  tmp_path,
                  expected,
                  monkeypatch,
                  runner)


def test_cache_deleted_profile(monkeypatch, runner, tmp_path, init_cache):
    input_file = "testfile_caching05_DeletedProfile.tsv"
    expected = pd.read_table("expected_clusters_caching05_dist1.tsv", sep="\t")
    do_cache_test(input_file,
                  init_cache,
                  tmp_path,
                  expected,
                  monkeypatch,
                  runner)


def test_cache_modified_seqs(monkeypatch, runner, tmp_path, init_cache):
    input_file = "testfile_caching06_ModifiedSeqs.tsv"
    expected = pd.read_table("expected_clusters_caching06_dist1.tsv", sep="\t")
    do_cache_test(input_file,
                  init_cache,
                  tmp_path,
                  expected,
                  monkeypatch,
                  runner)


def test_cache_multiple_tests(monkeypatch, runner, tmp_path, init_cache):
    input_file = "testfile_caching07_MultipleTests.tsv"
    expected = pd.read_table("expected_clusters_caching07_dist1.tsv", sep="\t")
    do_cache_test(input_file,
                  init_cache,
                  tmp_path,
                  expected,
                  monkeypatch,
                  runner)


def test_cache_no_changes(monkeypatch, runner, tmp_path, init_cache):
    input_file = "testfile_caching08_NoChanges.tsv"
    expected = pd.read_table("expected_clusters_caching08_dist1.tsv", sep="\t")
    do_cache_test(input_file,
                  init_cache,
                  tmp_path,
                  expected,
                  monkeypatch,
                  runner)
