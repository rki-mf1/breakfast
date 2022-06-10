from pathlib import Path

import click.testing
import pandas as pd
import pytest

from breakfast import console


@pytest.fixture
def runner():
    return click.testing.CliRunner()


def test_entrypoint(runner):
    result = runner.invoke(console.main, ["--help"])
    assert result.exit_code == 0


def test_simplerun(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main, ["--input-file", input_file, "--outdir", str(tmp_path)]
    )
    assert result.exit_code == 0


def test_file_generated(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    runner.invoke(console.main, ["--input-file", input_file, "--outdir", str(tmp_path)])
    assert (tmp_path / "clusters.tsv").exists


def test_dist0(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        ["--input-file", input_file, "--outdir", str(tmp_path), "--max-dist", "0"],
    )
    assert result.exit_code == 0
    expected_clustering = pd.read_table("expected_clusters_dist0.tsv", sep="\t")
    output_clustering = pd.read_table(tmp_path / "clusters.tsv", sep="\t")
    assert expected_clustering.equals(output_clustering)


def test_dist1(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        ["--input-file", input_file, "--outdir", str(tmp_path), "--max-dist", "1"],
    )
    assert result.exit_code == 0
    expected_clustering = pd.read_table("expected_clusters_dist1.tsv", sep="\t")
    output_clustering = pd.read_table(tmp_path / "clusters.tsv", sep="\t")
    assert expected_clustering.equals(output_clustering)


def test_dist1_noskipdel(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        [
            "--input-file",
            input_file,
            "--outdir",
            str(tmp_path),
            "--max-dist",
            "1",
            "--no-skip-del",
        ],
    )
    assert result.exit_code == 0
    expected_clustering = pd.read_table(
        "expected_clusters_dist1_noskipdel.tsv", sep="\t"
    )
    output_clustering = pd.read_table(tmp_path / "clusters.tsv", sep="\t")
    assert expected_clustering.equals(output_clustering)


def test_duplicate_ids(monkeypatch, runner, tmp_path):
    input_file = "duplicate-ids.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        ["--input-file", input_file, "--outdir", str(tmp_path), "--max-dist", "1"],
    )
    assert result.exit_code != 0


def test_missing_feature_column(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        [
            "--input-file",
            input_file,
            "--outdir",
            str(tmp_path),
            "--max-dist",
            "1",
            "--cluster-col",
            "somethingmissing",
        ],
    )
    assert result.exit_code != 0


def test_missing_id_column(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        [
            "--input-file",
            input_file,
            "--outdir",
            str(tmp_path),
            "--max-dist",
            "1",
            "--id-col",
            "somethingmissing",
        ],
    )
    assert result.exit_code != 0


def test_raw_file_generated(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    runner.invoke(
        console.main,
        [
            "--input-file",
            input_file,
            "--outdir",
            str(tmp_path),
            "--trim-start",
            "0",
            "--trim-end",
            "0",
            "--no-skip-del",
            "--no-skip-ins" "--var-type",
            "raw",
        ],
    )
    assert (tmp_path / "clusters.tsv").exists


def test_raw_dist1(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        [
            "--input-file",
            input_file,
            "--outdir",
            str(tmp_path),
            "--trim-start",
            "0",
            "--trim-end",
            "0",
            "--no-skip-del",
            "--no-skip-ins",
            "--var-type",
            "raw",
        ],
    )
    assert result.exit_code == 0
    expected_clustering = pd.read_table(
        "expected_clusters_dist1_noskipdel.tsv", sep="\t"
    )
    output_clustering = pd.read_table(tmp_path / "clusters.tsv", sep="\t")
    assert expected_clustering.equals(output_clustering)


def test_nextclade_file_generated(monkeypatch, runner, tmp_path):
    input_file = "testfile_nextclade.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    runner.invoke(
        console.main,
        [
            "--input-file",
            input_file,
            "--outdir",
            str(tmp_path),
            "--sep2",
            ",",
            "--id-col",
            "seqName",
            "--clust-col",
            "substitutions",
            "--var-type",
            "nextclade_dna",
        ],
    )
    assert (tmp_path / "clusters.tsv").exists


def test_nextclade_dist0(monkeypatch, runner, tmp_path):
    input_file = "testfile_nextclade.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        [
            "--input-file",
            input_file,
            "--outdir",
            str(tmp_path),
            "--sep2",
            ",",
            "--id-col",
            "seqName",
            "--clust-col",
            "substitutions",
            "--max-dist",
            "0",
            "--var-type",
            "nextclade_dna",
        ],
    )
    assert result.exit_code == 0
    expected_clustering = pd.read_table("expected_clusters_dist0.tsv", sep="\t")
    output_clustering = pd.read_table(tmp_path / "clusters.tsv", sep="\t")
    assert expected_clustering.equals(output_clustering)


def test_nextclade_dist1(monkeypatch, runner, tmp_path):
    input_file = "testfile_nextclade.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(
        console.main,
        [
            "--input-file",
            input_file,
            "--outdir",
            str(tmp_path),
            "--sep2",
            ",",
            "--id-col",
            "seqName",
            "--clust-col",
            "substitutions",
            "--max-dist",
            "1",
            "--var-type",
            "nextclade_dna",
        ],
    )
    assert result.exit_code == 0
    expected_clustering = pd.read_table("expected_clusters_dist1.tsv", sep="\t")
    output_clustering = pd.read_table(tmp_path / "clusters.tsv", sep="\t")
    assert expected_clustering.equals(output_clustering)
