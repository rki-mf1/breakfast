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
    result = runner.invoke(console.main,
                           ["--input-file", input_file, "--outdir", str(tmp_path)])
    assert result.exit_code == 0


def test_file_generated(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    runner.invoke(console.main, ["--input-file", input_file, "--outdir", str(tmp_path)])
    assert (tmp_path / "clusters.tsv").exists


def test_dist0(monkeypatch, runner, tmp_path):
    input_file = "testfile.tsv"
    monkeypatch.chdir(Path(__file__).parent)
    result = runner.invoke(console.main, [
                           "--input-file", input_file,
                           "--outdir", str(tmp_path),
                           "--max-dist", "0"])
    assert result.exit_code == 0
    print(result.output)
    expected_clustering = pd.read_table("expected_clusters_dist0.tsv", sep="\t")
    output_clustering = pd.read_table(tmp_path / "clusters.tsv", sep="\t")
    assert expected_clustering.equals(output_clustering)
