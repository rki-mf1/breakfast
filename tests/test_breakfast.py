from pathlib import Path

import click.testing
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
