import pytest
import os

def test_entrypoint():
    exit_status = os.system('python src/breakfast.py --help')
    assert exit_status == 0

def test_simplerun():
    exit_status = os.system('python src/breakfast.py --input-file test/testfile.tsv --outdir test/')
    assert exit_status == 0

def test_file_generated():
    os.system('python src/breakfast.py --input-file test/testfile.tsv --outdir test/')
    assert os.path.exists('test/clusters.tsv')

def test_clustering_dist0():
    os.system('python src/breakfast.py --input-file test/testfile.tsv --max-dist 0 --outdir test/')
    expected_file = 'test/expected_clusters_dist0.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

def test_clustering_dist1():
    os.system('python src/breakfast.py --input-file test/testfile.tsv --max-dist 1 --outdir test/')
    expected_file = 'test/expected_clusters_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

def test_clustering_dist1_noskipdel():
    os.system('python src/breakfast.py --input-file test/testfile.tsv --max-dist 1 --outdir test/ --no-skip-del')
    expected_file = 'test/expected_clusters_dist1_noskipdel.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()    