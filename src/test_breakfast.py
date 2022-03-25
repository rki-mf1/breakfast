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

# caching tests

def test_caching_export():
    os.system('python src/breakfast.py --input-file test/testfile.tsv --outdir test/ --max-dist 1 --output-cache test/testfile.json')
    assert exit_status == 0
    assert os.path.exists('test/testfile.json')

def test_caching_import():
    os.system('python src/breakfast.py --input-file test/testfile_caching01_AddedSeqs.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    assert exit_status == 0

def test_caching_01_AddingSequences():
    os.system('python src/breakfast.py --input-file test/testfile_caching01_AddedSeqs.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    expected_file = 'test/expected_clusters_caching01_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

def test_caching_02_DisorderedSequences():
    os.system('python src/breakfast.py --input-file test/testfile_caching02_DisorderedOnly.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    expected_file = 'test/expected_clusters_caching02_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

def test_caching_03_AddedAndDisorderedSequences():
    os.system('python src/breakfast.py --input-file test/testfile_caching03_AddedAndDisorderedSeqs.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    expected_file = 'test/expected_clusters_caching03_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

def test_caching_04_DeletingSingleSequences():
    os.system('python src/breakfast.py --input-file test/testfile_caching04_DeletedSingleSeq.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    expected_file = 'test/expected_clusters_caching04_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

# Complete Profile gets deleted
def test_caching_05_DeletedMultipleSequences():
    os.system('python src/breakfast.py --input-file test/testfile_caching05_DeletedProfile.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    expected_file = 'test/expected_clusters_caching05_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

def test_caching_06_ModifiedSequences():
    os.system('python src/breakfast.py --input-file test/testfile_caching06_ModifiedSeqs.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    expected_file = 'test/expected_clusters_caching06_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

def test_caching_07_MultipleTests():
    os.system('python src/breakfast.py --input-file test/testfile_caching07_MultipleTests.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    expected_file = 'test/expected_clusters_caching07_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()

def test_caching_08_MultipleTests():
    os.system('python src/breakfast.py --input-file test/testfile_caching08_NoChanges.tsv --outdir test/ --max-dist 1 --input-cache test/testfile.json')
    expected_file = 'test/expected_clusters_caching08_dist1.tsv'
    output_file = 'test/clusters.tsv'
    e = open(expected_file, "r")
    o = open(output_file, "r")
    assert e.read() == o.read()  