import os


# def test_entrypoint():
#     exit_status = os.system('breakfast --help')
#     assert exit_status == 0
# 
# 
# def test_simplerun():
#     exit_status = os.system('breakfast --input-file tests/testfile.tsv --outdir tests/')
#     assert exit_status == 0
# 
# 
# def test_file_generated():
#     os.system('breakfast --input-file tests/testfile.tsv --outdir tests/')
#     assert os.path.exists('tests/clusters.tsv')
# 
# 
# def test_clustering_dist0():
#     os.system('breakfast --input-file tests/testfile.tsv --max-dist 0 --outdir tests/')
#     expected_file = 'tests/expected_clusters_dist0.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# def test_clustering_dist1():
#     os.system('breakfast --input-file tests/testfile.tsv --max-dist 1 --outdir tests/')
#     expected_file = 'tests/expected_clusters_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# def test_clustering_dist1_noskipdel():
#     os.system('breakfast --input-file tests/testfile.tsv --max-dist 1 --outdir tests/ --no-skip-del')
#     expected_file = 'tests/expected_clusters_dist1_noskipdel.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# # caching tests
# def test_caching_export():
#     exit_status = os.system('breakfast --input-file tests/testfile.tsv --outdir tests/ --max-dist 1 --output-cache tests/testfile.json')
#     assert exit_status == 0
#     assert os.path.exists('tests/testfile.json')
# 
# 
# def test_caching_import():
#     exit_status = os.system('breakfast --input-file tests/testfile_caching01_AddedSeqs.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     assert exit_status == 0
# 
# 
# def test_caching_01_AddingSequences():
#     os.system('breakfast --input-file tests/testfile_caching01_AddedSeqs.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     expected_file = 'tests/expected_clusters_caching01_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# def test_caching_02_DisorderedSequences():
#     os.system('breakfast --input-file tests/testfile_caching02_DisorderedOnly.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     expected_file = 'tests/expected_clusters_caching02_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# def test_caching_03_AddedAndDisorderedSequences():
#     os.system('breakfast --input-file tests/testfile_caching03_AddedAndDisorderedSeqs.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     expected_file = 'tests/expected_clusters_caching03_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# def test_caching_04_DeletingSingleSequences():
#     os.system('breakfast --input-file tests/testfile_caching04_DeletedSingleSeq.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     expected_file = 'tests/expected_clusters_caching04_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# # Complete Profile gets deleted
# def test_caching_05_DeletedMultipleSequences():
#     os.system('breakfast --input-file tests/testfile_caching05_DeletedProfile.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     expected_file = 'tests/expected_clusters_caching05_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# def test_caching_06_ModifiedSequences():
#     os.system('breakfast --input-file tests/testfile_caching06_ModifiedSeqs.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     expected_file = 'tests/expected_clusters_caching06_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# def test_caching_07_MultipleTests():
#     os.system('breakfast --input-file tests/testfile_caching07_MultipleTests.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     expected_file = 'tests/expected_clusters_caching07_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
# 
# def test_caching_08_NoChanges():
#     os.system('breakfast --input-file tests/testfile_caching08_NoChanges.tsv --outdir tests/ --max-dist 1 --input-cache tests/testfile.json')
#     expected_file = 'tests/expected_clusters_caching08_dist1.tsv'
#     output_file = 'tests/clusters.tsv'
#     e = open(expected_file, "r")
#     o = open(output_file, "r")
#     assert e.read() == o.read()
# 
