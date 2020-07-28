import os
import unittest
import pandas as pd
import shutil
import sqlite3

from pandas.testing import assert_frame_equal
from utilities.predict_gRNA import Get_gRNA

grna = Get_gRNA()

def create_config_file():
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

	if "temp" not in os.listdir(root):
		os.mkdir(os.path.join(root, "temp"))
		tempdir = os.path.join(root, "temp")

	grna.generate_config(out_dir = tempdir)

	dataframe = pd.read_csv(filepath_or_buffer=os.path.join(tempdir, "config.txt"),
		sep="\t", index_col=False)
	
	shutil.rmtree(tempdir)

	return dataframe

def valid_config_file():
	df = {
	'organism': "Organism name",
	'circular_chromosome': "logical TRUE/FALSE",
	'input_file': "path to input fasta file",
	'gff_file': "Path to gff file",
	'find_gRNA_with_cutsites' : "logical TRUE/FALSE",
    'find_paired_gRNA' : "logical TRUE/FALSE",
    'BSgenome' : "Name of the BSgenome object (BSgenome.Mtuberculosis.VUmc.H37Rv for H37Rv)",
    'chromosomes_to_search' : "all",
    'min_gap' : 0,
    'max_gap' : 20, 
    'max_mismatch_gRNA' : 4,
    'PAM_sequence' : "AGAA",
    'gRNA_size' : 20,
    'scoring_method' : "CFDscore/Hsu-Zhang"
	}

	return pd.DataFrame(df, index=[0])

def files_created():
	known_files = ["input_fastaallgRNAs.fa", "OfftargetAnalysis.xls", "pairedgRNAs.xls", "REcutDetails.xls",
	"Rv0973c.gbk", "Summary.xls", "config.txt", "input_fasta.fa", "test.gff"]

	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

	gff_file = pd.read_sql(sql="SELECT * FROM gff_file", con=sqlite3.connect(os.path.join(root, "databases", "test.db")))

	if "temp" not in os.listdir(root):
		os.mkdir(os.path.join(root, "temp"))
		tempdir = os.path.join(root, "temp")


	gff_file.to_csv(path_or_buf = os.path.join(tempdir, "test.gff"), 
		header=False, index=False, sep="\t")

	querry_fasta = open(os.path.join(tempdir, "input_fasta.fa"), "a")
	querry_fasta.write(">Rv0973c\n")
	querry_fasta.write("atgggaatcactcgagtattggttgctaaccgcggcgagatcgcccggcgggtgttcgccacctgccgccggctggggctcggcaccgtcgccgtctacacag\n")
	querry_fasta.write("acccggatgccgcggcaccgcatgtcgccgaggccgacgcccgggtccggctgccgcagaccaccgactatctgaacgccgaggcgatcatcgcggccgcgcagg\n")
	querry_fasta.write("agccggagccgacgcggtgcatcccggctacggattcctctcggagaacgccgaattcgcggccgccgtgcaggaggccggcctaacctgggtcgggccaccggt\n")
	querry_fasta.write("gcagctgccggtgctggtgaaggcgtcggcgggcggtggcggtcgcggcatgcgagtggttcacgaattatcggccctgccggccgaagtcgaagccgcgcgacgt\n")
	querry_fasta.write("gaagcccaatccgcgttcggcgacccgaccgtattctgcgagcgctacctgcccaccgggcaccacgtcgaagtgcaagtcatggccgacacccatggcaccgtgt\n")
	querry_fasta.write("gggcggtcggggaacgggaatgctcgattcagcgccgccaccagaagatcatcgaagaggcaccgtcgccgctggtggaacgcgtaccggggatgcgggccaagct\n")
	querry_fasta.write("gttcgacgcggcccggctggcggccagcgcgatcggctacaccggggcaggcacggtggagttcctcgccgatgactcacctggccgggaaggtgagttctacttc\n")
	querry_fasta.write("ctggagatgaacacccggctacaggtcgagcacccggtcaccgaagagaccaccgggctggatctggtcgaactgcagctcatgattgccgactgcgggcgactcg\n")
	querry_fasta.write("acaccgaacctccccccgcccagggatattcgatcgaggcccgcctctacgccgaggatcccgcgcatggctggcagccacaggcaggcgtgatgcacacgattga\n")
	querry_fasta.write("ggttccgggggttcgggcgcagttcgactcgttgggacagcggaccggcatccggctggattccgggatcgtcgacggttccacagtgtcgatccactacgaccca\n")
	querry_fasta.write("atgctggccaaggtcgtctcctacggtgccacccgccggcaggccgcgcttgtgttggccgatgcgctggtacgcgcccggctgcacggtctgcgcaccaaccgtg\n")
	querry_fasta.write("agctcttggtcaacgtgctgcgtcatccggcgttcctcgacggcgccaccgacaccgggtttttcgacacgcacggcatggccgagttgtcgacaccgctggccga\n")
	querry_fasta.write("caccgcgaccctccggttgtcggcgatcgccgccgcgctggccgacgccgagcacaatcgggcgagcgcgggcgtgttcagctcgattcccagcggctggcgcaac\n")
	querry_fasta.write("ctggcctcgggctatcaggtcaagacctatcgtgacgacgcggacaccgaacaccgcgtcgaataccggttcaccagaacgggtctggcgcttcccggcgatccgg\n")
	querry_fasta.write("tggtacagctggtctcggctgacgtggaccaggtggtgctcgcccaggacggggtcgcacacggcttcacggttgcccgccacggccccgacgtctacgtcgactc\n")
	querry_fasta.write("ggcgcgcggacccgttcacctggtggcactgtcacgcttccccgagccgagctcggccgtcgagcaaggctcgctggtggcccccatgcccggcaacgtcatccgg\n")
	querry_fasta.write("atcggcgccgaggttggcgacacggtcacggccggtcagccgttgatctggctggaggccatgaagatggaacacaccatcgccgcgcctgccgacggcgtgctca\n")
	querry_fasta.write("cccacgtcagcgtcaacacgggtcaacaggtcgaagtaggcgccattctcgcacgagtagaagcacctcaaaatggcccagcagaaggagattcaccgtga")
	querry_fasta.close()

	querry_fasta_path = os.path.join(tempdir, "input_fasta.fa")
	testing_gff_path = os.path.join(tempdir, "test.gff")

	config = open(os.path.join(root, tempdir, "config.txt"), "a")
	config.write("organism\tcircular_chromosome\tinput_file\tgff_file\tfind_gRNA_with_cutsites\tfind_paired_gRNA\tBSgenome\tchromosomes_to_search\tmin_gap\tmax_gap\tmax_mismatch_gRNA\tPAM_sequence\tgRNA_size\tscoring_method\n")
	config.write(f"Mycobacterium tuberculosis\tTRUE\t{querry_fasta_path}\t{testing_gff_path}\tTRUE\tFALSE\tBSgenome.Mtuberculosis.VUmc.H37Rv\tall\t0\t20\t4\tAGAA\t20\tHsu-Zhang\n")
	config.close()

	grna.create_gRNA_library(out_dir=tempdir, 
	config_file=os.path.join(tempdir, "config.txt"))

	generated_files = os.listdir(tempdir)

	shutil.rmtree(tempdir)
	os.remove(os.path.join(root, "offTargets.RDS"))
	os.remove(os.path.join(root, "scores.RDS"))

	if set(known_files) == set(generated_files):
		return True
	else:
		return str("Files are missing")


class Predict_gRNA_test(unittest.TestCase):
	def test_config_file(self):
		try:
			assert_frame_equal(create_config_file(), valid_config_file())
		except AssertionError as e:
			raise self.failureException(e)

	def test_gRNA_files(self):
		self.assertEqual(True, files_created())

if __name__ == '__main__':
    unittest.main()

	