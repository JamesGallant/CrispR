import os
import unittest
import pandas as pd
import shutil

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
	'input_file': "path to fasta file",
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
	known_files = ["inputallgRNAs.fa", "OfftargetAnalysis.xls", "pairedgRNAs.xls", "REcutDetails.xls",
	"Rv0973c.gbk", "Summary.xls"]
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

	if "temp" not in os.listdir(root):
		os.mkdir(os.path.join(root, "temp"))
		tempdir = os.path.join(root, "temp")

	grna.create_gRNA_library(out_dir=tempdir, 
	config_file=os.path.join(os.getcwd(), "references", "gRNA", "config.txt"))

	generated_files = os.listdir(tempdir)

	shutil.rmtree(tempdir)

	if set(known_files) == set(generated_files):
		return True
	else:
		return False


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

	