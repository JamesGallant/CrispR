import os
import sqlite3
from datetime import date
import pandas as pd
from utilities.database_tools import Database

class CreateDatabaseMethod:
	def __init__(self, filepaths_dict, metadata_dict):
		self.root = os.path.dirname(os.path.abspath(__name__))
		self.tempdir = os.path.join(self.root, "temp")
		self.genome_fasta = filepaths_dict['genome']
		self.gene_fasta = filepaths_dict['multifasta']
		self.gff = filepaths_dict['gff']
		self.author = metadata_dict['author']
		self.department = metadata_dict['department']
		self.org_name = metadata_dict['org_name']
		self.common_name = metadata_dict['common_name']

	def create(self):
		version = str(date.today()).split("-")[0:2]
		version = "".join(version[0] + "." + version[1])

		genus, species, strain = self.org_name.split(" ")
		species = species.lower()

		if " " in self.department:
			self.department = self.department.replace(" ", "_")

		common_name = self.common_name
		if " " in common_name:
			common_name = common_name.replace(" ", "_")

		database_name = str(self.org_name).replace(" ", "_")

		db = Database(database=database_name)
		db.create_new_database(gff_file=self.gff, 
			genome=self.genome_fasta, 
			multifasta=self.gene_fasta)

		tempfa = pd.read_sql("SELECT * FROM genome", sqlite3.connect(os.path.join(os.getcwd(), "databases", f"{database_name}.db")))
		fasta_header = [header for header in tempfa['header']]
		sequence = [sequence for sequence in tempfa['sequence']]

		seqnames = ', '.join('"{0}"'.format(word) for word in fasta_header)
		seqnames = f'paste(c({seqnames}), sep="")'

		for header, seq in zip(fasta_header, sequence):

			with open(os.path.join(self.root, "temp", "fastafiles", f"{header}.fa"), 'w') as fasta:
				fasta.write(">" + header + "\n")
				fasta.write(seq)
				fasta.close()

		dcf_file = {
			'Package':	f'BSgenome.{common_name}.{self.department}.{strain}',
			'Title':	f'Genome sequence of {self.org_name}',	
			'Description':	f'BSgenome genome obtained from the genome sequence of {self.org_name}',
			'Version':	version,
			'Author':	self.author,
			'License':	"GPLv2",
			'organism':	self.org_name,
			'common_name':	common_name,
			'provider':	self.author,
			'provider_version':	f'{strain}.{str(date.today())}',
			'release_date':	str(date.today()),
			'release_name':	f'{self.department} {self.author} {strain}',
			'organism_biocview':	f'{genus}_{species}',
			'BSgenomeObjname':	common_name,
			'seqnames':	seqnames,
			'seqs_srcdir':	os.path.join(self.tempdir, "fastafiles")

			}

		with open(os.path.join(self.tempdir, "dcf_file.dcf"), 'w') as dcf:
			for k, v in dcf_file.items():
				dcf.write(f"{k}: {v} \n")





		