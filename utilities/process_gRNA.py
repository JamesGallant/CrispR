import os
import pandas as pd

from Bio import SeqIO

class Process:
	def __init__(self, gRNA_files_directory, config_file):
		self.input_dir = gRNA_files_directory
		self.config =pd.read_csv(filepath_or_buffer=config_file, sep="\t")
		self.fasta_file = [str(filepath) for filepath in self.config['input_file']][0]
		
		if "Summary.xls" not in os.listdir(self.input_dir):
			raise ValueError("Summary.xls not found")
		else:
			self.summary = pd.read_csv(filepath_or_buffer=os.path.join(self.input_dir, "Summary.xls"), 
				index_col = False, sep = "\t")

		if "OfftargetAnalysis.xls" in os.listdir(self.input_dir):
			self.offtarget = pd.read_csv(filepath_or_buffer=os.path.join(self.input_dir, "OfftargetAnalysis.xls"), 
				index_col = False, sep = "\t")
		else:
			self.offtarget = None

		self.valid_gRNA = {
		'name': [],
		'gRNA':	[]
		}

		self.invalid_gRNA = {
		'name': [],
		'gRNA': [],
		'notes': []	
		}


	def detected_genes(self):
		detected = self.summary['names'].tolist()
		detected = [ids.split("_")[0] for ids in detected]
		return list(set(detected))
		 
	def no_detection(self):
		"""
		Detects all the genes present within the input fasta file. This is cross referenced
		to a list of deduplicated list of gRNA names. The difference is the genes that did not produce a gRNA and 
		appended to a dictionary with some metadata
		"""
		out = {
		'name': [],
		'gRNA': [],
		'notes': []	
		}

		all_searched_genes = [fasta.id for fasta in SeqIO.parse(open(self.fasta_file), 'fasta')]
		not_detected = list(set(all_searched_genes) - set(self.detected_genes()))

		for genes in not_detected:
			out['name'].append(genes)
			out['gRNA'].append("NA")
			out['notes'].append(f"No gRNA detected based on PAM: {[str(filepath) for filepath in self.config['PAM_sequence']][0]}")

		return out

	def gRNA_first_pass(self):
		"""
		first pass iteration to filter for possible matches. This function will return three dictionaries:
		canidates which have passed constraints, backup candidates that can be accepted in follow up rounds and those that need to be ignored.
		Contstraints are gRNA has to be on the reverse strand and the 20th position is either a G/A
		"""

		candidates = {
		'name':	[],
		'gRNA':	[]
		}

		backup = {
		'name': [],
		'gRNA': [],
		'notes': []
		}

		dropped_gRNA = {
		'name':	[],
		'gRNA':	[],
		'notes':[]
		}

		for _, row in self.summary.iterrows():
			if row['names'][-1] == 'r':
				if row['gRNAsPlusPAM'][20] == "A" or row['gRNAsPlusPAM'][20] == "G":
					candidates['name'].append(row['names'])
					candidates['gRNA'].append(row['gRNAsPlusPAM'])
				else:
					backup['name'].append(row['names'])
					backup['gRNA'].append(row['gRNAsPlusPAM'])
					backup['notes'].append(f"position 20 is {row['gRNAsPlusPAM'][20]}")
			else:
				dropped_gRNA['name'].append(row['names'])
				dropped_gRNA['gRNA'].append(row['gRNAsPlusPAM'])
				dropped_gRNA['notes'].append(f"The gRNA {row['names']} is not on the reverse strand")

		return candidates, backup, dropped_gRNA


	def gRNA_second_pass(self):
		"""
		returns a dictionary

		"""
		candidates, backup, _ = self.gRNA_first_pass()



	def correct_gRNA_dictionary(self):
		"""
		returns a dictionary
		checks if gRNA's are in backup and not featured in the invalids or the valids. In this case the guide RNA from a specific
		gene needs to be moved for further processing. This function assumes that invalids are ground truth and if a gRNA is found 
		in backup and invalid the gRNA will be reassigned in invalid. If it not found in valids then it will be moved to valid 
		which disregards previous filtering criteria to reduce gRNA loss.
		"""
		genes_nogRNA_present = self.no_detection()
		candidates, backup, dropped_gRNA = self.gRNA_first_pass()

		valids = list(set([genes.split("_")[0] for genes in candidates['name']]))
		backup_genes = list(set([genes.split("_")[0] for genes in backup['name']]))

		invalid1 = list(set([genes.split("_")[0] for genes in genes_nogRNA_present['name']]))
		invalid2 = list(set([genes.split("_")[0] for genes in dropped_gRNA['name']]))

		invalids = list(set(invalid1+invalid2) - set(valids))

		
		#backup_genes = ["Rv0001", 'Rv3868', 'Rv3875']
		
		for genes in backup_genes:
			if genes not in valids:
				if genes not in invalids:
					print(f"{genes} move to candidates")
				else:
					print(f"{genes} move to invalid")
			else:
				print(f"{genes} stay as a backup")







