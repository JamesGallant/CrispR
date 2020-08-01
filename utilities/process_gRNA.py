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




	def correct_gRNA_dictionary(self, canidates, backup, dropped_gRNA):
		"""
		returns a mulitple dictionary
		checks if gRNA's are in backup and not featured in the valids. In this case the guide RNA from a specific
		gene needs to be moved for further processing. If it not found in valids then it will be moved to valid 
		which disregards previous filtering criteria to reduce gRNA loss.
		"""
		def _split(input_dict, splitter, idx):
			"""
			Takes a dictionary, returns a list
			"""
			return list(set([genes.split(splitter)[idx] for genes in input_dict['name']]))

		def _get_index(dataframe, identifier):
			"""
			takes a dictionary and a string. The identifier is a substring that will be matched to the main string in the 
			dataframe. retruns the index.
			"""
			for idx, value in backup.iterrows():
				if str(value['name']).find(identifier) == 0:
					return idx

		def _check_altered_dict(input_dict, identifier):
			"""
			Error handling step, if items are erroneously moved it will be caught here
			"""
			for _, items in input_dict.items():
				if str(items.values()).find(identifier) != -1:
					raise ValueError(f"Critical error, dictionary correction failed. identifier {identifier} is present in wrong dictionary")


		#candidates, backup, dropped_gRNA = self.gRNA_first_pass()

		valids = _split(input_dict=candidates, splitter="_", idx=0)
		backup_genes = _split(input_dict=backup, splitter="_", idx=0)

		genes_nogRNA_present = self.no_detection()
		invalid1 = list(set([genes.split("_")[0] for genes in genes_nogRNA_present['name']]))
		invalid2 = list(set([genes.split("_")[0] for genes in dropped_gRNA['name']]))
		invalids = list(set(invalid1+invalid2) - set(valids))

		for genes in backup_genes:
			if genes not in valids:
				if genes not in invalids:
					#moves backup genes back to candidates
					backup = pd.DataFrame(backup_dict)
					backup_row_idx = _get_index(dataframe=backup, identifier=genes)
					candidates['name'].append(backup.iloc[backup_row_idx, 0])
					candidates['gRNA'].append(backup.iloc[backup_row_idx, 1])
					backup.drop(backup.index[[1]], inplace=True)
					backup = backup.to_dict()
					_check_altered_dict(input_dict=backup, identifier=genes)
				else:
					#Moves guide RNA's between backup and invalid think if this makes sense
					pass

		return candidates, backup, dropped_gRNA


	def check_offtarget_effects(self):
		"""
		In this function the canidate gRNA's will be cross referenced with the off target effects. candidates will be moved between
		backup and valids iteratively
		"""
		pass







