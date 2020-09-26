import pandas as pd

class RefineCripri:
	"""
	Dataframe is a pandas dataframe
	"""
	def __init__(self, grna_dataframe, strand, fasta_dataframe):
		self.grna_dataframe = grna_dataframe
		self.strand = strand
		self.fasta = fasta_dataframe


	def detected_genes(self):
		detected = self.grna_dataframe['names'].tolist()
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
		'notes': [],
		'PAM': []
		}

		all_searched_genes = [header for header in self.fasta['header']]
		not_detected = list(set(all_searched_genes) - set(self.detected_genes()))
		for genes in not_detected:
			out['name'].append(genes)
			out['gRNA'].append("NA")
			out['PAM'].append("NA")
			out['notes'].append(f"No gRNA detected based on PAM: {set(self.grna_dataframe['PAM'])}")

		return out


	def initial_filter(self):
		"""
		first pass iteration to filter for possible matches. This function will return three dictionaries:
		canidates which have passed constraints, backup candidates that can be accepted in follow up rounds and those that need to be ignored.
		Contstraints are gRNA has to be on the reverse strand and the 20th position is either a G/A
		"""

		candidates = {
		'name':	[],
		'gRNA':	[],
		'PAM': [],
		'notes': []
		}

		backup = {
		'name': [],
		'gRNA': [],
		'PAM': [],
		'notes': []
		}

		dropped_gRNA = {
		'name':	[],
		'gRNA':	[],
		'PAM': [],
		'notes':[]
		}

		for _, row in self.grna_dataframe.iterrows():
			if row['names'][-1] == self.strand:
				if row['gRNAsPlusPAM'][0] == "A" or row['gRNAsPlusPAM'][0] == "G":
					candidates['name'].append(row['names'])
					candidates['gRNA'].append(row['gRNAsPlusPAM'])
					candidates['PAM'].append(row['PAM'])
					candidates['notes'].append("PASS")
				else:
					backup['name'].append(row['names'])
					backup['gRNA'].append(row['gRNAsPlusPAM'])
					backup['PAM'].append(row['PAM'])
					backup['notes'].append(f"position 20 is {row['gRNAsPlusPAM'][0]}")
			else:
				dropped_gRNA['name'].append(row['names'])
				dropped_gRNA['gRNA'].append(row['gRNAsPlusPAM'])
				dropped_gRNA['PAM'].append(row['PAM'])
				dropped_gRNA['notes'].append(f"The gRNA {row['names']} is not on the reverse strand")

		return candidates, backup, dropped_gRNA


	def correct_gRNA_dictionary(self, candidates, backup, dropped_gRNA):
		"""
		returns a mulitple dictionary
		checks if gRNA's are in backup and not featured in the valids. In this case the guide RNA from a specific
		gene needs to be moved for further processing. If it not found in valids then it will be moved to valid 
		which disregards previous filtering criteria to reduce gRNA loss.
		"""
		backup2 = {
		'name': [],
		'gRNA': [],
		'PAM': [],
		'notes': []
		}

		candidate_genes = [names.split("_")[0] for names in candidates['name']]
		backup_genes = [names.split("_")[0] for names in backup['name']]
		genes_to_move = set(backup_genes) - set(candidate_genes)

		for genes in genes_to_move:
			for idx, items in enumerate(backup['name']):
				if genes in items:
					candidates['name'].append(backup['name'][idx])
					candidates['gRNA'].append(backup["gRNA"][idx])
					candidates['PAM'].append(backup["PAM"][idx])
					candidates['notes'].append(backup["notes"][idx])	


		candidate_genes_new = [names for names in candidates['name']]
		for idx, value in enumerate(backup['name']):
			if value not in candidate_genes_new:
				backup2['name'].append(backup['name'][idx])
				backup2['gRNA'].append(backup['gRNA'][idx])
				backup2['PAM'].append(backup['PAM'][idx])
				backup2['notes'].append(backup['notes'][idx])

		return candidates, backup2, dropped_gRNA


	def has_offtarget(self, candidates, backup, dropped_gRNA):
		"""
		Takes the three input dictionaries, performs the function and returns them again. This function will check if there are
		off targets present, if not it will be moved to candidates. If off targets are present it will be moved to backup for
		follow up. Executes a correction before returning the dataframe
		"""
		for identifiers in candidates["name"]:
			cross_reference = self.grna_dataframe.loc[self.grna_dataframe['names'] == identifiers].to_dict()
			off_target_presence = str(list(cross_reference['top5OfftargetTotalScore'].values())[0])
			if off_target_presence != "nan":
				backup['name'].append(identifiers)
				backup['gRNA'].append(str(list(cross_reference['gRNAsPlusPAM'].values())[0]))
				backup['PAM'].append(str(list(cross_reference['PAM'].values())[0]))
				backup['notes'].append("Has off target")
				candidates = pd.DataFrame(candidates)
				candidates = candidates[candidates['name'] != identifiers]
				candidates.reset_index(drop=True, inplace=True)
				candidates = candidates.to_dict()

		candidates['name'] = list(candidates['name'].values())
		candidates['gRNA'] = list(candidates['gRNA'].values())
		candidates['PAM'] = list(candidates['PAM'].values())
		candidates['notes'] = list(candidates['notes'].values())

		candidates, backup, dropped_gRNA = self.correct_gRNA_dictionary(candidates=candidates, backup=backup, dropped_gRNA=dropped_gRNA)

		return candidates, backup, dropped_gRNA


	def cripr_interference(self):
		candidates, backup, dropped = self.has_offtarget(candidates=self.correct_gRNA_dictionary(candidates=self.initial_filter()[0], 
			backup=self.initial_filter()[1], 
			dropped_gRNA=self.initial_filter()[2])[0],
			backup=self.correct_gRNA_dictionary(candidates=self.initial_filter()[0], 
			backup=self.initial_filter()[1], 
			dropped_gRNA=self.initial_filter()[2])[1],
			dropped_gRNA=self.correct_gRNA_dictionary(candidates=self.initial_filter()[0], 
			backup=self.initial_filter()[1], 
			dropped_gRNA=self.initial_filter()[2])[2])
		
		return candidates, backup, dropped

