import os
import pandas as pd
from collections import Counter

from PyQt5 import QtCore
from utilities.build_reference import BSgenome
from utilities.predict_gRNA import Get_gRNA
from utilities.annotate_off_target import Annotater
from utilities.cripr import RefineCripri
from utilities.database_tools import SQL


class BSgenome_worker(QtCore.QRunnable):
	def __init__(self):
		super().__init__()
		self.bs = BSgenome()
		self.root = os.path.dirname(os.path.abspath(__name__))

	@ QtCore.pyqtSlot()
	def run(self):
		self.bs.bsgenome_from_seed(seedfile=os.path.join(self.root, "temp", "dcf_file.dcf"))


class FindgRNA_worker(QtCore.QRunnable):
	def __init__(self):
		super().__init__()
		self.root = os.path.dirname(os.path.abspath(__name__))
		self.tempdir = os.path.join(self.root, "temp")
		self.global_gRNA = os.path.join(self.root, "temp", "global_gRNA")
		self.config_file = os.path.join(self.tempdir, "config.txt")
		self.offtarget_file = os.path.join(self.global_gRNA, "OfftargetAnalysis.xls")

		self.gRNA = Get_gRNA()


	@ QtCore.pyqtSlot()
	def run(self):
		"""
		Creates a gRNA dataset from the CrisprSeek library in R, then adds it to the organism DB
		"""
		self.gRNA.create_gRNA_library(config_file=self.config_file,
			out_dir = self.global_gRNA)


		annotater = Annotater(config_file=pd.read_csv(self.config_file, sep="\t"),
		 offtarget_file=pd.read_csv(self.offtarget_file, index_col=False, sep="\t"))
		offtarget_file = annotater.annotate()

		offtarget_file.to_csv(path_or_buf=os.path.join(self.global_gRNA, "OfftargetAnalysis.txt"), index=False, header=True, sep="\t")


class CrisprInterference_worker(QtCore.QRunnable):
	def __init__(self, database, mismatch, strand, max_grna, genes_masks, max_primer_size):
		super().__init__()
		self.root = os.path.dirname(os.path.abspath(__name__))
		self.tempdir = os.path.join(self.root, "temp")
		self.database = database
		self.mismatch = mismatch
		self.strand = strand
		self.max_grna = max_grna
		self.gene_mask_dict = genes_masks
		self.max_primer_size = max_primer_size


	def scan_maxmismatches(self, candidates, backup):
		"""
		Function takes dictionaries and checks if the number of guide RNAs if more than the max it removes some candidates
		Needs to consider A/G and PAMS preferentially Does not consider PAMS atm
		"""
		candidates2 = {
			'name': [],
			'gRNA': [],
			'PAM': [],
			'notes': []
		}

		number_of_genes = [genes.split("_")[0] for genes in candidates['name']]
		out = Counter(number_of_genes)

		for genes, gene_counts in out.items():
			get_gene_idx = [idx for idx, value in enumerate(candidates['name']) if genes in value]
			if gene_counts > self.max_grna:
				# count the number of genes to move
				no_genes_to_remove = len(get_gene_idx) - self.max_grna
				# get the index that correpsonds to gRNA's with either an A/G in all samples
				get_at_idx = [idx for idx, value in enumerate(candidates['notes']) if "Position 20 is" in value]
				# get the index of A/G that is specific to our target genes
				genes_that_has_at = list(set(get_gene_idx) & set(get_at_idx))
				# Find the genes that pass as well
				genes_that_pass = list(set(get_gene_idx) - set(get_at_idx))
				# Create a new list of indices to that will be used to pop and append
				indices_to_move = genes_that_has_at 
				if len(indices_to_move) > no_genes_to_remove:
					for remove_items in range(no_genes_to_remove):
						indices_to_move.pop(remove_items)
				else:
					move_grna_count = no_genes_to_remove - len(genes_that_has_at)
					for add_items in range(move_grna_count):
						indices_to_move.append(genes_that_pass[add_items])


				indices_to_keep = list(set(get_gene_idx) - set(indices_to_move))

				# move to backup
				for geneitems in indices_to_move:
					backup['name'].append(candidates['name'][geneitems])
					backup['gRNA'].append(candidates['gRNA'][geneitems])
					backup['PAM'].append(candidates['PAM'][geneitems])
					backup['notes'].append(candidates['notes'][geneitems])


				for geneitems in indices_to_keep:
					candidates2['name'].append(candidates['name'][geneitems])
					candidates2['gRNA'].append(candidates['gRNA'][geneitems])
					candidates2['PAM'].append(candidates['PAM'][geneitems])
					candidates2['notes'].append(candidates['notes'][geneitems])

			else:
				for items in get_gene_idx:
					candidates2['name'].append(candidates['name'][items])
					candidates2['gRNA'].append(candidates['gRNA'][items])
					candidates2['PAM'].append(candidates['PAM'][items])
					candidates2['notes'].append(candidates['notes'][items])

		return candidates2, backup

	def grab_targets(self, candidates, backup, dropped, gene_mask_dict):
		"""
		gets the genes correpsonding to the targets as given by the user returns the altered dictionaries
		"""
		user_candidates = {
			'name': [],
			'gRNA': [],
			'PAM': [],
			'notes': []
		}

		user_backup = {
			'name': [],
			'gRNA': [],
			'PAM': [],
			'notes': []
		}

		user_dropped = {
			'name': [],
			'gRNA': [],
			'PAM': [],
			'notes': []
		}


		for genes in gene_mask_dict['genes']:
			get_candidate_idx = [idx for idx, value in enumerate(candidates['name']) if genes in value]
			get_backup_idx = [idx for idx, value in enumerate(backup['name']) if genes in value]
			get_dropped_idx = [idx for idx, value in enumerate(dropped['name']) if genes in value]

			for items in get_candidate_idx:
				user_candidates['name'].append(candidates['name'][items])
				user_candidates['gRNA'].append(candidates['gRNA'][items])
				user_candidates['PAM'].append(candidates['PAM'][items])
				user_candidates['notes'].append(candidates['notes'][items])		

			for items in get_backup_idx:
				user_backup['name'].append(backup['name'][items])
				user_backup['gRNA'].append(backup['gRNA'][items])
				user_backup['PAM'].append(backup['PAM'][items])
				user_backup['notes'].append(backup['notes'][items])

			for items in get_dropped_idx:
				user_dropped['name'].append(dropped['name'][items])
				user_dropped['gRNA'].append(dropped['gRNA'][items])
				user_dropped['PAM'].append(dropped['PAM'][items])
				user_dropped['notes'].append(dropped['notes'][items])

		return user_candidates, user_backup, user_dropped



	def grab_offtargets(self, query, offtargets, offtarget_ids):
		"""
		query corresponds to the candidates/backup dicts and off targets is the global offtarget files. Will return a 
		pandas dataframe
		"""
		out = pd.DataFrame()

		for items in query['name']:
			if items in offtarget_ids:
				out = out.append(offtargets[offtargets['name'].str.match(items)], ignore_index=True)

		return out

	def list_comparison(self, list1, list2):
		"""
		check if list1 has a element present in list 2

		""" 
		set1 = set(list1) 
		set2 = set(list2) 
		if set1.intersection(set2): 
			return True 
		else: 
			return False


	def move_grna_by_offtargets(self, grna_dataframe, dropped_dataframe, offtarget_dataframe, masks):
		"""
		This function will do the masking, if a grna is present in grna_dataframe['name'] that is also present in gene_mask_dict['masks']
		and that can be found in the offtarget database based on the function grab_offtargets it will be moved to dropped
		Returns:
		new dropped dataframe, new grna dataframe
		use for both canididates and backups
		"""
		temp = {
		'name': [],
		'gRNA': [],
		'PAM': [],
		'notes': []
		}

		# gets the indices where offtarget genes are located in masking list from user, map to grna names to cross reference later and dedup
		drop_offtarget_idx = [idx for idx, value in enumerate(offtarget_dataframe['annotation']) if value in masks]
		grna_to_move = [offtarget_dataframe['name'][i] for i in drop_offtarget_idx]
		grna_to_move = list(set(grna_to_move))

		#map grna to grna dataframe and find incices where they coincide
		grna_idx = [idx for idx, value in enumerate(grna_dataframe['name']) if value in grna_to_move]

		# append the information to dropped data
		for grna_items in grna_idx:
			dropped_dataframe['name'].append(grna_dataframe['name'][grna_items])
			dropped_dataframe['gRNA'].append(grna_dataframe['gRNA'][grna_items])
			dropped_dataframe['PAM'].append(grna_dataframe['PAM'][grna_items])
			dropped_dataframe['notes'].append("Has off target that was blocked")


		# checks if its neccesary to remake the candidates and does so
		if len(grna_dataframe['name']) > len(grna_idx):
			remake_idx = [idx for idx, value in enumerate(grna_dataframe['name']) if idx not in grna_idx]
			for remake_index in remake_idx:
				temp['name'].append(grna_dataframe['name'][remake_index])
				temp['gRNA'].append(grna_dataframe['gRNA'][remake_index])
				temp['PAM'].append(grna_dataframe['PAM'][remake_index])
				temp['notes'].append(grna_dataframe['notes'][remake_index])
		else:
			temp = dict.fromkeys(grna_dataframe, [])

		return temp, dropped_dataframe


	def force_max_grna_in_candidates(self, candidates, backup, max_grna):
		"""
		Assumes that all backup genes are valid and can be used, this must run AFTER masking algorithms
		This algorithm does not discriminate between offtargets and position 20 grna's in backup
		Takes in a max grna from user as well as candidates and backup in either dict of pandas dataframe form
		"""

		# This algorithm does not discriminate between offtargets and position 20 grna's in backup
		candidates_genes_count = Counter(candidates['genes'])
		backup_genes_count = Counter(backup['genes'])
		backup_drop_idx = []
		for genes, counts in candidates_genes_count.items():
			if counts < max_grna:
				genes_in_backup = backup_genes_count.get(genes, 0)
				if genes_in_backup > 0:
					max_number_to_move = max_grna - counts
					moved = 0
					moving_idx = [idx for idx, value in enumerate(backup['genes']) if value in genes]
					if len(moving_idx) > max_number_to_move:
						moving_idx_counter = 0
						while moving_idx_counter < max_number_to_move:
							for backup_indices in moving_idx:
								candidates['name'].append(backup['name'][backup_indices])
								candidates['gRNA'].append(backup['gRNA'][backup_indices])
								candidates['PAM'].append(backup['PAM'][backup_indices])
								candidates['notes'].append(backup['notes'][backup_indices])
								candidates['genes'].append(backup['genes'][backup_indices])
								backup_drop_idx.append(backup_indices)
								moving_idx_counter += 1


					if len(moving_idx) > 0 and len(moving_idx) < max_number_to_move:
						for backup_indices in moving_idx:
							candidates['name'].append(backup['name'][backup_indices])
							candidates['gRNA'].append(backup['gRNA'][backup_indices])
							candidates['PAM'].append(backup['PAM'][backup_indices])
							candidates['notes'].append(backup['notes'][backup_indices])
							candidates['genes'].append(backup['genes'][backup_indices])
							backup_drop_idx.append(backup_indices)


		backup = pd.DataFrame(backup)
		backup.drop(backup.index[backup_drop_idx], inplace=True)
		backup.reset_index(drop=True, inplace=True)
		backup = backup.to_dict('list')
		return candidates, backup


	def force_ag_base(self, dataframe, max_primer_size):
		"""
		Algorithm matches the grna to the gene and increments the bases if they are not a/g this modifies the index of gRNA in place
		"""
		sqlrunner = SQL(database=os.path.join(self.root, "databases", self.database))
		for idx, genes in enumerate(dataframe['genes']):
			gene_sequence = str(sqlrunner.get_gene_sequence(gene=genes)).lower()
		
			
			complement_strand_dict = {
			'g': 'c',
			'G': 'C',
			'a': 't',
			'A': 'T',
			'c': 'g',
			'C': 'G',
			't': 'a',
			'T': 'A'
			}

			grna = str(dataframe['gRNA'][idx]).lower()
			fails_contstraint = False if grna[0] == 'a' or grna[0] == 'g' else True

			sequence_swapped = ""

			if fails_contstraint:
				for base in gene_sequence:
					sequence_swapped += complement_strand_dict.get(base, "N")

				grna = grna[::-1]
				pam_len = len(str(dataframe['PAM'][idx]).lower())
				loc_in_gene = sequence_swapped.find(grna)
				primer_wo_pam_start = loc_in_gene + pam_len
				primer_wo_pam_stop = primer_wo_pam_start + 20 #20 is the primer length, fixed value

				
				if primer_wo_pam_stop + max_primer_size < len(gene_sequence):
					counter = 1
					while counter <= max_primer_size:
						target_base_location = primer_wo_pam_stop + counter
						target_base = sequence_swapped[target_base_location]
						if target_base == "a" or target_base == "g":
							grna_out = sequence_swapped[loc_in_gene:target_base_location]
							dataframe['gRNA'][idx] = grna_out.upper()

							if "position 20 is" in dataframe['notes'][idx]:
								dataframe['notes'][idx] = "PASS"

							break
						else:
							counter += 1
				else:
					counter = 1
					while counter < len(gene_sequence):
						target_base_location = primer_wo_pam_stop + counter
						target_base = sequence_swapped[target_base_location]
						if target_base == "a" or target_base == "g":
							grna_out = sequence_swapped[loc_in_gene:target_base_location]
							dataframe['gRNA'][idx] = grna_out.upper()

							if "position 20 is" in dataframe['notes'][idx]:
								dataframe['notes'][idx] = "PASS"
							break
						else:
							counter += 1

		return dataframe


	def calculate_gc_content(self, dataframe):
		"""
		Calculates the total GC content given a dataframe, returns the dataframe with an extra column named gc_content
		"""
		dataframe = pd.DataFrame(dataframe)
		gc_content_list = []

		for seqs in dataframe['gRNA']:
			seqs_normalised = seqs.lower()
			basecount = Counter(seqs_normalised)
			guanine = basecount['g']
			cytosine = basecount['c']
			gc_content = str(round(((guanine + cytosine)/len(seqs)) * 100,3))
			gc_content_list.append(gc_content)

		dataframe['gc_content'] = gc_content_list

		dataframe.to_dict()
		return dataframe


	def calculate_primer_len(self, dataframe):
		"""
		Calculates the length of a primer
		"""
		dataframe = pd.DataFrame(dataframe)
		seqlen_list = []
		for sequence in dataframe['gRNA']:
			seqlen_list.append(len(sequence))

		dataframe['primer_length'] = seqlen_list

		dataframe.to_dict()
		return dataframe



	@ QtCore.pyqtSlot()
	def run(self):
		sqlrunner = SQL(database=os.path.join(self.root, "databases", self.database))
		gRNA_db = sqlrunner.get_global_gRNA(mismatch=str(self.mismatch))
		multifasta = sqlrunner.get_gene_multifasta()

		gRNA_runner = RefineCripri(grna_dataframe=gRNA_db,
			strand=self.strand, 
			fasta_dataframe=multifasta)


		candidates, backup, dropped = gRNA_runner.cripr_interference()

		if bool(self.gene_mask_dict['genes']):
			candidates, backup, dropped = self.grab_targets(candidates=candidates, 
				backup=backup, 
				dropped=dropped, 
				gene_mask_dict=self.gene_mask_dict)


		candidates, backup = self.scan_maxmismatches(candidates=candidates, backup=backup)

		offtargets = sqlrunner.get_offtargets_by_mismatch(mismatch=self.mismatch)
		#assuming we can tolerate mismatches on +
		offtargets.dropna(subset = ['annotation'], inplace=True)
		offtargets = offtargets[offtargets['strand'] != '+']
		offtargets['annotation'] = offtargets['annotation'].apply(lambda x: x.replace("_", "") if isinstance(x, str) else x)
		offtargets = offtargets.query("gene != annotation")
		offtargets.reset_index(drop=True, inplace=True)
		#Evaluate if we need to calculate off targets
		offtarget_ids = list(set(offtargets['name']))
		candidates_has_offtargets = self.list_comparison(list1=candidates['name'], list2=offtarget_ids)
		backup_has_offtargets = self.list_comparison(list1=backup['name'], list2=offtarget_ids)

		if candidates_has_offtargets:
			candidates_offtargets = self.grab_offtargets(query=candidates, 
				offtargets=offtargets, 
				offtarget_ids=offtarget_ids)

			candidates, dropped = self.move_grna_by_offtargets(grna_dataframe=candidates, 
				dropped_dataframe=dropped, 
				offtarget_dataframe=candidates_offtargets,
				masks = self.gene_mask_dict['masks'])
		else:
			offtargets = offtargets.to_dict()
			candidates_offtargets = dict.fromkeys(offtargets, [])
			offtargets = pd.DataFrame(offtargets)


		if backup_has_offtargets:
			backup_offtargets = self.grab_offtargets(query=backup, offtargets=offtargets, offtarget_ids=offtarget_ids)

			backup, dropped = self.move_grna_by_offtargets(grna_dataframe=backup, 
				dropped_dataframe=dropped, offtarget_dataframe=backup_offtargets,
				masks=self.gene_mask_dict['masks'])

		else:
			offtargets = offtargets.to_dict()
			backup_offtargets = dict.fromkeys(offtargets, [])
			offtargets = pd.DataFrame(offtargets)

		# move from backup back to candidates here if we want to force more, other function will only work on "PASS" grna's
		candidates['genes'] = [genes.split("_")[0] for genes in candidates['name']]
		backup['genes'] = [genes.split("_")[0] for genes in backup['name']]

		candidates, backup = self.force_max_grna_in_candidates(candidates=candidates, backup=backup, max_grna=self.max_grna)

		candidates = self.force_ag_base(dataframe=candidates, max_primer_size=self.max_primer_size)

		backup = self.force_ag_base(dataframe=backup, max_primer_size=self.max_primer_size)

		candidates, backup, dropped = map(self.calculate_primer_len, [candidates, backup, dropped])

		candidates, backup, dropped = map(self.calculate_gc_content, [candidates, backup, dropped])

		candidates, backup, dropped, candidates_offtargets, backup_offtargets = map(pd.DataFrame,
		 [candidates, backup, dropped, candidates_offtargets, backup_offtargets])
		
		offtargets = candidates_offtargets
		offtargets['from'] = "candidates"
		backup_offtargets['from'] = "backup"
		offtargets.append(backup_offtargets, ignore_index=True)
		candidates.to_csv(os.path.join(self.root, "temp", "candidates.txt"), header=True, index=False, sep=",")
		backup.to_csv(os.path.join(self.root, "temp", "backup.txt"), header=True, index=False, sep=",")
		dropped.to_csv(os.path.join(self.root, "temp", "dropped.txt"), header=True, index=False, sep=",")
		offtargets.to_csv(os.path.join(self.root, "temp", "offtargets.txt"), header=True, index=False, sep=",")


class CustomSQL_worker(QtCore.QRunnable):
	def __init__(self,database, sql_query):
		super().__init__()
		self.root = os.path.dirname(os.path.abspath(__name__))

		if not "temp" in os.listdir(self.root):
			os.mkdir(os.path.join(self.root, "temp"))

		self.database = database
		self.sql_query = sql_query
		

	@ QtCore.pyqtSlot()
	def run(self):
		sqlrunner = SQL(database=self.database)
		out = sqlrunner.custom_sql(statement=self.sql_query)
		out.to_csv(os.path.join(self.root, "temp", "query.txt"), header=True, index=False, sep=",")

		