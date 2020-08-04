import os
import pandas as pd
import sqlite3

from Bio import SeqIO
from .database_tools import Database

class Process:
	def __init__(self, gRNA_files_directory, config_file):
		self.root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		self.input_dir = gRNA_files_directory
		self.config =pd.read_csv(filepath_or_buffer=config_file, sep="\t")
		self.fasta_file = [str(filepath) for filepath in self.config['input_file']][0]
		self.db_tools = Database()
		
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


	def _config_file_parameters(self, column_name):
		"""
		The config file needs specific parsing, well handle it here given the column name.
		returns a string
		"""
		return [str(filepath) for filepath in self.config[column_name]][0]

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
			out['notes'].append(f"No gRNA detected based on PAM: {self._config_file_parameters(column_name='PAM_sequence')}")

		return out


	def gRNA_first_pass(self):
		"""
		first pass iteration to filter for possible matches. This function will return three dictionaries:
		canidates which have passed constraints, backup candidates that can be accepted in follow up rounds and those that need to be ignored.
		Contstraints are gRNA has to be on the reverse strand and the 20th position is either a G/A
		"""

		candidates = {
		'name':	[],
		'gRNA':	[],
		'notes': []
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
				if row['gRNAsPlusPAM'][19] == "A" or row['gRNAsPlusPAM'][19] == "G":
					candidates['name'].append(row['names'])
					candidates['gRNA'].append(row['gRNAsPlusPAM'])
					candidates['notes'].append("PASS")
				else:
					backup['name'].append(row['names'])
					backup['gRNA'].append(row['gRNAsPlusPAM'])
					backup['notes'].append(f"position 20 is {row['gRNAsPlusPAM'][19]}")
			else:
				dropped_gRNA['name'].append(row['names'])
				dropped_gRNA['gRNA'].append(row['gRNAsPlusPAM'])
				dropped_gRNA['notes'].append(f"The gRNA {row['names']} is not on the reverse strand")

		return candidates, backup, dropped_gRNA


	def correct_gRNA_dictionary(self, candidates, backup, dropped_gRNA):
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
					backup = pd.DataFrame(backup)
					backup_row_idx = _get_index(dataframe=backup, identifier=genes)
					candidates['name'].append(backup.iloc[backup_row_idx, 0])
					candidates['gRNA'].append(backup.iloc[backup_row_idx, 1])
					candidates['notes'].append(backup.iloc[backup_row_idx, 2])
					backup.drop(backup_row_idx, inplace=True)
					backup.reset_index(drop=True, inplace=True)
					backup = backup.to_dict()
					_check_altered_dict(input_dict=backup, identifier=genes)
				else:
					#Moves guide RNA's between backup and invalid think if this makes sense
					pass

		return candidates, backup, dropped_gRNA


	def has_offtarget(self, candidates, backup, dropped_gRNA):
		"""
		Takes the three input dictionaries, performs the function and returns them again. This function will check if there are
		off targets present, if not it will be moved to candidates. If off targets are present it will be moved to backup for
		follow up.
		"""
		candidates, backup, dropped_gRNA = self.correct_gRNA_dictionary(candidates=self.gRNA_first_pass()[0],
			backup=self.gRNA_first_pass()[1], dropped_gRNA=self.gRNA_first_pass()[2])

		for identifiers in candidates["name"]:
			cross_reference = self.summary.loc[self.summary['names'] == identifiers].to_dict()
			off_target_presence = str(list(cross_reference['top5OfftargetTotalScore'].values())[0])
			if off_target_presence != "nan":
				backup['name'].append(identifiers)
				backup['gRNA'].append(str(list(cross_reference['gRNAsPlusPAM'].values())[0]))
				backup['notes'].append("Has off target")
				candidates = pd.DataFrame(candidates)
				candidates = candidates[candidates['name'] != identifiers]
				candidates.reset_index(drop=True, inplace=True)
				candidates = candidates.to_dict()
				

		candidates['name'] = list(candidates['name'].values())
		candidates['gRNA'] = list(candidates['gRNA'].values())
		candidates['notes'] = list(candidates['notes'].values())

		candidates, backup, dropped_gRNA = self.correct_gRNA_dictionary(candidates=candidates, backup=backup, dropped_gRNA=dropped_gRNA)

		return candidates, backup, dropped_gRNA


	def annotate_off_target(self):
		"""
		Annotate the off target with the gff file. This will add a column on the self.offtarget variable data frame. Use this with threading.

		"""

		def _database_query(position, query, chromosome):
			"""
			grabs the rows that corresponds to the start and stop positions per row in the self.offtarget dataframe
			"""
			if query.nunique()[0] > 1: # Chromosome in first column for this query
				target = query.loc[(int(position) >= query['start']) & (int(position) <= query['stop'])]
				target = target[target.iloc[:, 0] == chromosome].to_dict()
			else:
				target = query.loc[(int(position) >= query['start']) & (int(position) <= query['stop'])].to_dict()

			return target


		start_anno = []
		stop_anno = []

		availible_databases = [db.split(".")[0] for db in os.listdir(os.path.join(self.root, "databases"))]
		current_org = self._config_file_parameters(column_name='organism')

		if " " in current_org:
			current_org = current_org.replace(" ", "_")

		if current_org not in availible_databases:
			raise ValueError(f"No database found for {current_org}. availible_databases: {self.db_tools.list_databases()}")

		query_db = pd.read_sql(sql = "SELECT * FROM gff_file", con=sqlite3.connect(os.path.join(self.root, "databases", f"{current_org}.db")))
	

		for _, row in self.offtarget.iterrows():
			start = row['chromStart']
			stop = row['chromEnd']
			chromosome = row['chrom']

			get_left_genelist = _database_query(position=start, query=query_db, chromosome=chromosome)
			get_right_genelist = _database_query(position=start, query=query_db, chromosome=chromosome)

			start_anno.append(str(list(get_left_genelist['attribute'].values())[0]).split(";")[0].split("=")[1])
			stop_anno.append(str(list(get_right_genelist['attribute'].values())[0]).split(";")[0].split("=")[1])


		self.offtarget['start_annotation'] = start_anno
		self.offtarget['stop_annotation'] = stop_anno

		return self.offtarget


	def create_gRNA_db(self, candidates, dropped_gRNA):
		"""
		Creates a sql database for the guide RNA's found for a specific PAM. In this case backup should not be neccesary so we 
		disregard it
		"""	

		no_detection = self.no_detection()

		#candidates, _, dropped_gRNA = self.has_offtarget(candidates=self.correct_gRNA_dictionary(candidates=self.gRNA_first_pass()[0],
		#	backup=self.gRNA_first_pass()[1],
		#	dropped_gRNA=self.gRNA_first_pass()[2])[0],
		#	backup=self.correct_gRNA_dictionary(candidates=self.gRNA_first_pass()[0],
		#	backup=self.gRNA_first_pass()[1],
		#	dropped_gRNA=self.gRNA_first_pass()[2])[1],
		#	dropped_gRNA=self.correct_gRNA_dictionary(candidates=self.gRNA_first_pass()[0],
		#	backup=self.gRNA_first_pass()[1],
		#	dropped_gRNA=self.gRNA_first_pass()[2])[2])


		PAM = self._config_file_parameters(column_name="PAM_sequence")
		DB_file = self._config_file_parameters(column_name="organism")

		if " " in DB_file:
			DB_file.replace(" ", "_")

		availible_databases = [db.split(".")[0] for db in os.listdir(os.path.join(self.root, "databases"))]
		if f"{DB_file}.db" not in availible_databases:
			raise ValueError("SQL database not found, make sure config file organisms parameter matches the database name")

		sql_dataframes = [candidates, dropped_gRNA, no_detection]

		switch = {
		0: 'valid_gRNA',
		1: 'invalid_gRNA',
		2: 'no_gRNA_detected'
		}

		for idx, df in enumerate(sql_dataframes):
			df = pd.DataFrame(df)
			gene_names = []
	
			for _, row in df.iterrows():
				gene_names.append(row['name'].split("_")[0])

			df['genes'] = gene_names
			df['PAM'] = PAM
			gene_names.clear()

			df.to_sql(name=switch.get(idx), 
				con=sqlite3.connect(os.path.join(self.root, "databases", f"{DB_file}.db")),
				if_exists="append", 
				index=False)


	def create_offtarget_db(self):
		# needs work
		"""
		In this function the canidate gRNA's with off target effects will be cross referenced with the off target genes. candidates will be moved between
		backup and valids iteratively if matched with a masked gene
		"""
		self.annotate_off_target()

		candidates, backup, dropped_gRNA = self.has_offtarget(candidates=self.correct_gRNA_dictionary(candidates=self.gRNA_first_pass()[0],
			backup=self.gRNA_first_pass()[1],
			dropped_gRNA=self.gRNA_first_pass()[2])[0],
			backup=self.correct_gRNA_dictionary(candidates=self.gRNA_first_pass()[0],
			backup=self.gRNA_first_pass()[1],
			dropped_gRNA=self.gRNA_first_pass()[2])[1],
			dropped_gRNA=self.correct_gRNA_dictionary(candidates=self.gRNA_first_pass()[0],
			backup=self.gRNA_first_pass()[1],
			dropped_gRNA=self.gRNA_first_pass()[2])[2])

		candidates = pd.DataFrame(backup)
		off_target_represented = candidates[candidates.notes == "Has off target"]
		off_target_represented.reset_index(drop=True, inplace=True)

		off_target_df ={
		'name': [],
		'gRNA': [],
		'off_target_sequence': [],
		'flanking_sequence':	[],
		'strand':	[],
		'mismatches': [],
		'chromosome':	[],
		'start':	[],
		'stop':	[],
		'gene':	[]
		}
		#handle off targets in list of genes to be masked
		#needs a function to process the the pandas dataframe for appending to a dictionary
		if len(off_target_represented.index) > 1:
			for _, row in off_target_represented.iterrows():
				gRNA_name = row['name']
				gene_name = gRNA_name.split("_")[0]
				off_target_info = pd.DataFrame(self.offtarget.query('`name` == @gRNA_name'))
				off_target_info.query('start_annotation != @gene_name', inplace=True)
				off_target_info.query('stop_annotation != @gene_name', inplace=True)
				off_target_info = off_target_info[off_target_info['strand'] == "-"]
				off_target_info.reset_index(drop=True, inplace=True)
				print(off_target_info)
				

		else:
			gRNA_name = off_target_represented['name'].to_string(index=False).split(" ")[-1]
			gene_name = gRNA_name.split("_")[0]
			off_target_info = pd.DataFrame(self.offtarget.query('name == @gRNA_name'))
			off_target_info.query('start_annotation != @gene_name', inplace=True)
			off_target_info.query('stop_annotation != @gene_name', inplace=True)
			off_target_info = off_target_info[off_target_info['strand'] == "-"]
			off_target_info.reset_index(drop=True, inplace=True)
			print(off_target_info['name'])
			