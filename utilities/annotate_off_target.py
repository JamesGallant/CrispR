import sqlite3
import os

import pandas as pd

class Annotater:
	def __init__(self, config_file, offtarget_file):
		self.root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		self.config = config_file
		self.offtarget = offtarget_file

	def _pandas_parser(self, column_name):
		return [str(value) for value in self.config[column_name]][0]


	def _find_sql_database(self):
		availible_databases = [db.split(".")[0] for db in os.listdir(os.path.join(self.root, "databases"))]
		db_name = self._pandas_parser(column_name='organism')

		if " " in db_name:
			db_name = db_name.replace(" ", "_")

		if db_name not in availible_databases:
			raise ValueError(f"{db_name} database not found. availible_databases: {availible_databases}")

		return db_name


	def _database_query(self, position, query, chromosome):
		"""
		grabs the rows that corresponds to the start and stop positions per row in the self.offtarget dataframe
		"""
		if query.nunique()[0] > 1: # Chromosome in first column for this query
			target = query.loc[(int(position) >= query['start']) & (int(position) <= query['stop'])]
			target = target[target.iloc[:, 0] == chromosome].to_dict()
		else:
			target = query.loc[(int(position) >= query['start']) & (int(position) <= query['stop'])].to_dict()

		return target


	def annotate(self):
		"""
		Requires an input vector and a gff_file 
		"""

		anno = []

		current_org = self._find_sql_database()

		query_db = pd.read_sql(sql = "SELECT * FROM gff_file", con=sqlite3.connect(os.path.join(self.root, "databases", f"{current_org}.db")))
	

		for _, row in self.offtarget.iterrows():
			start = row['chromStart']
			stop = row['chromEnd']
			chromosome = row['chrom']

			get_genelist = self._database_query(position=start, query=query_db, chromosome=chromosome)
			
			if bool(get_genelist['attribute']):
				anno.append(str(list(get_genelist['attribute'].values())).split(";")[0].split("=")[1])
			else:
				anno.append("NA")
				

		self.offtarget['annotation'] = anno

		return self.offtarget