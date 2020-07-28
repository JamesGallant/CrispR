import os
import sqlite3
from sqlite3 import Error
import pandas as pd


def _create_connection(db_file):
	conn = None
	try:
		conn = sqlite3.connect(db_file)
		return conn
	except Error as e:
		print(e)

	return conn

def _create_table(connection, sql_statement):
    try:
        c = connection.cursor()
        c.execute(sql_statement)
    except Error as e:
        print(e)


class Database:
	def __init__(self):
		self.root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		self.gff_headers = ["seqname", "source", "feature", "start", "stop", "score", "strand", "frame", "attribute"]


	def create_new_database(self, name, gff_file):
		"""
		requires a name for the database and the associated gff file
		no extention needed 
		"""
		connection = _create_connection(db_file=os.path.join(self.root, "databases", f"{name}.db"))		

		gff = pd.read_csv(filepath_or_buffer=gff_file, sep="\t", names=self.gff_headers)

		gff.to_sql(name = "gff_file", con=connection, if_exists='replace', index = False)


	def list_databases(self):
		return [databases for databases in os.listdir(os.path.join(self.root, "databases"))]
		

