import os
import pandas as pd
from Bio import SeqIO
import sqlite3

from sqlite3 import Error


def list_databases(self):
        return [databases for databases in os.listdir(os.path.join(self.root, "databases"))]

def _create_table(connection, sql_statement):
    try:
        c = connection.cursor()
        c.execute(sql_statement)
    except Error as e:
        print(e)


class Database:
    """
    Database tools, will be able to create and organism db and populate its global gRNA's
    """
    def __init__(self, database):
        self.root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.database = os.path.join(self.root, "databases", f"{database}.db")
        self.connection = sqlite3.connect(self.database)
        self.DB = database



    def create_new_database(self, gff_file, genome, multifasta):
        """
        requires a name for the database, the associated gff file path and a path to a fasta file with the all the genes
        no extention needed  for name
        """
        gff_headers = ["seqname", "source", "feature", "start", "stop", "score", "strand", "frame", "attribute"]

        gff = pd.read_csv(filepath_or_buffer=gff_file, sep="\t", names=gff_headers)

        genome_dict = {
        'header': [str(fasta.id).split("|")[0] for fasta in SeqIO.parse(open(genome), 'fasta')],
        'sequence': [str(fasta.seq) for fasta in SeqIO.parse(open(genome), 'fasta')],
        'length':	[len(fasta.seq) for fasta in SeqIO.parse(open(genome), 'fasta')]
        }

        mulitfasta_header = []
        def _gene_name_format(string):
            if "_" in string:
                return	string.replace("_", "")
            else:
                return string

        for elements in [str(fasta.id) for fasta in SeqIO.parse(open(multifasta), 'fasta')]:
                if "Locus_tag" in elements:
                    target_string = elements[elements.find("locus_tag="):]
                    target = _gene_name_format(string=target_string[:target_string.find("]")].split("=")[1])
                    mulitfasta_header.append(target)
                else:
                    target = _gene_name_format(string=elements.split("|")[0])
                    mulitfasta_header.append(target)


        multifasta_dict = {
        'header': mulitfasta_header,
        'sequence': [str(fasta.seq) for fasta in SeqIO.parse(open(multifasta), 'fasta')],
        'gene_length':	[len(fasta.seq) for fasta in SeqIO.parse(open(multifasta), 'fasta')]
        }

        genome_dataframe = pd.DataFrame(genome_dict)
        multifasta_dataframe = pd.DataFrame(multifasta_dict)
        gff.to_sql(name = "gff_file", con=self.connection, if_exists='replace', index = False)
        genome_dataframe.to_sql(name = "genome", con=self.connection, if_exists='replace', index = False)
        multifasta_dataframe.to_sql(name = "genes", con=self.connection, if_exists='replace', index = False)


    def create_gRNA_database(self, summary, offtarget, config_file):
        """
        Creates a global gRNA database that can be querried after creation using the predict gRNA tools.
        This will have a global gRNA table as well as its associated off targets and meta data. Requires pandas dataframe
        """
        summary_file = pd.read_csv(filepath_or_buffer=summary,
                index_col = False, sep = "\t")

        offtarget_file = pd.read_csv(filepath_or_buffer=offtarget,
                index_col = False, sep = "\t")

        config_file = pd.read_csv(filepath_or_buffer=config_file,
            sep="\t")

        summary_file = summary_file[['names', 'gRNAsPlusPAM', 'top5OfftargetTotalScore']]

        summary_file['PAM'] = ["".join(config_file['PAM_sequence']) for _ in range(len(summary_file.index))]
        summary_file['primer_length'] = [int(config_file['gRNA_size']) for _ in range(len(summary_file.index))]
        summary_file['mismatch_length'] = [int(config_file['max_mismatch_gRNA']) for _ in range(len(summary_file.index))]

        offtarget_file['PAM'] = ["".join(config_file['PAM_sequence']) for _ in range(len(offtarget_file.index))]
        offtarget_file['primer_length'] = [int(config_file['gRNA_size']) for _ in range(len(offtarget_file.index))]
        offtarget_file['mismatch_length'] = [int(config_file['max_mismatch_gRNA']) for _ in range(len(offtarget_file.index))]

        offtarget_file['gene'] = [gene.split("_")[0] for gene in offtarget_file['name']]

        offtarget_file.drop(columns=['inExon', 'inIntron', 'entrez_id', 'isCanonicalPAM', 'forViewInUCSC'], inplace=True)

        summary_file.to_sql(name = "global_gRNA", con=self.connection, if_exists='append', index = False)
        offtarget_file.to_sql(name = "global_offtarget", con=self.connection, if_exists='append', index = False)



class SQL:
    def __init__(self, database):
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.conn = sqlite3.connect(os.path.join(self.root, "databases", f"{database}.db"))
        self.pandasconn = sqlite3.connect(os.path.join(self.root, "databases", f"{database}.db"))
        self.conn.row_factory = lambda cursor, row: row[0]
        self.cursor = self.conn.cursor()

    def list_pams(self):
        try:
            return self.cursor.execute("SELECT DISTINCT PAM FROM global_gRNA").fetchall()
        except:
            return [""]

    def get_all_offtargets(self):
        try:
            return self.cursor.execute("SELECT DISTINCT annotation FROM global_offtarget").fetchall()
        except:
            return ""

    def list_mismatches(self, pam):
        try:
            out = self.cursor.execute("SELECT DISTINCT mismatch_length FROM global_gRNA WHERE PAM = '{pam}'".format(pam=pam)).fetchall()
            return [str(items) for items in out]
        except:
            return [""]

    def list_mismatches_processing(self):
        try:
            out = self.cursor.execute("SELECT DISTINCT mismatch_length FROM global_gRNA").fetchall()
            return [str(items) for items in out]
        except:
            return [""]

    def gene_count(self):
        try:
            return str(self.cursor.execute("SELECT COUNT(*) FROM genes").fetchone())
        except:
            return "No database detected"

    def gRNA_count(self):
        try:
            return str(self.cursor.execute("SELECT COUNT(*) FROM global_gRNA").fetchone())
        except:
            return "No database detected"

    def gRNA_count_selection(self, pam, mismatch):
        try:
            return str(self.cursor.execute("SELECT COUNT(*) FROM global_gRNA WHERE PAM = '{pam}' AND mismatch_length = '{mismatch}'".format(pam=pam,
                mismatch=mismatch)).fetchone())
        except:
            return "None detected"

    def genes_with_gRNA(self):
        try:
            gRNAs = pd.read_sql(sql="SELECT names FROM global_gRNA", con=self.pandasconn)
            gene_count = [genes.split("_")[0] for genes in gRNAs['names']]
            return str(len(set(gene_count)))
        except:
            return "None detected"

    def genes_with_gRNA_selection(self, pam, mismatch):
        try:
            gRNAs = pd.read_sql(sql="SELECT names FROM global_gRNA WHERE PAM = '{pam}' AND mismatch_length = '{mismatch}'".format(pam=pam,
                mismatch=mismatch),
                con=self.pandasconn)
            gene_count = [genes.split("_")[0] for genes in gRNAs['names']]
            return str(len(set(gene_count)))
        except:
            return "None detected"

    def gRNA_onReverse_selection(self, pam, mismatch):
        try:
            gRNAs = pd.read_sql(sql="SELECT names FROM global_gRNA WHERE PAM = '{pam}' AND mismatch_length = '{mismatch}'".format(pam=pam,
                mismatch=mismatch),
                con=self.pandasconn)
            out = [genes for genes in gRNAs['names'] if genes[-1] == 'r']
            return str(len(out))
        except:
            return "None detected"

    def get_global_gRNA(self, mismatch):
        try:
            return pd.read_sql(sql="SELECT * FROM global_gRNA WHERE mismatch_length = '{mismatch}'".format(mismatch=mismatch),
                con=self.pandasconn)
        except AssertionError as e:
            print(str(e))

    def get_gene_multifasta(self):
        try:
            return pd.read_sql(sql="SELECT header, sequence FROM genes",
                con=self.pandasconn)
        except:
            return "No files"

    def get_offtargets_by_mismatch(self, mismatch):
        try:
            out = pd.read_sql(sql="SELECT name, gRNAPlusPAM, OffTargetSequence, [n.mismatch], strand, chrom, chromStart, chromEnd, annotation, PAM, mismatch_length, gene FROM global_offtarget WHERE mismatch_length = '{mismatch}'".format(mismatch=mismatch),
                con=self.pandasconn)
            return out
        except:
            return "No database"

    def get_gene_sequence(self, gene):
        try:
            return str(self.cursor.execute("SELECT sequence FROM genes WHERE header = '{gene}'".format(gene=gene)).fetchone())
        except:
            return "None detected"

    def custom_sql(self, statement):
        out = pd.read_sql(sql="{statement}".format(statement=statement),
                con=self.pandasconn)
        return out
