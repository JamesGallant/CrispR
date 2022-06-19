import os
import pandas as pd
from datetime import date
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
        self.root = os.path.dirname(os.path.abspath(__name__))
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
            'length': [len(fasta.seq) for fasta in SeqIO.parse(open(genome), 'fasta')]
        }

        mulitfasta_header = []

        def _gene_name_format(string):
            if "_" in string:
                return string.replace("_", "")
            else:
                return string

        for elements in [str(fasta.id) for fasta in SeqIO.parse(open(multifasta), 'fasta')]:
            if "Locus_tag" in elements:
                target_string = elements[elements.find("locus_tag="):]
                target = _gene_name_format(string=target_string[:target_string.find("]")].split("=")[1])
                target = target.lower()
                mulitfasta_header.append(target)
            else:
                target = _gene_name_format(string=elements.split("|")[0])
                target = target.lower()
                mulitfasta_header.append(target)

        multifasta_dict = {
            'header': mulitfasta_header,
            'sequence': [str(fasta.seq) for fasta in SeqIO.parse(open(multifasta), 'fasta')],
            'gene_length': [len(fasta.seq) for fasta in SeqIO.parse(open(multifasta), 'fasta')]
        }

        genome_dataframe = pd.DataFrame(genome_dict)
        multifasta_dataframe = pd.DataFrame(multifasta_dict)
        gff.to_sql(name="gff_file", con=self.connection, if_exists='replace', index=False)
        genome_dataframe.to_sql(name="genome", con=self.connection, if_exists='replace', index=False)
        multifasta_dataframe.to_sql(name="genes", con=self.connection, if_exists='replace', index=False)

    def create_gRNA_database(self, summary, offtarget, config_file):
        """
        Creates a global gRNA database that can be querried after creation using the predict gRNA tools.
        This will have a global gRNA table as well as its associated off targets and meta data. Requires pandas dataframe
        """
        summary_file = pd.read_csv(filepath_or_buffer=summary,
                                   index_col=False, sep="\t")

        offtarget_file = pd.read_csv(filepath_or_buffer=offtarget,
                                     index_col=False, sep="\t")

        config_file = pd.read_csv(filepath_or_buffer=config_file,
                                  sep="\t")

        summary_file = summary_file[['names', 'gRNAsPlusPAM', 'top5OfftargetTotalScore']]

        summary_file['PAM'] = ["".join(config_file['PAM_sequence']) for _ in range(len(summary_file.index))]
        summary_file['primer_length'] = [int(config_file['gRNA_size']) for _ in range(len(summary_file.index))]
        summary_file['mismatch_length'] = [int(config_file['max_mismatch_gRNA']) for _ in
                                           range(len(summary_file.index))]

        offtarget_file['PAM'] = ["".join(config_file['PAM_sequence']) for _ in range(len(offtarget_file.index))]
        offtarget_file['primer_length'] = [int(config_file['gRNA_size']) for _ in range(len(offtarget_file.index))]
        offtarget_file['mismatch_length'] = [int(config_file['max_mismatch_gRNA']) for _ in
                                             range(len(offtarget_file.index))]

        offtarget_file['gene'] = [gene.split("_")[0] for gene in offtarget_file['name']]

        offtarget_file.drop(columns=['inExon', 'inIntron', 'entrez_id', 'isCanonicalPAM', 'forViewInUCSC'],
                            inplace=True)

        summary_file.to_sql(name="global_gRNA", con=self.connection, if_exists='append', index=False)
        offtarget_file.to_sql(name="global_offtarget", con=self.connection, if_exists='append', index=False)

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

        tempfa = pd.read_sql("SELECT * FROM genome",
                             sqlite3.connect(os.path.join(os.getcwd(), "databases", f"{database_name}.db")))
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
            'Package': f'BSgenome.{common_name}.{self.department}.{strain}',
            'Title': f'Genome sequence of {self.org_name}',
            'Description': f'BSgenome genome obtained from the genome sequence of {self.org_name}',
            'Version': version,
            'Author': self.author,
            'License': "GPLv2",
            'organism': self.org_name,
            'common_name': common_name,
            'provider': 'NCBI',
            'genome': 'placeholder',
            'circ_seqs': 'character(0)',
            'source_url': 'www.example.com',
            'release_date': str(date.today()),
            'organism_biocview': f'{genus}_{species}',
            'BSgenomeObjname': common_name,
            'seqnames': seqnames,
            'seqs_srcdir': os.path.join(self.tempdir, "fastafiles")

        }

        with open(os.path.join(self.tempdir, "dcf_file.dcf"), 'w') as dcf:
            for k, v in dcf_file.items():
                dcf.write(f"{k}: {v} \n")

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
            out = self.cursor.execute(
                "SELECT DISTINCT mismatch_length FROM global_gRNA WHERE PAM = '{pam}'".format(pam=pam)).fetchall()
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
            return str(self.cursor.execute(
                "SELECT COUNT(*) FROM global_gRNA WHERE PAM = '{pam}' AND mismatch_length = '{mismatch}'".format(
                    pam=pam,
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
            gRNAs = pd.read_sql(
                sql="SELECT names FROM global_gRNA WHERE PAM = '{pam}' AND mismatch_length = '{mismatch}'".format(
                    pam=pam,
                    mismatch=mismatch),
                con=self.pandasconn)
            gene_count = [genes.split("_")[0] for genes in gRNAs['names']]
            return str(len(set(gene_count)))
        except:
            return "None detected"

    def gRNA_onReverse_selection(self, pam, mismatch):
        try:
            gRNAs = pd.read_sql(
                sql="SELECT names FROM global_gRNA WHERE PAM = '{pam}' AND mismatch_length = '{mismatch}'".format(
                    pam=pam,
                    mismatch=mismatch),
                con=self.pandasconn)
            out = [genes for genes in gRNAs['names'] if genes[-1] == 'r']
            return str(len(out))
        except:
            return "None detected"

    def get_global_gRNA(self, mismatch):
        try:
            return pd.read_sql(
                sql="SELECT * FROM global_gRNA WHERE mismatch_length = '{mismatch}'".format(mismatch=mismatch),
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
            out = pd.read_sql(
                sql="SELECT name, gRNAPlusPAM, OffTargetSequence, [n.mismatch], strand, chrom, chromStart, chromEnd, annotation, PAM, mismatch_length, gene FROM global_offtarget WHERE mismatch_length = '{mismatch}'".format(
                    mismatch=mismatch),
                con=self.pandasconn)
            return out
        except:
            return "No database"

    def get_gene_sequence(self, gene):
        try:
            return str(
                self.cursor.execute("SELECT sequence FROM genes WHERE header = '{gene}'".format(gene=gene)).fetchone())
        except:
            return "None detected"

    def custom_sql(self, statement):
        out = pd.read_sql(sql="{statement}".format(statement=statement),
                          con=self.pandasconn)
        return out

