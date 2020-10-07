import os
import pandas as pd
import time
import json
import multiprocessing

from PyQt5 import QtWidgets, QtGui, QtCore

from src.py.database_tools import SQL
from src.py.pandashandler import PandasModel
from src.py.workers import CustomSQL_worker
from src.py.API_calls import TaxonomyAPI
from src.py.misc_functions import possible_pams_ranked


class ImportGenes(QtWidgets.QDialog):

    def __init__(self, *args, **kwargs):
        super(ImportGenes, self).__init__(*args, **kwargs)

        self.root = os.path.dirname(os.path.abspath(__name__))
        self.setWindowTitle("Import gene lists")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))

        exitbtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel

        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.accepted.connect(self.accept)
        self.exit.rejected.connect(self.reject)

        self.gene_mask_filepath = {
            'genes': "",
            'masks': ""
        }

        self.initiateDialogUI()

    def initiateDialogUI(self):

        self.targetGenes_textinput = QtWidgets.QTextEdit()

        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.import_controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def import_controls(self):
        groupbox = QtWidgets.QGroupBox("Import controls")

        self.genes_fileimport = QtWidgets.QPushButton('Upload genes', self)
        self.genes_fileimport.setToolTip("Upload list of genes, one per line")
        self.genes_fileimport.clicked.connect(lambda: self.file_import_function(dictkey="genes"))

        self.genes_textimport_label = QtWidgets.QLabel("Paste gene list")
        self.genes_textimport_label.setStyleSheet("font-weight: bold")
        self.genes_textimport = QtWidgets.QTextEdit()

        self.mask_fileimport = QtWidgets.QPushButton('Upload masks', self)
        self.mask_fileimport.setToolTip("Upload list of genes to mask, one per line")
        self.mask_fileimport.clicked.connect(lambda: self.file_import_function(dictkey="masks"))

        self.masks_textimport_label = QtWidgets.QLabel("Paste mask list")
        self.masks_textimport_label.setStyleSheet("font-weight: bold")
        self.masks_textimport = QtWidgets.QTextEdit()

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.genes_fileimport, 0, 0)
        grid.addWidget(self.mask_fileimport, 0, 1)
        grid.addWidget(self.genes_textimport_label, 1, 0)
        grid.addWidget(self.genes_textimport, 2, 0)
        grid.addWidget(self.masks_textimport_label, 1, 1)
        grid.addWidget(self.masks_textimport, 2, 1)

        groupbox.setLayout(grid)
        return groupbox

    def file_import_function(self, dictkey):
        """
		Handles the file upload dialog
		"""
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', '.')
        self.gene_mask_filepath.update({dictkey: filename[0]})

        if dictkey == "genes":
            if bool(self.gene_mask_filepath['genes']):
                self.genes_fileimport.setText("success")
                self.genes_fileimport.setStyleSheet("background: #A5D6A7; color: black")

                genefile = open(self.gene_mask_filepath['genes'], 'r')
                genefile_lines = genefile.readlines()
                for items in genefile_lines:
                    self.genes_textimport.append(items)

        else:
            if bool(self.gene_mask_filepath['masks']):
                self.mask_fileimport.setText("success")
                self.mask_fileimport.setStyleSheet("background: #A5D6A7; color: black")

                genefile = open(self.gene_mask_filepath['masks'])
                genefile_lines = genefile.readlines()
                for items in genefile_lines:
                    self.masks_textimport.append(items)

    def out(self):
        """
		Return important things here for proceesing in main window
		"""
        out = {
            'genes': [str(items) for items in self.genes_textimport.toPlainText().splitlines() if items != ""],
            'masks': [str(items) for items in self.masks_textimport.toPlainText().splitlines() if items != ""]
        }

        return out


class ExportDatabase(QtWidgets.QDialog):
    def __init__(self, *args, **kwargs):
        super(ExportDatabase, self).__init__(*args, **kwargs)
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.setWindowTitle("Export Database")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))
        exitbtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel

        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.accepted.connect(self.accept)
        self.exit.rejected.connect(self.reject)

        self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                    os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
        self.availible_databases.append("all")
        self.initiateDialogUI()

    def initiateDialogUI(self):
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.export_controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def export_controls(self):
        groupbox = QtWidgets.QGroupBox("Export controls")
        self.databases = QtWidgets.QComboBox()
        self.databases.setToolTip("Choose database to export")
        self.databases.addItems(self.availible_databases)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.databases, 0, 0)
        groupbox.setLayout(grid)
        return groupbox

    def files(self):
        # destination = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory"))
        return str(self.databases.currentText()).replace(" ", "_")


class ExportGuideRNA(QtWidgets.QDialog):
    def __init__(self, candidates, backup, dropped, offtargets):
        super(ExportGuideRNA, self).__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.setWindowTitle("Export guide RNA data")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))
        exitbtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel

        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.accepted.connect(self.accept)
        self.exit.rejected.connect(self.reject)

        self.candidates = pd.DataFrame(candidates)
        self.backup = pd.DataFrame(backup)
        self.dropped = pd.DataFrame(dropped)
        self.offtargets = pd.DataFrame(offtargets)

        self.initiateDialogUI()

    def initiateDialogUI(self):
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.export_controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def export_controls(self):
        groupbox = QtWidgets.QGroupBox("Export controls")
        self.dataframes = QtWidgets.QComboBox()

        self.export_to_excel = QtWidgets.QRadioButton("Excel")
        self.export_to_excel.setChecked(True)
        self.export_to_excel.toggled.connect(lambda: self.exporter_toggle(self.export_to_excel))
        self.export_to_txt = QtWidgets.QRadioButton("Text")
        self.export_to_txt.toggled.connect(lambda: self.exporter_toggle(self.export_to_txt))

        self.exporter_group = QtWidgets.QButtonGroup()
        self.exporter_group.addButton(self.export_to_excel)
        self.exporter_group.addButton(self.export_to_txt)
        self.exporter_group.setExclusive(True)

        self.delimeter_dict = {
            'space': " ",
            'tab': "\t",
            'comma': ",",
            'semi-colon': ';',
            'colon': ":"
        }

        self.delimeter_list = [delim for delim in self.delimeter_dict.keys()]
        self.delimeter = QtWidgets.QComboBox()
        self.delimeter.addItems(self.delimeter_list)
        self.delimeter.setEnabled(False)

        self.dataframes.setToolTip("Choose data to export")
        self.availible_dataframes = ["All", "Candidates", "Backup", "Dropped", "Offtargets"]
        self.dataframes.addItems(self.availible_dataframes)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.export_to_excel, 0, 0)
        grid.addWidget(self.export_to_txt, 0, 1)
        grid.addWidget(self.dataframes, 1, 0)
        grid.addWidget(self.delimeter, 1, 1)
        groupbox.setLayout(grid)
        return groupbox

    def exporter_toggle(self, radio):
        if radio.text() == "Excel":
            self.delimeter.setEnabled(False)
        else:
            self.delimeter.setEnabled(True)

    def accept(self):
        destination = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory"))
        files_to_export = str(self.dataframes.currentText())

        if self.export_to_txt.isChecked():
            if files_to_export == "All":
                self.candidates.to_csv(os.path.join(destination, "candidates.txt"),
                                       header=True,
                                       index=False,
                                       sep=self.delimeter_dict.get(self.delimeter.currentText()))

                self.backup.to_csv(os.path.join(destination, "Backup.txt"),
                                   header=True,
                                   index=False,
                                   sep=self.delimeter_dict.get(self.delimeter.currentText()))

                self.dropped.to_csv(os.path.join(destination, "Dropped.txt"),
                                    header=True,
                                    index=False,
                                    sep=self.delimeter_dict.get(self.delimeter.currentText()))

                self.offtargets.to_csv(os.path.join(destination, "Offtargets.txt"),
                                       header=True,
                                       index=False,
                                       sep=self.delimeter_dict.get(self.delimeter.currentText()))

            if files_to_export == "Candidates":
                self.candidates.to_csv(os.path.join(destination, "candidates.txt"),
                                       header=True,
                                       index=False,
                                       sep=self.delimeter_dict.get(self.delimeter.currentText()))

            if files_to_export == "Backup":
                self.backup.to_csv(os.path.join(destination, "Backup.txt"),
                                   header=True,
                                   index=False,
                                   sep=self.delimeter_dict.get(self.delimeter.currentText()))

            if files_to_export == "Dropped":
                self.dropped.to_csv(os.path.join(destination, "Dropped.txt"),
                                    header=True,
                                    index=False,
                                    sep=self.delimeter_dict.get(self.delimeter.currentText()))

            if files_to_export == "Offtargets":
                self.offtargets.to_csv(os.path.join(destination, "Offtargets.txt"),
                                       header=True,
                                       index=False,
                                       sep=self.delimeter_dict.get(self.delimeter.currentText()))

        if self.export_to_excel.isChecked():
            if files_to_export == "All":
                self.candidates.to_excel(os.path.join(destination, "Candidates.xlsx"),
                                         sheet_name="Sheet1",
                                         header=True,
                                         index=False)

                self.backup.to_excel(os.path.join(destination, "Backup.xlsx"),
                                     sheet_name="Sheet1",
                                     header=True,
                                     index=False)

                self.dropped.to_excel(os.path.join(destination, "Dropped.xlsx"),
                                      sheet_name="Sheet1",
                                      header=True,
                                      index=False)

                self.offtargets.to_excel(os.path.join(destination, "Offtarget.xlsx"),
                                         sheet_name="Sheet1",
                                         header=True,
                                         index=False)

            if files_to_export == "Candidates":
                self.candidates.to_excel(os.path.join(destination, "Candidates.xlsx"),
                                         sheet_name="Sheet1",
                                         header=True,
                                         index=False)

            if files_to_export == "Backup":
                self.backup.to_excel(os.path.join(destination, "Backup.xlsx"),
                                     sheet_name="Sheet1",
                                     header=True,
                                     index=False)

            if files_to_export == "Dropped":
                self.dropped.to_excel(os.path.join(destination, "Dropped.xlsx"),
                                      sheet_name="Sheet1",
                                      header=True,
                                      index=False)

            if files_to_export == "Offtargets":
                self.offtargets.to_excel(os.path.join(destination, "Offtarget.xlsx"),
                                         sheet_name="Sheet1",
                                         header=True,
                                         index=False)
        self.close()


class CreateDatabaseDialog(QtWidgets.QDialog):
    def __init__(self):
        super(CreateDatabaseDialog, self).__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))
        exitbtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        self.setWindowTitle("Create new database")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))
        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.accepted.connect(self.accept)
        self.exit.rejected.connect(self.reject)

        self.filepaths = {
            'genome': "",
            'multifasta': "",
            'gff': "",
        }

        self.initiateDialogUI()

    def initiateDialogUI(self):
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def controls(self):
        groupBox = QtWidgets.QGroupBox("Create new database")

        self.author_name = QtWidgets.QLabel("Author")
        self.author_input = QtWidgets.QLineEdit("your name")
        self.author_input.setToolTip("Pretty self explanatory I think")

        self.company_name = QtWidgets.QLabel("Department")
        self.company_input = QtWidgets.QLineEdit("vumc")
        self.company_input.setToolTip("The company or department you work for")

        self.organism_name = QtWidgets.QLabel("Organism name")
        self.organism_input = QtWidgets.QLineEdit("Mycobacterium tuberculosis H37Rv")
        self.organism_input.setToolTip("This must have the convention: Genus species strain")

        self.common_name = QtWidgets.QLabel("Common name")
        self.common_input = QtWidgets.QLineEdit("Mtb")
        self.common_input.setToolTip("What do you call it in day to day life")

        self.genome_fasta_uploadlabel = QtWidgets.QLabel("Genome fasta")
        self.genome_fasta_uploaddescription = QtWidgets.QLabel("No file loaded")
        self.genome_fasta_uploadbtn = QtWidgets.QPushButton('upload', self)
        self.genome_fasta_uploadbtn.setToolTip('upload a genome fasta file')
        self.genome_fasta_uploadbtn.clicked.connect(lambda: self._open_file(dictkey='genome'))

        self.gene_fasta_uploadlabel = QtWidgets.QLabel("Gene fasta")
        self.gene_fasta_uploaddescription = QtWidgets.QLabel("No file loaded")
        self.gene_fasta_uploadbtn = QtWidgets.QPushButton('upload', self)
        self.gene_fasta_uploadbtn.setToolTip('Upload a multifasta file containing the genes')
        self.gene_fasta_uploadbtn.clicked.connect(lambda: self._open_file(dictkey='multifasta'))

        self.gff_uploadlabel = QtWidgets.QLabel("GFF file")
        self.gff_uploaddescription = QtWidgets.QLabel("No file loaded")
        self.gff_uploadbtn = QtWidgets.QPushButton('upload', self)
        self.gff_uploadbtn.setToolTip('Upload a GFF file')
        self.gff_uploadbtn.clicked.connect(lambda: self._open_file(dictkey='gff'))

        grid = QtWidgets.QGridLayout()

        # magicNumbers: rowidx, colidx, rowspan, colspan
        grid.addWidget(self.author_name, 0, 0, )
        grid.addWidget(self.author_input, 0, 1)
        grid.addWidget(self.company_name, 0, 2)
        grid.addWidget(self.company_input, 0, 3)
        grid.addWidget(self.organism_name, 1, 0)
        grid.addWidget(self.organism_input, 1, 1)
        grid.addWidget(self.common_name, 1, 2)
        grid.addWidget(self.common_input, 1, 3)
        grid.addWidget(self.genome_fasta_uploadlabel, 2, 0)
        grid.addWidget(self.genome_fasta_uploadbtn, 2, 1)
        grid.addWidget(self.genome_fasta_uploaddescription, 2, 2, 1, 4)
        grid.addWidget(self.gene_fasta_uploadlabel, 3, 0)
        grid.addWidget(self.gene_fasta_uploadbtn, 3, 1)
        grid.addWidget(self.gene_fasta_uploaddescription, 3, 2, 1, 4)
        grid.addWidget(self.gff_uploadlabel, 4, 0)
        grid.addWidget(self.gff_uploadbtn, 4, 1)
        grid.addWidget(self.gff_uploaddescription, 4, 2, 1, 4)

        groupBox.setLayout(grid)

        return groupBox

    def _open_file(self, dictkey):
        """
		Changes the upload button colour based on incroporation of paths in a dictionary
		green is active
		"""
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', '.')
        self.filepaths.update({dictkey: filename[0]})

        if bool(self.filepaths["genome"]):
            genome_filepath = self.filepaths["genome"]
            genome_file, genome_extention = os.path.splitext(genome_filepath)
            genome_file = genome_file.split("/")[-1]

            self.genome_fasta_uploaddescription.setText(f"{genome_file} is loaded")

            if genome_extention == ".fa" or genome_extention == ".fasta":
                self.genome_fasta_uploadbtn.setText("success")
                self.genome_fasta_uploadbtn.setStyleSheet("background: #A5D6A7; color: black")
            else:
                self.genome_fasta_uploadbtn.setText("Failed")
                self.genome_fasta_uploadbtn.setStyleSheet("background: #FFAB91; color: black")

        if bool(self.filepaths["multifasta"]):
            gene_filepath = self.filepaths["multifasta"]
            gene_file, gene_extention = os.path.splitext(gene_filepath)
            gene_file = gene_file.split("/")[-1]

            self.gene_fasta_uploaddescription.setText(f"{gene_file} is loaded")

            if gene_extention == ".fa" or gene_extention == ".fasta":
                self.gene_fasta_uploadbtn.setText("success")
                self.gene_fasta_uploadbtn.setStyleSheet("background: #A5D6A7; color: black")
            else:
                self.gene_fasta_uploadbtn.setText("Failed")
                self.gene_fasta_uploadbtn.setStyleSheet("background: #FFAB91; color: black")

        if bool(self.filepaths["gff"]):
            gff_filepath = self.filepaths["gff"]
            gff_file, gff_extention = os.path.splitext(gff_filepath)
            gff_file = gff_file.split("/")[-1]

            self.gff_uploaddescription.setText(f"{gff_file} is loaded")

            if gff_extention == ".gff":
                self.gff_uploadbtn.setText("success")
                self.gff_uploadbtn.setStyleSheet("background: #A5D6A7; color: black")
            else:
                self.gff_uploadbtn.setText("Failed")
                self.gff_uploadbtn.setStyleSheet("background: #FFAB91; color: black")

    def filepaths_out(self):
        return self.filepaths

    def metadata_out(self):
        return {
            'author': self.author_input.text(),
            'department': self.company_input.text(),
            'org_name': self.organism_input.text(),
            'common_name': self.common_input.text()
        }


class DeleteDatabaseDialog(QtWidgets.QDialog):
    def __init__(self, database_list, bsgenome_list):
        super(DeleteDatabaseDialog, self).__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))
        exitbtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        self.setWindowTitle("Delete database")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))
        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.accepted.connect(self.accept)
        self.exit.rejected.connect(self.reject)

        self.databases = database_list
        self.bsgenome = bsgenome_list
        self.initiateDialogUI()

    def initiateDialogUI(self):
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def controls(self):
        groupBox = QtWidgets.QGroupBox("Delete database")

        label = QtWidgets.QLabel("<b>Deleting databases is permanent <\b>")
        self.database_box = QtWidgets.QComboBox()
        self.database_box.setToolTip("Choose database to delete")
        self.database_box.addItems(self.databases)
        self.database_box.currentIndexChanged.connect(self.database_box_change)

        self.bsgenome_box = QtWidgets.QComboBox()
        self.bsgenome_box.setToolTip("Choose bsgenome object to delete, must match the database")
        self.bsgenome_box.addItems(self.bsgenome)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(label, 0, 0, 1, 1)
        grid.addWidget(self.database_box, 1, 0)
        grid.addWidget(self.bsgenome_box, 1, 1)
        groupBox.setLayout(grid)

        return groupBox

    def database_box_change(self, value):
        self.bsgenome_box.setCurrentIndex(value)

    def out(self):
        return {
            'database': str(self.database_box.currentText()),
            'bsgenome': str(self.bsgenome_box.currentText())
        }


class SQLDialog(QtWidgets.QDialog):
    def __init__(self, database_list):
        super(SQLDialog, self).__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.setMinimumSize(630, 150)
        self.threadingPool = QtCore.QThreadPool()
        self.dataframe = pd.DataFrame()

        if not "temp" in os.listdir(self.root):
            os.mkdir(os.path.join(self.root, "temp"))

        self.setWindowTitle("Query SQL database")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))
        exitbtn = QtWidgets.QDialogButtonBox.Close
        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.rejected.connect(self.reject)

        self.database_list = database_list
        self.table_list = ['genome', 'genes', 'gff_file', 'global_gRNA', 'global_offtarget']
        self.initiateDialogUI()

    def initiateDialogUI(self):
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def controls(self):
        groupBox = QtWidgets.QGroupBox("SQL query editor")

        database_label = QtWidgets.QLabel("Database")
        self.database_combobox = QtWidgets.QComboBox()
        self.database_combobox.addItems(self.database_list)

        query_label = QtWidgets.QLabel("SQL statement")
        self.query_editor = QtWidgets.QLineEdit()

        self.table_combobox = QtWidgets.QComboBox()
        self.table_combobox.addItems(self.table_list)
        self.table_combobox.activated.connect(self.update_editor)
        self.table_combobox.setToolTip("This is only to show the tables for SQL statement")

        submit_btn = QtWidgets.QPushButton("Submit")
        submit_btn.clicked.connect(self.run_sql)

        self.progressbar = QtWidgets.QProgressBar(self)
        self.progressbar_val = 0
        self.progressbar.setValue(self.progressbar_val)
        self.progressbar.setMaximum(100)

        self.tableView = QtWidgets.QTableView()

        grid = QtWidgets.QGridLayout()
        grid.addWidget(database_label, 0, 0)
        grid.addWidget(query_label, 0, 1)
        grid.addWidget(self.database_combobox, 1, 0)
        grid.addWidget(self.query_editor, 1, 1)
        grid.addWidget(self.table_combobox, 2, 0)
        grid.addWidget(submit_btn, 2, 1)
        grid.addWidget(self.progressbar, 3, 0, 1, 0)
        grid.addWidget(self.tableView, 4, 0, 1, 0)
        groupBox.setLayout(grid)

        return groupBox

    def update_editor(self):
        self.query_editor.clear()
        sql_statement = f"SELECT * FROM {self.table_combobox.currentText()}"
        self.query_editor.setText(sql_statement)

    def run_sql(self):
        database_file = str(self.database_combobox.currentText()).replace(" ", "_")
        statement = self.query_editor.text()

        worker = CustomSQL_worker(database=database_file, sql_query=statement)
        self.threadingPool.start(worker)

        while self.threadingPool.activeThreadCount() == 1:
            QtWidgets.QApplication.processEvents()

            if self.progressbar_val < 90:
                self.progressbar_val += 1
                self.progressbar.setValue(self.progressbar_val)
                time.sleep(0.5)

        if self.threadingPool.waitForDone():
            self.dataframe = pd.read_csv(filepath_or_buffer=os.path.join(self.root, "temp", "query.txt"), sep=",")

            while self.progressbar_val < 100:
                self.progressbar_val += 1
                self.progressbar.setValue(self.progressbar_val)

        self.progressbar_val = 0
        self.progressbar.setValue(self.progressbar_val)
        model = PandasModel(self.dataframe)
        self.tableView.setModel(model)


class StatsDialog(QtWidgets.QDialog):
    def __init__(self, database_list):
        super(StatsDialog, self).__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))

        self.setWindowTitle("Guide RNA statistics")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))
        exitbtn = QtWidgets.QDialogButtonBox.Close
        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.rejected.connect(self.reject)

        self.availible_databases = database_list

        try:
            self.sql_runner = SQL(database=self.availible_databases[0].replace(" ", "_"))
            self.availible_pams = self.sql_runner.list_pams()
            self.metadata_mismatches_value = self.sql_runner.list_mismatches(pam=self.availible_pams[0])
            self.processeed_mismatches_value = self.sql_runner.list_mismatches_processing()

            self.metadata1_value = self.sql_runner.gene_count()
            self.metadata2_value = self.sql_runner.gRNA_count()
            self.metadata3_value = self.sql_runner.genes_with_gRNA()
            self.metadata4_value = f"{str(round(int(int(self.metadata3_value) / int(self.metadata1_value) * 100), 2))} %"

            self.metadata5_value = self.sql_runner.gRNA_count_selection(pam=self.availible_pams[0],
                                                                        mismatch=self.metadata_mismatches_value[0])

            self.metadata6_value = self.sql_runner.gRNA_onReverse_selection(pam=self.availible_pams[0],
                                                                            mismatch=self.metadata_mismatches_value[0])

            self.metadata7_value = self.sql_runner.genes_with_gRNA_selection(pam=self.availible_pams[0],
                                                                             mismatch=self.metadata_mismatches_value[0])

            self.metadata8_value = f"{str(round(int(int(self.metadata7_value) / int(self.metadata1_value) * 100), 2))} %"

        except:
            self.metadata_mismatches_value = []
            self.availible_pams = []
            self.metadata1_value = "No database"
            self.metadata2_value = "No database"
            self.metadata3_value = "No database"
            self.metadata4_value = "No database"

            self.metadata5_value = "No database"
            self.metadata6_value = "No database"
            self.metadata7_value = "No database"

            self.metadata8_value = "No database"
            self.processeed_mismatches_value = [""]

        self.initiateDialogUI()

    def initiateDialogUI(self):
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def controls(self):
        groupBox = QtWidgets.QGroupBox("Guide RNA statistics")

        self.current_databases_label = QtWidgets.QLabel("organisms")
        self.current_databases = QtWidgets.QComboBox()
        self.current_databases.addItems(self.availible_databases)
        self.current_databases.setToolTip("Choose your organism")
        self.current_databases.activated.connect(self.update_by_current_database)

        self.current_pams_label = QtWidgets.QLabel("pams")
        self.current_pams = QtWidgets.QComboBox()
        self.current_pams.addItems(self.availible_pams)
        self.current_pams.activated.connect(self.update_by_current_pams)

        self.metadata_mismatches_label = QtWidgets.QLabel("mismatches")
        self.metadata_mismatches = QtWidgets.QComboBox()
        self.metadata_mismatches.setToolTip("Choose mismatch length")
        self.metadata_mismatches.addItems(self.metadata_mismatches_value)
        self.metadata_mismatches.activated.connect(self.update_by_current_mismatch)

        self.metadata_global_header = QtWidgets.QLabel("global stats")
        self.metadata_global_header.setAlignment(QtCore.Qt.AlignCenter)
        self.metadata_selected_header = QtWidgets.QLabel("stats based on selection")
        self.metadata_selected_header.setAlignment(QtCore.Qt.AlignCenter)

        self.metadata_label1 = QtWidgets.QLabel("Total genes:")
        self.metadata_label1.setStyleSheet("font-weight: bold")
        self.metadata_label1.setToolTip("Total number of genes present")

        self.metadata1 = QtWidgets.QLabel(self.metadata1_value)
        self.metadata1.setFrameShape(QtWidgets.QFrame.Panel)
        self.metadata1.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.metadata1.setLineWidth(3)

        self.metadata_label2 = QtWidgets.QLabel("gRNA count:")
        self.metadata_label2.setStyleSheet("font-weight: bold")
        self.metadata_label2.setToolTip("Total number of guide RNA's present")

        self.metadata2 = QtWidgets.QLabel(self.metadata2_value)
        self.metadata2.setFrameShape(QtWidgets.QFrame.Panel)
        self.metadata2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.metadata2.setLineWidth(3)

        self.metadata_label3 = QtWidgets.QLabel("Mapped genes:")
        self.metadata_label3.setStyleSheet("font-weight: bold")
        self.metadata_label3.setToolTip(f"Genes with guide RNAs based on all metrics")

        self.metadata3 = QtWidgets.QLabel(self.metadata3_value)
        self.metadata3.setFrameShape(QtWidgets.QFrame.Panel)
        self.metadata3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.metadata3.setLineWidth(3)

        self.metadata_label4 = QtWidgets.QLabel("Genome coverage:")
        self.metadata_label4.setStyleSheet("font-weight: bold")
        self.metadata_label4.setToolTip(f"Percentage of the genome covered based on all metrics")

        self.metadata4 = QtWidgets.QLabel(self.metadata4_value)
        self.metadata4.setFrameShape(QtWidgets.QFrame.Panel)
        self.metadata4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.metadata4.setLineWidth(3)

        self.metadata_label5 = QtWidgets.QLabel("gRNA count:")
        self.metadata_label5.setStyleSheet("font-weight: bold")
        self.metadata_label5.setToolTip(f"Number of gRNA's based on selection")

        self.metadata5 = QtWidgets.QLabel(self.metadata5_value)
        self.metadata5.setFrameShape(QtWidgets.QFrame.Panel)
        self.metadata5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.metadata5.setLineWidth(3)

        self.metadata_label6 = QtWidgets.QLabel("gRNA's on reverse:")
        self.metadata_label6.setStyleSheet("font-weight: bold")
        self.metadata_label6.setToolTip(f"Number of gRNA's on the reverse strand based on selection")

        self.metadata6 = QtWidgets.QLabel(self.metadata6_value)
        self.metadata6.setFrameShape(QtWidgets.QFrame.Panel)
        self.metadata6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.metadata6.setLineWidth(3)

        self.metadata_label7 = QtWidgets.QLabel("Mapped genes:")
        self.metadata_label7.setStyleSheet("font-weight: bold")
        self.metadata_label7.setToolTip(f"Number of genes with PAMS based on selection")
        self.metadata7 = QtWidgets.QLabel(self.metadata7_value)
        self.metadata7.setFrameShape(QtWidgets.QFrame.Panel)
        self.metadata7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.metadata7.setLineWidth(3)

        self.metadata_label8 = QtWidgets.QLabel("Genome coverage:")
        self.metadata_label8.setStyleSheet("font-weight: bold")
        self.metadata_label8.setToolTip(f"Percentage of genomes covered based on selection")

        self.metadata8 = QtWidgets.QLabel(self.metadata8_value)
        self.metadata8.setFrameShape(QtWidgets.QFrame.Panel)
        self.metadata8.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.metadata8.setLineWidth(3)

        if not self.availible_pams:
            current_pams_message = "No guide RNA's have been predicted yet"
        else:
            current_pams_message = "Choose the PAM"

        self.current_pams.setToolTip(current_pams_message)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.current_databases_label, 0, 0)
        grid.addWidget(self.current_pams_label, 0, 1)
        grid.addWidget(self.metadata_mismatches_label, 0, 3)
        grid.addWidget(self.current_databases, 1, 0)
        grid.addWidget(self.current_pams, 1, 1, 1, 2)
        grid.addWidget(self.metadata_mismatches, 1, 3)
        grid.addWidget(self.metadata_global_header, 2, 0)
        grid.addWidget(self.metadata_selected_header, 2, 2)
        grid.addWidget(self.metadata_label1, 3, 0)
        grid.addWidget(self.metadata1, 3, 1)
        grid.addWidget(self.metadata_label2, 4, 0)
        grid.addWidget(self.metadata2, 4, 1)
        grid.addWidget(self.metadata_label3, 5, 0)
        grid.addWidget(self.metadata3, 5, 1)
        grid.addWidget(self.metadata_label4, 6, 0)
        grid.addWidget(self.metadata4, 6, 1)
        grid.addWidget(self.metadata_selected_header, 2, 2)
        grid.addWidget(self.metadata_label5, 3, 2)
        grid.addWidget(self.metadata5, 3, 3)
        grid.addWidget(self.metadata_label6, 4, 2)
        grid.addWidget(self.metadata6, 4, 3)
        grid.addWidget(self.metadata_label7, 5, 2)
        grid.addWidget(self.metadata7, 5, 3)
        grid.addWidget(self.metadata_label8, 6, 2)
        grid.addWidget(self.metadata8, 6, 3)
        groupBox.setLayout(grid)

        return groupBox

    def _metadata_helper(self, database, user_pam, user_mismatch):
        """
        Helper function to udpate metadata
        :return: Void
        """
        # widget clear
        self.current_pams.clear()
        self.metadata_mismatches.clear()
        self.metadata1.clear()
        self.metadata2.clear()
        self.metadata3.clear()
        self.metadata4.clear()
        self.metadata5.clear()
        self.metadata6.clear()
        self.metadata7.clear()
        self.metadata8.clear()

        # list clear
        self.availible_pams.clear()
        self.metadata_mismatches_value.clear()

        # calculate
        self.availible_pams = database.list_pams()
        if user_pam:
            self.metadata_mismatches_value = database.list_mismatches(pam=user_pam)
        else:
            self.metadata_mismatches_value = database.list_mismatches(pam=self.availible_pams[0])

        self.metadata1_value = database.gene_count()
        self.metadata2_value = database.gRNA_count()
        self.metadata3_value = database.genes_with_gRNA()

        self.metadata5_value = self.sql_runner.gRNA_count_selection(pam=user_pam,
                                                                    mismatch=user_mismatch)

        self.metadata6_value = self.sql_runner.gRNA_onReverse_selection(pam=user_pam,
                                                                        mismatch=user_mismatch)

        self.metadata7_value = self.sql_runner.genes_with_gRNA_selection(pam=user_pam,
                                                                         mismatch=user_mismatch)

        if self.availible_pams[0] == "":
            self.metadata4_value = "None detected"
            self.metadata8_value = "None detected"
        else:
            self.metadata4_value = f"{str(round(int(int(self.metadata3_value) / int(self.metadata1_value) * 100), 2))} %"
            self.metadata8_value = f"{str(round(int(int(self.metadata7_value) / int(self.metadata1_value) * 100), 2))} %"

        # populate widgets
        self.current_pams.addItems(self.availible_pams)
        self.metadata_mismatches.addItems(self.metadata_mismatches_value)

        pams_idx = self.current_pams.findText(user_pam, QtCore.Qt.MatchFixedString)
        mismatch_idx = self.metadata_mismatches.findText(str(user_mismatch), QtCore.Qt.MatchFixedString)

        if pams_idx == -1:
            self.current_pams.setCurrentIndex(0)
        else:
            self.current_pams.setCurrentIndex(pams_idx)

        if mismatch_idx == -1:
            self.metadata_mismatches.setCurrentIndex(0)
        else:
            self.metadata_mismatches.setCurrentIndex(mismatch_idx)

        self.metadata1.setText(self.metadata1_value)
        self.metadata2.setText(self.metadata2_value)
        self.metadata3.setText(self.metadata3_value)
        self.metadata4.setText(self.metadata4_value)
        self.metadata5.setText(self.metadata5_value)
        self.metadata6.setText(self.metadata6_value)
        self.metadata7.setText(self.metadata7_value)
        self.metadata8.setText(self.metadata8_value)

    def update_by_current_database(self, value):
        """
        updates the metadata when user chooses new organism
        """
        pams_value = str(self.current_pams.currentText())
        mismatch_value = str(self.metadata_mismatches.currentText())

        self.sql_runner = SQL(database=str(self.current_databases.itemText(value)).replace(" ", "_"))

        self._metadata_helper(database=self.sql_runner, user_pam=pams_value,
                              user_mismatch=mismatch_value)

    def update_by_current_pams(self, value):
        """
        Update metadata based on pams
        """
        organism = str(self.current_databases.currentText())
        pams_value = self.availible_pams[value]
        mismatch_value = str(self.metadata_mismatches.currentText())

        self.sql_runner = SQL(database=organism.replace(" ", "_"))

        self._metadata_helper(database=self.sql_runner, user_pam=pams_value,
                              user_mismatch=mismatch_value)

    def update_by_current_mismatch(self, value):
        """
        Update metadata based on mismatches
        """
        organism = str(self.current_databases.currentText())
        pams_value = str(self.current_pams.currentText())
        mismatch_value = self.metadata_mismatches_value[value]

        self.sql_runner = SQL(database=organism.replace(" ", "_"))

        self._metadata_helper(database=self.sql_runner, user_pam=pams_value,
                              user_mismatch=mismatch_value)


class SearchGrnaDialog(QtWidgets.QDialog):
    def __init__(self, database_list, cas9_list):
        super(SearchGrnaDialog, self).__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.availible_cas9_orgs = cas9_list
        self.pams_by_org_dict = possible_pams_ranked(cas9=self.availible_cas9_orgs[0])
        self.pams_by_org = list(self.pams_by_org_dict.keys())
        self.availible_databases = database_list

        self.setWindowTitle("Search new guide RNA's")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))

        exitbtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.accepted.connect(self.accept)
        self.exit.rejected.connect(self.reject)

        self.initiateDialogUI()

    def initiateDialogUI(self):
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def controls(self):
        groupBox = QtWidgets.QGroupBox("Search guide RNA's")

        self.cas9_origin_label = QtWidgets.QLabel("Cas9 origin")
        self.cas9_origin = QtWidgets.QComboBox()
        self.cas9_origin.addItems(self.availible_cas9_orgs)
        self.cas9_origin.setToolTip("The organism from which the cas9 allele was cloned")
        self.cas9_origin.activated.connect(self.update_pams)

        self.PAM_label = QtWidgets.QLabel("PAM")
        self.PAM_input = QtWidgets.QComboBox()
        self.PAM_input.addItems(self.pams_by_org)
        self.PAM_input.setToolTip("Choose the PAM")

        self.mismatch_label = QtWidgets.QLabel("Mismatches")
        self.mismatch = QtWidgets.QSpinBox()
        self.mismatch.setMinimum(0)
        self.mismatch.setMaximum(15)
        self.mismatch.setValue(4)
        self.mismatch.setToolTip(
            "Number of base mismathces in the off targets, use 0 to find guide RNA's with no off targets")

        self.cores_label = QtWidgets.QLabel("Cores")
        self.cores = QtWidgets.QSpinBox()
        self.cores.setMinimum(1)
        self.cores.setMaximum(multiprocessing.cpu_count())
        self.cores.setValue(int(multiprocessing.cpu_count()/2))
        self.setToolTip(f"Number of cores to use for the calculations. Max is {multiprocessing.cpu_count()}")

        self.circ_chromosome_label = QtWidgets.QLabel("Circular")
        self.circ_chromosome = QtWidgets.QLineEdit("TRUE")
        self.circ_chromosome.setToolTip("Is the crhomosome cirular, <b>TRUE</b> for bacteria")

        self.search_databases_label = QtWidgets.QLabel("Organism")
        self.search_databases = QtWidgets.QComboBox()
        self.search_databases.addItems(self.availible_databases)
        self.search_databases.setToolTip("Choose your organism")

        self.taxononmy_id_input_label = QtWidgets.QLabel("Taxonomy ID")
        self.taxononmy_id_input = QtWidgets.QLineEdit()
        self.taxononmy_id_input.setToolTip("Taxonomy id for your organism")

        self.lookup_tax_id_label = QtWidgets.QLabel("Lookup taxonomy ID")
        self.lookup_btn = QtWidgets.QPushButton("Lookup")
        self.lookup_btn.setToolTip("Find id corresponding to organism, requires internet connection. Its not always possible")
        self.lookup_btn.clicked.connect(self.lookup_tax_id_function)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.cas9_origin_label, 0, 0)
        grid.addWidget(self.cas9_origin, 0, 1)
        grid.addWidget(self.PAM_label, 1, 0)
        grid.addWidget(self.PAM_input, 1, 1)
        grid.addWidget(self.mismatch_label, 2, 0)
        grid.addWidget(self.mismatch, 2, 1)
        grid.addWidget(self.cores_label, 3, 0)
        grid.addWidget(self.cores, 3, 1)
        grid.addWidget(self.circ_chromosome_label, 4, 0)
        grid.addWidget(self.circ_chromosome, 4, 1)
        grid.addWidget(self.search_databases_label, 5, 0)
        grid.addWidget(self.search_databases, 5, 1)
        grid.addWidget(self.taxononmy_id_input_label, 6, 0)
        grid.addWidget(self.taxononmy_id_input, 6, 1)
        grid.addWidget(self.lookup_tax_id_label, 7, 0)
        grid.addWidget(self.lookup_btn, 7, 1)
        groupBox.setLayout(grid)
        return groupBox

    def lookup_tax_id_function(self):
        genus, species, strain = str(self.search_databases.currentText()).split(" ")
        tax_lookup = TaxonomyAPI(genus=genus, species=species)
        tax_id_raw = tax_lookup.get_taxon()
        tax_id_decoded = json.loads(tax_id_raw.replace("\'", "\""))
        tax_id = tax_id_decoded['taxId']
        self.taxononmy_id_input.setText(tax_id)

    def update_pams(self):
        print("change")
        self.pams_by_org_dict = possible_pams_ranked(cas9=self.cas9_origin.currentText())
        self.pams_by_org = list(self.pams_by_org_dict.keys())
        self.PAM_input.clear()
        self.PAM_input.addItems(self.pams_by_org)

    def out(self):
        return {
            'pam': self.PAM_input.currentText(),
            'mismatch': self.mismatch.value(),
            'chromosome': self.circ_chromosome.text(),
            'organism': self.search_databases.currentText(),
            'taxonomy_id': self.taxononmy_id_input.text(),
            'cores': self.cores.value()
        }


class DisplayGuideRNA(QtWidgets.QDialog):
    def __init__(self, database_list, cas9_list):
        super(DisplayGuideRNA, self).__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.setMinimumSize(630, 50)
        self.availible_databases = database_list
        self.availible_cas9 = cas9_list
        sqlrunner = SQL(database=str(self.availible_databases[0]).replace(" ", "_"))
        self.processeed_mismatches_value = sqlrunner.list_mismatches_processing()

        self.gene_mask_dict = {
            'genes': [],
            'masks': []
        }

        self.gene_mask_filepath = {
            'genes': "",
            'masks': ""
        }

        self.setWindowTitle("Display guide RNA's")
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))

        exitbtn = QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        self.exit = QtWidgets.QDialogButtonBox(exitbtn)
        self.exit.accepted.connect(self.accept)
        self.exit.rejected.connect(self.reject)

        self.initiateDialogUI()

    def initiateDialogUI(self):
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.controls())
        self.layout.addWidget(self.exit)
        self.setLayout(self.layout)

    def controls(self):
        groupBox = QtWidgets.QGroupBox("Search guide RNA's")

        self.genes_fileimport = QtWidgets.QPushButton('Upload genes', self)
        self.genes_fileimport.setToolTip("Upload list of genes, one per line")
        self.genes_fileimport.clicked.connect(lambda: self.file_import_function(dictkey="genes"))

        self.genes_textimport_label = QtWidgets.QLabel("Paste gene list")
        self.genes_textimport_label.setStyleSheet("font-weight: bold")
        self.genes_textimport = QtWidgets.QTextEdit()

        self.mask_fileimport = QtWidgets.QPushButton('Upload masks', self)
        self.mask_fileimport.setToolTip("Upload list of genes to mask, one per line")
        self.mask_fileimport.clicked.connect(lambda: self.file_import_function(dictkey="masks"))

        self.masks_textimport_label = QtWidgets.QLabel("Paste mask list")
        self.masks_textimport_label.setStyleSheet("font-weight: bold")
        self.masks_textimport = QtWidgets.QTextEdit()

        self.processed_primerlen_label = QtWidgets.QLabel("Primer length")
        self.processed_primerlen = QtWidgets.QSpinBox()
        self.processed_primerlen.setMinimum(0)
        self.processed_primerlen.setValue(5)
        self.processed_primerlen.setToolTip(
            "Increments the guide RNA lenght to find a a/g if not present, this is the max bases to increment")

        self.processed_addnucleotides_5_label = QtWidgets.QLabel("Add bases 5'")
        self.processed_addnucleotides_3_label = QtWidgets.QLabel("Add bases 3'")
        self.processed_addnucleotides_5 = QtWidgets.QLineEdit()
        self.processed_addnucleotides_3 = QtWidgets.QLineEdit()
        self.processed_addnucleotides_5.setToolTip("Add bases on the 5' of the primer: eg. NNN-primer")
        self.processed_addnucleotides_3.setToolTip("Add bases on the 3' of the primer: eg. primer-NNN")

        self.processed_organism_label = QtWidgets.QLabel("Organism")
        self.processed_organism = QtWidgets.QComboBox()
        self.processed_organism.addItems(self.availible_databases)
        self.processed_organism.setToolTip("Choose organism")
        self.processed_organism.activated.connect(self.update_processed_mismatches)

        self.processed_maxMismatch_label = QtWidgets.QLabel("Max mismatches")
        self.processed_maxMismatch = QtWidgets.QComboBox()
        self.processed_maxMismatch.addItems(self.processeed_mismatches_value)
        self.processed_maxMismatch.setToolTip("Maximum mismatches to consider")

        self.processed_maxPrimers_label = QtWidgets.QLabel("Max gRNA's")
        self.processed_maxPrimers = QtWidgets.QSpinBox()
        self.processed_maxPrimers.setMinimum(1)
        self.processed_maxPrimers.setValue(3)
        self.processed_maxPrimers.setToolTip("Maximum guide RNA's to add to database")

        self.cas9_label = QtWidgets.QLabel("Cas9 origin")
        self.cas9 = QtWidgets.QComboBox()
        self.cas9.addItems(self.availible_cas9)
        self.cas9.setToolTip("The organism from whence the cas9 came")

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.genes_fileimport, 0, 0)
        grid.addWidget(self.mask_fileimport, 0, 1)
        grid.addWidget(self.genes_textimport_label, 1, 0)
        grid.addWidget(self.masks_textimport_label, 1, 1)
        grid.addWidget(self.genes_textimport, 2, 0)
        grid.addWidget(self.masks_textimport, 2, 1)
        grid.addWidget(self.processed_addnucleotides_5_label, 3, 0)
        grid.addWidget(self.processed_addnucleotides_3_label, 3, 1)
        grid.addWidget(self.processed_addnucleotides_5, 4, 0)
        grid.addWidget(self.processed_addnucleotides_3, 4, 1)
        grid.addWidget(self.processed_organism_label, 5, 0)
        grid.addWidget(self.processed_maxMismatch_label, 5, 1)
        grid.addWidget(self.processed_organism, 6, 0)
        grid.addWidget(self.processed_maxMismatch, 6, 1)
        grid.addWidget(self.processed_maxPrimers_label, 7, 0)
        grid.addWidget(self.processed_primerlen_label, 7, 1)
        grid.addWidget(self.processed_maxPrimers, 8, 0)
        grid.addWidget(self.processed_primerlen, 8, 1)
        grid.addWidget(self.cas9_label, 9, 0)
        grid.addWidget(self.cas9, 10, 0)

        groupBox.setLayout(grid)
        return groupBox

    def file_import_function(self, dictkey):
        """
        Handles the file upload dialog
        """
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', '.')
        self.gene_mask_filepath.update({dictkey: filename[0]})

        if dictkey == "genes":
            if bool(self.gene_mask_filepath['genes']):
                self.genes_fileimport.setText("success")
                self.genes_fileimport.setStyleSheet("background: #A5D6A7; color: black")

                genefile = open(self.gene_mask_filepath['genes'], 'r')
                genefile_lines = genefile.readlines()
                for items in genefile_lines:
                    self.genes_textimport.append(items)

        else:
            if bool(self.gene_mask_filepath['masks']):
                self.mask_fileimport.setText("success")
                self.mask_fileimport.setStyleSheet("background: #A5D6A7; color: black")

                genefile = open(self.gene_mask_filepath['masks'])
                genefile_lines = genefile.readlines()
                for items in genefile_lines:
                    self.masks_textimport.append(items)

    def update_processed_mismatches(self, value):
        self.processed_maxMismatch.clear()
        self.processeed_mismatches_value.clear()
        sqlrunner = SQL(database=str(self.processed_organism.itemText(value).replace(" ", "_")))
        self.processeed_mismatches_value = sqlrunner.list_mismatches_processing()
        self.processed_maxMismatch.addItems(self.processeed_mismatches_value)

    def out(self):
        return {
            'genes': [str(items) for items in self.genes_textimport.toPlainText().splitlines() if items != ""],
            'masks': [str(items) for items in self.masks_textimport.toPlainText().splitlines() if items != ""],
            'max_primer_len': self.processed_primerlen.value(),
            'max_grna_count': self.processed_maxPrimers.value(),
            'max_mismatch': self.processed_maxMismatch.currentText(),
            'organism': self.processed_organism.currentText(),
            'nucleotides_5': self.processed_addnucleotides_5.text(),
            'nucleotides_3': self.processed_addnucleotides_3.text(),
            'cas9': self.cas9.currentText()
        }
