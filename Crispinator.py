import sys
import os
import qdarkstyle
import shutil
import sqlite3
import webbrowser

from src.py.workers import *
from src.py.dialogs import *
from src.py.pandashandler import PandasModel
from src.py.database_tools import Database, SQL, CreateDatabaseMethod

# from PyQt5 import QtWidgets, QtGui, QtCore

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        self.root = os.path.dirname(os.path.abspath(__file__))
        super().__init__()
        self.threadingPool = QtCore.QThreadPool()
        self.candidate_gRNA_df = None
        self.backup_gRNA_df = None
        self.dropped_gRNA_df = None
        self.offtarget_df = None
        self.database_querried = False

        if "temp" in os.listdir(self.root):
            shutil.rmtree(os.path.join(self.root, "temp"))

        try:
            self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                        os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
            self.availible_bsgenome = [str(os.path.splitext(items)[0]) for items in
                                       os.listdir(os.path.join(self.root, "src", "local_packages"))]

            self.availible_cas9 = ["Streptococcus thermophilus"]
        except:
            self.availible_databases = [""]
            self.availible_bsgenome = [""]
            self.availible_cas9 = [""]

        self.initiateUI()
        self.create_menubar()
        self.create_toolbar()

    # display widgets
    def display_gRNA(self):
        groupBox = QtWidgets.QGroupBox()

        grid = QtWidgets.QGridLayout()

        self.display_candidates_label = QtWidgets.QLabel("Candidate guide RNA's")
        self.display_candidates = QtWidgets.QTableView()
        self.display_backup_label = QtWidgets.QLabel("Backup guide RNA's")
        self.display_backup = QtWidgets.QTableView()
        self.display_dropped_label = QtWidgets.QLabel("Dropped guide RNA's")
        self.display_dropped = QtWidgets.QTableView()
        self.display_offtargets_label = QtWidgets.QLabel("Offtargets")
        self.display_offtargets = QtWidgets.QTableView()

        grid.addWidget(self.display_candidates_label, 0, 0)
        grid.addWidget(self.display_candidates, 1, 0)
        grid.addWidget(self.display_backup_label, 0, 1)
        grid.addWidget(self.display_backup, 1, 1)

        grid.addWidget(self.display_dropped_label, 2, 0)
        grid.addWidget(self.display_dropped, 3, 0)
        grid.addWidget(self.display_offtargets_label, 2, 1)
        grid.addWidget(self.display_offtargets, 3, 1)

        groupBox.setLayout(grid)

        return groupBox

    # main window
    def initiateUI(self):
        self.setWindowTitle("Crispinator: The MMI guide RNA database")
        self.statusBar().showMessage('Ready')
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "icon.png")))
        QtWidgets.QToolTip.setFont(QtGui.QFont('Times New Roman', 12))

        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)

        layout = QtWidgets.QGridLayout(central_widget)

        layout.addWidget(self.display_gRNA(), 0, 0)

        self.main_progressbar = QtWidgets.QProgressBar()
        self.main_progressbar.setMaximum(100)
        self.main_progressbar_value = 0
        self.main_progressbar.setValue(self.main_progressbar_value)
        self.statusBar().addPermanentWidget(self.main_progressbar)
        self.showMaximized()
        self.show()

    # functions
    def create_menubar(self):
        """
        actions for the menubar
        """
        menu = self.menuBar()
        # Database
        databasemenu = menu.addMenu('Database')

        createDatabaseDialog = QtWidgets.QAction("Create database", self)
        createDatabaseDialog.setShortcut("Ctrl+A")
        createDatabaseDialog.triggered.connect(lambda: self.database_functions(trigger="create_db"))

        deleteDatabaseDialog = QtWidgets.QAction("Remove database", self)
        deleteDatabaseDialog.setShortcut("Ctrl+S")
        deleteDatabaseDialog.triggered.connect(lambda: self.database_functions(trigger="delete_db"))

        importDatabaseAction = QtWidgets.QAction('Database', self)
        importDatabaseAction.setShortcut("Ctrl+T")
        importDatabaseAction.triggered.connect(self.import_database_function)

        exportDatabasesAction = QtWidgets.QAction('Databases', self)
        exportDatabasesAction.setShortcut("Ctrl+Y")
        exportDatabasesAction.triggered.connect(lambda: self.exportfiles(trigger="database"))

        sqlDialog = QtWidgets.QAction("SQL interface", self)
        sqlDialog.setShortcut("Ctrl+D")
        sqlDialog.triggered.connect(lambda: self.database_functions(trigger="sql"))

        databasemenu.addAction(createDatabaseDialog)
        databasemenu.addAction(deleteDatabaseDialog)
        databasemenu.addAction(importDatabaseAction)
        databasemenu.addAction(exportDatabasesAction)
        databasemenu.addAction(sqlDialog)

        # Guide RNA
        grnamenu = menu.addMenu('Guide RNA')
        showGrnaStatsDialog = QtWidgets.QAction("Statistics", self)
        showGrnaStatsDialog.setShortcut("Ctrl+G")
        showGrnaStatsDialog.triggered.connect(lambda: self.grna_functions(trigger="stats"))

        predictGrnaDialog = QtWidgets.QAction("Search", self)
        predictGrnaDialog.setShortcut("Ctrl+H")
        predictGrnaDialog.triggered.connect(lambda: self.grna_functions(trigger="search"))

        displayGrnaDialog = QtWidgets.QAction("Display", self)
        displayGrnaDialog.setShortcut("Ctrl+J")
        displayGrnaDialog.triggered.connect(lambda: self.grna_functions(trigger="display"))

        exportGuideRNA = QtWidgets.QAction("Guide RNA's", self)
        exportGuideRNA.setShortcut("Ctrl+U")
        exportGuideRNA.triggered.connect(self.export_guide_rna)

        grnamenu.addAction(showGrnaStatsDialog)
        grnamenu.addAction(predictGrnaDialog)
        grnamenu.addAction(displayGrnaDialog)
        grnamenu.addAction(exportGuideRNA)

        # Navigate
        navMenu = menu.addMenu('External')
        mycobrowser = QtWidgets.QAction("Mycobrowser", self)
        mycobrowser.triggered.connect(lambda: self.navigate_functions(trigger="mycobrowser"))
        chopchop = QtWidgets.QAction("ChopChop", self)
        chopchop.triggered.connect(lambda: self.navigate_functions(trigger="chopchop"))
        ncbi = QtWidgets.QAction("NCBI", self)
        ncbi.triggered.connect(lambda: self.navigate_functions(trigger="ncbi"))
        git = QtWidgets.QAction("Source code", self)
        git.triggered.connect(lambda: self.navigate_functions(trigger="git"))

        navMenu.addAction(mycobrowser)
        navMenu.addAction(chopchop)
        navMenu.addAction(ncbi)
        navMenu.addAction(git)

        # about
        aboutMenu = menu.addMenu('About')
        showLicenceDialog = QtWidgets.QAction("licence", self)
        showLicenceDialog.triggered.connect(lambda: self.about_functions(trigger="license"))

        aboutMenu.addAction(showLicenceDialog)

    def create_toolbar(self):
        toolbar = QtWidgets.QToolBar("Main")
        toolbar.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon | QtCore.Qt.AlignLeading)
        toolbar.setIconSize(QtCore.QSize(48, 48))
        self.addToolBar(toolbar)

        showGrnaStatsIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "metadata.png"))
        showGrnaStats = QtWidgets.QAction(showGrnaStatsIcon, "View metadata", self)
        showGrnaStats.setToolTip("Displays data associated with the guide RNA's currently in the database")
        showGrnaStats.setShortcut("Ctrl+G")
        showGrnaStats.triggered.connect(lambda: self.grna_functions(trigger="stats"))

        predictNewGrnaIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "search.png"))
        predictNewGrna = QtWidgets.QAction(predictNewGrnaIcon, "Predict gRNA", self)
        predictNewGrna.setToolTip("Predict guide RNA's")
        predictNewGrna.setShortcut("Ctrl+H")
        predictNewGrna.triggered.connect(lambda: self.grna_functions(trigger="search"))

        displayGrnaIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "display.png"))
        displayGrna = QtWidgets.QAction(displayGrnaIcon, "CrisprI", self)
        displayGrna.setToolTip("Query database and display guide RNA's for Crispr Interference")
        displayGrna.setShortcut("Ctrl+J")
        displayGrna.triggered.connect(lambda: self.grna_functions(trigger="display"))

        cleargRNATablesIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "clear.png"))
        cleargRNATables = QtWidgets.QAction(cleargRNATablesIcon, "Clear tables", self)
        cleargRNATables.setToolTip("Clear all the tables")
        cleargRNATables.triggered.connect(lambda: self.grna_functions(trigger="clear"))

        exportGrnaIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "Export_grna.png"))
        exportGrna = QtWidgets.QAction(exportGrnaIcon, "Export gRNA", self)
        exportGrna.setToolTip("Export guide RNA and related info")
        exportGrna.setShortcut("Ctrl+U")
        exportGrna.triggered.connect(self.export_guide_rna)

        addDataBaseIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "add_database.png"))
        addDataBase = QtWidgets.QAction(addDataBaseIcon, "Add orgnanism", self)
        addDataBase.setToolTip("Add new organism database, required to search guide RNA's")
        addDataBase.setShortcut("Ctrl+A")
        addDataBase.triggered.connect(lambda: self.database_functions(trigger="create_db"))

        databaseImportIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "import_database.png"))
        databaseImport = QtWidgets.QAction(databaseImportIcon, "Import database", self)
        databaseImport.setToolTip("Import externally generated database")
        databaseImport.setShortcut("Ctrl+T")
        databaseImport.triggered.connect(self.import_database_function)

        removeDataBaseIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "remove_database.png"))
        removeDataBase = QtWidgets.QAction(removeDataBaseIcon, "Remove orgnanism", self)
        removeDataBase.setToolTip("Removes an organism database")
        removeDataBase.setShortcut("Ctrl+S")
        removeDataBase.triggered.connect(lambda: self.database_functions(trigger="delete_db"))

        databaseExportIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "export_database.png"))
        databaseExport = QtWidgets.QAction(databaseExportIcon, "Export database", self)
        databaseExport.setToolTip("Export database, give it to a friend")
        databaseExport.setShortcut("Ctrl+Y")
        databaseExport.triggered.connect(lambda: self.exportfiles(trigger="database"))

        SQLQueryIcon = QtGui.QIcon(os.path.join(self.root, "gui", "Icons", "sql.png"))
        SQLQuery = QtWidgets.QAction(SQLQueryIcon, "SQL", self)
        SQLQuery.setToolTip("Query database directly using SQL")
        SQLQuery.setShortcut("Ctrl+D")
        SQLQuery.triggered.connect(lambda: self.database_functions(trigger="sql"))

        toolbar.addAction(showGrnaStats)
        toolbar.addAction(predictNewGrna)
        toolbar.addAction(displayGrna)
        toolbar.addAction(exportGrna)
        toolbar.addAction(cleargRNATables)
        toolbar.addAction(addDataBase)
        toolbar.addAction(removeDataBase)
        toolbar.addAction(databaseImport)
        toolbar.addAction(databaseExport)
        toolbar.addAction(SQLQuery)

    # UI behaviour menubar
    def fileclose(self):
        self.close()

    def import_database_function(self):
        importer = ImportDatabase()

        if importer.exec_():
            files = importer.out()
            sanity = [files['database'], files['bsgenome']]
            if all(sanity):
                shutil.copy(files['database'], os.path.join(self.root, "databases"))
                self.availible_databases.clear()
                self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                            os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]

                shutil.copy(files['bsgenome'], os.path.join(self.root, "src", "local_packages"))


                self.availible_bsgenome.clear()
                self.availible_bsgenome = [items for items in
                                           os.listdir(os.path.join(self.root, "src", "local_packages"))]


                worker = Install_New_R_Package(package_name=os.path.basename(files['bsgenome']))
                self.threadingPool.start(worker)

                while self.threadingPool.activeThreadCount() == 1:
                    QtWidgets.QApplication.processEvents()
                    self.statusBar().showMessage("Creating BSgenome R package ...")
                    QtWidgets.QApplication.processEvents()
                    if self.main_progressbar_value < 99:
                        self.main_progressbar_value += 1
                        self.main_progressbar.setValue(self.main_progressbar_value)
                        time.sleep(0.5)

                if self.threadingPool.waitForDone():
                    while self.main_progressbar_value < 100:
                        self.main_progressbar_value += 1
                        self.main_progressbar.setValue(self.main_progressbar_value)

                    QtWidgets.QMessageBox.about(self, "Info", "Files copied and installed successfully")
                    self.main_progressbar_value = 0
                    self.main_progressbar.setValue(self.main_progressbar_value)

            else:
                if not sanity[0]:
                    QtWidgets.QMessageBox.about(self, "Error", "database file missing")
                    return None

                if not sanity[1]:
                    QtWidgets.QMessageBox.about(self, "Error", "bsgenome file missing")
                    return None

    def exportfiles(self, trigger):
        if trigger == "database":
            exporter = ExportDatabase(self)
            if exporter.exec_():
                destination = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory"))
                database, bsgenome = exporter.files()
                shutil.copy(os.path.join(self.root, "databases", f"{database}.db"), destination)
                shutil.copy(os.path.join(self.root, "src", "local_packages", f"{bsgenome}"), destination)
                QtWidgets.QMessageBox.about(self, "Info", "Files copied successfully")

    def export_guide_rna(self):
        if self.database_querried:
            grna_exporter = ExportGuideRNA(candidates=self.candidate_gRNA_df,
                                           backup=self.backup_gRNA_df,
                                           dropped=self.dropped_gRNA_df,
                                           offtargets=self.offtarget_df)

            grna_exporter.exec_()
        else:
            QtWidgets.QMessageBox.about(self, "Error", "No data to export")

    def create_database_runner(self):
        """
        Everything neccesary to trigger the database dialog and run the functions to create a new database
        """
        create_db_dialog = CreateDatabaseDialog()
        if create_db_dialog.exec_():

            self.statusBar().showMessage("Entering database dialog")

            filepaths = create_db_dialog.filepaths_out()
            bsgenome_metadata = create_db_dialog.metadata_out()

            os.makedirs(os.path.join(self.root, "temp", "fastafiles"), exist_ok=True)

            self.statusBar().showMessage("Checking common errors")
            if not bool(filepaths["genome"]):
                QtWidgets.QMessageBox.about(self, "Error", "No genome fasta file detected")
                return None

            if not bool(filepaths["multifasta"]):
                QtWidgets.QMessageBox.about(self, "Error", "No multigene fasta file detected")
                return None

            if not bool(filepaths["gff"]):
                QtWidgets.QMessageBox.about(self, "Error", "No gff file detected")
                return None

            org = str(bsgenome_metadata['org_name']).split(" ")
            if len(org) < 3:
                QtWidgets.QMessageBox.about(self, "Error",
                                            "The organism needs to have 3 fields: <b>genus</b> <b>species</b> <b>strain</b> separated by spaces")
                return None

            self.statusBar().showMessage("Creating SQL database ...")
            self.main_progressbar_value = 10
            self.main_progressbar.setValue(self.main_progressbar_value)
            create_db_runner = CreateDatabaseMethod(filepaths_dict=filepaths, metadata_dict=bsgenome_metadata)
            create_db_runner.create()
            self.main_progressbar_value += 10
            self.main_progressbar.setValue(self.main_progressbar_value)

            bsgenome_worker = BSgenome_worker()
            self.threadingPool.start(bsgenome_worker)

            while self.threadingPool.activeThreadCount() == 1:
                QtWidgets.QApplication.processEvents()
                self.statusBar().showMessage("Creating BSgenome R package .")
                self.statusBar().showMessage("Creating BSgenome R package ..")
                self.statusBar().showMessage("Creating BSgenome R package ...")
                if self.main_progressbar_value < 99:
                    self.main_progressbar_value += 1
                    self.main_progressbar.setValue(self.main_progressbar_value)
                    time.sleep(0.9)

            if self.threadingPool.waitForDone():
                self.availible_databases.clear()
                self.availible_bsgenome.clear()

                self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                            os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
                self.availible_bsgenome = [str(os.path.splitext(items)[0]) for items in
                                           os.listdir(os.path.join(self.root, "src", "local_packages"))]

                while self.main_progressbar_value < 100:
                    self.main_progressbar_value += 1
                    self.main_progressbar.setValue(self.main_progressbar_value)
                    time.sleep(0.01)

                self.main_progressbar_value = 0
                self.main_progressbar.setValue(self.main_progressbar_value)

                if "temp" in os.listdir(self.root):
                    shutil.rmtree(os.path.join(self.root, "temp"))

                self.statusBar().showMessage("Ready")

    def delete_database_runner(self):
        """
        Opens and runs dialog to delete database
        """
        delete_db_dialog = DeleteDatabaseDialog(database_list=self.availible_databases,
                                                bsgenome_list=self.availible_bsgenome)

        if delete_db_dialog.exec_():
            self.statusBar().showMessage("Deleting database ...")
            items_to_delete = delete_db_dialog.out()
            self.availible_databases.clear()
            self.availible_bsgenome.clear()

            database = items_to_delete.get('database', None)
            database = database.replace(" ", "_")
            bsgenome = items_to_delete.get('bsgenome', None)

            if ".tar" in bsgenome:
                bsgenome = bsgenome.replace(".tar", "")

            os.remove(os.path.join(self.root, "databases", f"{database}.db"))
            os.remove(os.path.join(self.root, "src", "local_packages", f"{bsgenome}.tar.gz"))

            self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                        os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
            self.availible_bsgenome = [str(os.path.splitext(items)[0]) for items in
                                       os.listdir(os.path.join(self.root, "src", "local_packages"))]

            while self.main_progressbar_value < 100:
                self.main_progressbar_value += 1
                self.main_progressbar.setValue(self.main_progressbar_value)

            self.main_progressbar_value = 0
            self.main_progressbar.setValue(self.main_progressbar_value)
            self.statusBar().showMessage("Ready")

    def sql_dialog_function(self):
        """
        Opens a dialog that can be used to querry the SQL database directly
        """
        if not "temp" in os.listdir(self.root):
            tempdir = os.path.join(self.root, "temp")
            os.mkdir(tempdir)
        else:
            tempdir = os.path.join(self.root, "temp")

        sql_dialog = SQLDialog(database_list=self.availible_databases)

        if sql_dialog.exec_():
            self.statusBar().showMessage("SQL editor active")

        shutil.rmtree(tempdir)
        self.statusBar().showMessage("Ready")

    def database_functions(self, trigger):
        """
        Better as switch statement?
        """
        if trigger == "create_db":
            self.create_database_runner()

        if trigger == "delete_db":
            self.delete_database_runner()

        if trigger == "sql":
            self.sql_dialog_function()

    def navigate_functions(self, trigger):
        if trigger == "mycobrowser":
            webbrowser.open_new_tab("https://mycobrowser.epfl.ch/")

        if trigger == "chopchop":
            webbrowser.open_new_tab("https://chopchop.cbu.uib.no/")

        if trigger == "ncbi":
            webbrowser.open_new_tab("https://www.ncbi.nlm.nih.gov/")

        if trigger == "git":
            webbrowser.open_new_tab("https://github.com/JamesGallant/CrispR")

    ############ RUNNER FUNCTIONS #########

    def grna_stats_runner(self):
        """Display stats
        """
        grna_stats = StatsDialog(database_list=self.availible_databases)
        grna_stats.exec()

    def grna_search_runner(self):
        """
        search guide RNA's
        :return:
        """
        search_grna = SearchGrnaDialog(database_list=self.availible_databases, cas9_list=self.availible_cas9)
        if search_grna.exec_():
            user_data = search_grna.out()

            chosen_org = user_data['organism']
            is_circular = user_data['chromosome']
            mismatch = user_data['mismatch']
            pam = user_data['pam']
            tax_id = user_data['taxonomy_id']
            cores = user_data['cores']

            if is_circular not in {"TRUE", "FALSE"}:
                QtWidgets.QMessageBox.about(self, "Error", "chromosome needs to be TRUE or FALSE")
                return None

            sql_runner = SQL(database=chosen_org.replace(" ", "_"))
            searched_pams = sql_runner.list_pams()

            if pam in searched_pams:
                QtWidgets.QMessageBox.about(self, "Error",
                                            f"PAM {pam} has already been searched, choose a different pam")
                return None

            if tax_id == "":
                QtWidgets.QMessageBox.about(self, "Error",
                                            f"Taxonomy id required")
                return None

            self.statusBar().showMessage("Creating guide RNA database...")
            self.main_progressbar_value += 1
            self.main_progressbar.setValue(self.main_progressbar_value)

            if not "temp" in os.listdir(self.root):
                os.mkdir(os.path.join(self.root, "temp"))

            tempdir = os.path.join(self.root, "temp")

            if "global_gRNA" not in tempdir:
                os.mkdir(os.path.join(tempdir, "global_gRNA"))

            global_gRNA = os.path.join(tempdir, "global_gRNA")

            chosen_org_modified = chosen_org.replace(" ", "_")
            strain = chosen_org.split(" ")[-1]

            self.main_progressbar_value += 1
            self.main_progressbar.setValue(self.main_progressbar_value)

            tempfa = pd.read_sql("SELECT * FROM genome", sqlite3.connect(os.path.join(self.root,
                                                                                      "databases",
                                                                                      f"{chosen_org_modified}.db")))

            tempgff = pd.read_sql("SELECT * FROM gff_file", sqlite3.connect(os.path.join(self.root,
                                                                                         "databases",
                                                                                         f"{chosen_org_modified}.db")))

            tempgenes = pd.read_sql("SELECT * FROM genes", sqlite3.connect(os.path.join(self.root,
                                                                                        "databases",
                                                                                        f"{chosen_org_modified}.db")))

            tempgff.to_csv(f"{os.path.join(tempdir, chosen_org_modified)}.gff", header=False, index=False, sep="\t")

            fasta_header = [header for header in tempfa['header']]
            fasta_sequence = [sequence for sequence in tempfa['sequence']]

            self.main_progressbar_value += 1
            self.main_progressbar.setValue(self.main_progressbar_value)

            with open(os.path.join(tempdir, f"{chosen_org_modified}.fasta"), 'w') as input_fasta:
                for header, sequence in zip(fasta_header, fasta_sequence):
                    input_fasta.write(">" + header + '\n')
                    input_fasta.write(sequence + '\n')

            self.main_progressbar_value += 1
            self.main_progressbar.setValue(self.main_progressbar_value)

            fasta_header.clear()
            fasta_sequence.clear()

            fasta_header = [header for header in tempgenes['header']]
            fasta_sequence = [sequence for sequence in tempgenes['sequence']]

            with open(os.path.join(tempdir, f"{chosen_org_modified}_genes.fasta"), 'w') as input_fasta:
                for header, sequence in zip(fasta_header, fasta_sequence):
                    input_fasta.write(">" + header + '\n')
                    input_fasta.write(sequence + '\n')

            self.main_progressbar_value += 1
            self.main_progressbar.setValue(self.main_progressbar_value)

            fasta_header.clear()
            fasta_sequence.clear()

            for objects in self.availible_bsgenome:
                detected_package = objects.split("_")[0]
                detected_package = detected_package.split(".")[-1]
                if strain == detected_package:
                    bsgenome_package = str(objects.split("_")[0])

            genus, species, strain = chosen_org.split(" ")
            species = species.lower()

            config_file = {
                'organism': [f"{genus} {species} {strain}"],
                'taxonomy_id': [tax_id],
                'circular_chromosome': [is_circular],
                'input_file': [os.path.join(tempdir, f"{chosen_org_modified}_genes.fasta")],
                'gff_file': [os.path.join(tempdir, f"{chosen_org_modified}.gff")],
                'find_gRNA_with_cutsites': ["FALSE"],
                'find_paired_gRNA': ["FALSE"],
                'BSgenome': [bsgenome_package],
                'chromosomes_to_search': ["all"],
                'min_gap': [0],
                'max_gap': [20],
                'gRNA_size': [20],
                'max_mismatch_gRNA': [mismatch],
                'PAM_sequence': [pam],
                'PAM_length': [len(pam)],
                'n.cores': [cores],
                'scoring_method': ["CFDscore"]
            }

            config_file = pd.DataFrame(config_file)
            config_file.to_csv(os.path.join(tempdir, "config.txt"), index=False, sep="\t")

            findgRNA_worker = FindgRNA_worker()

            self.threadingPool.start(findgRNA_worker)

            while self.threadingPool.activeThreadCount() == 1:
                self.statusBar().showMessage("predicting all guide RNA's...")
                QtWidgets.QApplication.processEvents()
                if self.main_progressbar_value <= 90:
                    self.main_progressbar_value += 1
                    self.main_progressbar.setValue(self.main_progressbar_value)
                    time.sleep(2)

            if self.threadingPool.waitForDone():
                self.statusBar().showMessage("buidling global gRNA database...")
                database = Database(database=chosen_org_modified)
                database.create_gRNA_database(summary=os.path.join(global_gRNA, "Summary.xls"),
                                              offtarget=os.path.join(global_gRNA, "OfftargetAnalysis.txt"),
                                              config_file=os.path.join(tempdir, "config.txt"))

                while self.main_progressbar_value <= 100:
                    self.main_progressbar_value += 1
                    self.main_progressbar.setValue(self.main_progressbar_value)

                shutil.rmtree(tempdir)
                self.main_progressbar_value = 0
                self.main_progressbar.setValue(self.main_progressbar_value)
                self.statusBar().showMessage("Ready")

    def grna_display_runner(self):
        display_grna = DisplayGuideRNA(database_list=self.availible_databases, cas9_list=self.availible_cas9)

        if display_grna.exec_():

            holder = PandasModel(pd.DataFrame({'': []}))
            self.display_candidates.setModel(holder)
            self.display_backup.setModel(holder)
            self.display_dropped.setModel(holder)
            self.display_offtargets.setModel(holder)
            self.database_querried = False

            self.statusBar().showMessage("Preparing ...")
            self.main_progressbar_value += 1
            self.main_progressbar.setValue(self.main_progressbar_value)

            user_options = display_grna.out()
            if "temp" not in os.listdir(self.root):
                tempdir = os.path.join(self.root, "temp")
                os.mkdir(os.path.join(tempdir))
            else:
                tempdir = os.path.join(self.root, "temp")

            database = str(user_options['organism']).replace(" ", "_")
            mismatch = user_options['max_mismatch']
            max_grna = user_options['max_grna_count']
            max_primer_len = user_options['max_primer_len']
            user_cas9 = user_options['cas9']
            user_pam_tolerance = user_options['pam_tolerance']
            user_fiveprime = user_options['nucleotides_5']
            user_threeprime = user_options['nucleotides_3']

            gene_mask_dictionary = {
                'genes': [items.replace("_", "").lower() if "_" in items else items.lower() for items in
                          user_options['genes']],
                'masks': [items.replace("_", "").lower() if "_" in items else items.lower() for items in
                          user_options['masks']]
            }

            sqlrunner = SQL(database=database)
            headers = sqlrunner.custom_sql("SELECT header FROM genes").to_dict('list')

            gene_check = [True if gene in headers['header'] else False for gene in gene_mask_dictionary['genes']]
            mask_check = [True if gene in headers['header'] else False for gene in gene_mask_dictionary['masks']]

            for idx, val in enumerate(gene_check):
                if not val:
                    db = database.replace("_", " ")
                    QtWidgets.QMessageBox.about(self, "Error",
                                                f"{gene_mask_dictionary['genes'][idx]} was not found in {db}")

                    self.main_progressbar_value = 0
                    self.main_progressbar.setValue(self.main_progressbar_value)
                    return None

            for idx, val in enumerate(mask_check):
                if not val:
                    db = database.replace("_", " ")
                    QtWidgets.QMessageBox.about(self, "Error",
                                                f"{gene_mask_dictionary['masks'][idx]} was not found in {db}")

                    self.main_progressbar_value = 0
                    self.main_progressbar.setValue(self.main_progressbar_value)
                    return None

            if mismatch == "":
                QtWidgets.QMessageBox.about(self, "Error",
                                            "First search guide RNA's")

                self.main_progressbar_value = 0
                self.main_progressbar.setValue(self.main_progressbar_value)
                return None

            # Strand is r for reverse
            worker = CrisprInterference_worker(database=database,
                                               mismatch=mismatch,
                                               strand='r',
                                               max_grna=max_grna,
                                               genes_masks=gene_mask_dictionary,
                                               max_primer_size=max_primer_len,
                                               cas9_organism=user_cas9,
                                               pam_tolerance=user_pam_tolerance,
                                               fiveprime_nucleotides=user_fiveprime,
                                               threeprime_nucleotides=user_threeprime)

            self.threadingPool.start(worker)

            while self.threadingPool.activeThreadCount() == 1:
                self.statusBar().showMessage("Gathering guide RNA's...")
                QtWidgets.QApplication.processEvents()

                if self.main_progressbar_value < 90:
                    self.main_progressbar_value += 1
                    self.main_progressbar.setValue(self.main_progressbar_value)
                    time.sleep(0.8)

            if self.threadingPool.waitForDone():
                self.statusBar().showMessage("Gathering data ...")
                self.candidate_gRNA_df = pd.read_csv(
                    filepath_or_buffer=os.path.join(self.root, "temp", "candidates.txt"),
                    sep=",")
                self.backup_gRNA_df = pd.read_csv(filepath_or_buffer=os.path.join(self.root, "temp", "backup.txt"),
                                                  sep=",")
                self.dropped_gRNA_df = pd.read_csv(filepath_or_buffer=os.path.join(self.root, "temp", "dropped.txt"),
                                                   sep=",")
                self.offtarget_df = pd.read_csv(filepath_or_buffer=os.path.join(self.root, "temp", "offtargets.txt"),
                                                sep=",")

                cand_model, backup_model, dropped_model, offtargets_model = map(PandasModel, [self.candidate_gRNA_df,
                                                                                              self.backup_gRNA_df,
                                                                                              self.dropped_gRNA_df,
                                                                                              self.offtarget_df])

                while self.main_progressbar_value < 100:
                    self.main_progressbar_value += 1
                    self.statusBar().showMessage("Formatting for display...")
                    self.main_progressbar.setValue(self.main_progressbar_value)
                    time.sleep(0.01)

                self.display_candidates.setModel(cand_model)
                self.display_backup.setModel(backup_model)
                self.display_dropped.setModel(dropped_model)
                self.display_offtargets.setModel(offtargets_model)
                self.database_querried = True
                self.main_progressbar_value = 0
                self.main_progressbar.setValue(self.main_progressbar_value)
                self.statusBar().showMessage("Ready")

                hits = [genes for genes in self.candidate_gRNA_df['genes']]
                missed = list(set(gene_mask_dictionary['genes']) - set(hits))

                EOSpopup(missed_genes=missed).exec_()

            shutil.rmtree(tempdir)

    def grna_clear_runner(self):
        def _worker(dataframe=pd.DataFrame()):
            if dataframe is not None:
                empty = pd.DataFrame(dict.fromkeys(dataframe, []))
                return PandasModel(empty)

        self.display_candidates.setModel(_worker(self.candidate_gRNA_df))
        self.display_backup.setModel(_worker(self.backup_gRNA_df))
        self.display_dropped.setModel(_worker(self.dropped_gRNA_df))
        self.display_offtargets.setModel(_worker(self.offtarget_df))

    def grna_functions(self, trigger):
        switch = {
            'stats': self.grna_stats_runner,
            'search': self.grna_search_runner,
            'display': self.grna_display_runner,
            'clear': self.grna_clear_runner
        }
        switch.get(trigger, None)()

    def show_license_func(self):
        showLicense().exec_()

    def about_functions(self, trigger):
        if trigger == "license":
            self.show_license_func()

    # internal event modifiers here
    def closeEvent(self, event):
        user_response = QtWidgets.QMessageBox.question(self, 'Message', "Are you sure you would like to quit",
                                                       QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)

        if user_response == QtWidgets.QMessageBox.Yes:
            files_in_root = os.listdir(os.path.dirname(os.path.abspath(__name__)))
            if "temp" in files_in_root:
                shutil.rmtree(os.path.join(os.path.dirname(os.path.abspath(__name__)), "temp"))

            if "scores.RDS" in files_in_root:
                os.remove(os.path.join(os.path.dirname(os.path.abspath(__name__)), "scores.RDS"))

            if "offTargets.RDS" in files_in_root:
                os.remove(os.path.join(os.path.dirname(os.path.abspath(__name__)), "offTargets.RDS"))
            event.accept()
        else:
            event.ignore()


def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())

    gui = MainWindow()
    app.exec_()
    sys.exit()


if __name__ == '__main__':
    main()
