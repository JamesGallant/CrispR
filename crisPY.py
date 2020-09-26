# Statistics tab interactions are buggy, split helper function to account for different user choices
import sys
import qdarkstyle
import shutil
import sqlite3
from datetime import date

from gui.components.workers import FindgRNA_worker, BSgenome_worker, CrisprInterference_worker
from gui.components.dialogs import *
from gui.components.pandashandler import PandasModel
from gui.behaviour.main_functions import CreateDatabaseMethod
from utilities.database_tools import Database, SQL


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        self.root = os.path.dirname(os.path.abspath(__name__))
        super().__init__()
        self.progressbar_counter = 0
        self.gRNA_progressbar_counter = 0
        self.main_progressbar_value = 0
        self.gene_mask_dict = {
            'genes': [],
            'masks': []
        }

        self.filepaths = {
            'genome': "",
            'multifasta': "",
            'gff': "",
        }
        self.threadingPool = QtCore.QThreadPool()
        self.candidate_gRNA_df = None
        self.backup_gRNA_df = None
        self.dropped_gRNA_df = None
        self.offtarget_df = None
        self.database_querried = False

        try:
            self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                        os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
            self.availible_bsgenome = [str(os.path.splitext(items)[0]) for items in
                                       os.listdir(os.path.join(self.root, "R", "local_packages"))]
        except:
            self.availible_databases = [""]
            self.availible_bsgenome = [""]

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

        self.initiateUI()
        self.create_menubar()

    # main window
    def initiateUI(self):
        self.setWindowTitle("CrisPY: The MMI guide RNA database")
        self.statusBar().showMessage('Ready')
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.root, "gui", "metadata", "icon.png")))
        QtWidgets.QToolTip.setFont(QtGui.QFont('Times New Roman', 12))

        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)

        layout = QtWidgets.QGridLayout(central_widget)

        # layout.addWidget(self.availibleData_controls(), 0, 0)
        # layout.addWidget(self.create_new_db_controls(), 0, 1)
        # layout.addWidget(self.predict_gRNA_controls(), 1, 0)
        # layout.addWidget(self.process_gRNA_data(), 1, 1)
        layout.addWidget(self.display_gRNA(), 0, 0)

        self.main_progressbar = QtWidgets.QProgressBar()
        self.main_progressbar.setMaximum(100)
        self.main_progressbar_value = 0
        self.main_progressbar.setValue(self.main_progressbar_value)
        self.statusBar().addPermanentWidget(self.main_progressbar)
        self.showMaximized()
        self.show()

    # UI menubar
    def create_menubar(self):
        """
        actions for the menubar
        """
        menu = self.menuBar()
        # File
        filemenu = menu.addMenu('File')
        importDatabaseAction = QtWidgets.QAction('Database', self)
        importDatabaseAction.setShortcut("Ctrl+T")
        importDatabaseAction.triggered.connect(self.import_database_function)

        exportDatabasesAction = QtWidgets.QAction('Databases', self)
        exportDatabasesAction.setShortcut("Ctrl+Y")
        exportDatabasesAction.triggered.connect(lambda: self.exportfiles(trigger="database"))
        exportGuideRNA = QtWidgets.QAction("Guide RNA's", self)
        exportGuideRNA.setShortcut("Ctrl+U")
        exportGuideRNA.triggered.connect(self.export_guide_rna)

        exitAction = QtWidgets.QAction('Exit', self)
        exitAction.setShortcut("esc")
        exitAction.triggered.connect(self.fileclose)

        import_group = filemenu.addMenu('Import')
        import_group.addAction(importDatabaseAction)
        export_group = filemenu.addMenu('Export')
        export_group.addAction(exportDatabasesAction)
        export_group.addAction(exportGuideRNA)
        filemenu.addAction(exitAction)

        # Database
        databasemenu = menu.addMenu('Database')

        createDatabaseDialog = QtWidgets.QAction("Create database", self)
        createDatabaseDialog.setShortcut("Ctrl+A")
        createDatabaseDialog.triggered.connect(lambda: self.database_functions(trigger="create_db"))

        deleteDatabaseDialog = QtWidgets.QAction("Remove database", self)
        deleteDatabaseDialog.setShortcut("Ctrl+S")
        deleteDatabaseDialog.triggered.connect(lambda: self.database_functions(trigger="delete_db"))

        sqlDialog = QtWidgets.QAction("SQL interface", self)
        sqlDialog.setShortcut("Ctrl+D")
        sqlDialog.triggered.connect(lambda: self.database_functions(trigger="sql"))

        databasemenu.addAction(createDatabaseDialog)
        databasemenu.addAction(deleteDatabaseDialog)
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

        grnamenu.addAction(showGrnaStatsDialog)
        grnamenu.addAction(predictGrnaDialog)
        grnamenu.addAction(displayGrnaDialog)

        # Help
        helpmenu = menu.addMenu('Help')

    # UI behaviour menubar
    def fileclose(self):
        self.close()

    def import_database_function(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, "Add database", ".", "Database files (*.db)")
        file_to_copy = filename[0]
        shutil.copy(file_to_copy, os.path.join(self.root, "databases"))
        self.availible_databases.clear()
        self.current_databases.clear()
        self.search_databases.clear()
        self.processed_organism.clear()

        self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                    os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
        self.current_databases.addItems(self.availible_databases)
        self.search_databases.addItems(self.availible_databases)
        self.processed_organism.addItems(self.availible_databases)

    def exportfiles(self, trigger):
        if trigger == "database":
            exporter = ExportDatabase(self)
            if exporter.exec_():
                # self.statusBar.showMessage("Copying files...")
                destination = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory"))
                files_to_copy = exporter.files()
                if files_to_copy == "all":
                    for files in files_to_copy:
                        shutil.copy(os.path.join(self.root, "databases", f"{files}.db"), destination)
                else:
                    shutil.copy(os.path.join(self.root, "databases", f"{files_to_copy}.db"), destination)

    def export_guide_rna(self):
        if self.database_querried:
            grna_exporter = ExportGuideRNA(candidates=self.candidate_gRNA_df,
                                           backup=self.backup_gRNA_df,
                                           dropped=self.dropped_gRNA_df,
                                           offtargets=self.offtarget_df)

            if grna_exporter.exec_():
                to_export = grna_exporter
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

            try:
                genus, species, strain = str(bsgenome_metadata['org_name']).split(" ")
                del genus
                del species
                del strain
            except:
                QtWidgets.QMessageBox.about(self, "Error",
                                            "The organism needs to have 3 fields: <b>genus</b> <b>species</b> <b>strain</b> separated by spaces")

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
                self.availible_pams.clear()
                self.current_databases.clear()
                self.search_databases.clear()
                self.availible_databases.clear()
                self.availible_bsgenome.clear()

                self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                            os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
                self.availible_bsgenome = [str(os.path.splitext(items)[0]) for items in
                                           os.listdir(os.path.join(self.root, "R", "local_packages"))]
                self.current_databases.addItems(self.availible_databases)
                self.search_databases.addItems(self.availible_databases)

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
            self.current_databases.clear()
            self.search_databases.clear()
            self.availible_databases.clear()
            self.availible_bsgenome.clear()

            database = items_to_delete.get('database', None)
            database = database.replace(" ", "_")
            bsgenome = items_to_delete.get('bsgenome', None)

            if ".tar" in bsgenome:
                bsgenome = bsgenome.replace(".tar", "")

            os.remove(os.path.join(self.root, "databases", f"{database}.db"))
            os.remove(os.path.join(self.root, "R", "local_packages", f"{bsgenome}.tar.gz"))

            self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                        os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
            self.availible_bsgenome = [str(os.path.splitext(items)[0]) for items in
                                       os.listdir(os.path.join(self.root, "R", "local_packages"))]
            self.current_databases.addItems(self.availible_databases)
            self.search_databases.addItems(self.availible_databases)

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
        search_grna = SearchGrnaDialog(database_list=self.availible_databases)
        if search_grna.exec_():
            user_data = search_grna.out()

            chosen_org = user_data['organism']
            is_circular = user_data['chromosome']
            mismatch = user_data['mismatch']
            pam = user_data['pam']

            if is_circular not in {"TRUE", "FALSE"}:
                QtWidgets.QMessageBox.about(self, "Error", "chromosome needs to be TRUE or FALSE")
                return None

            sql_runner = SQL(database=chosen_org.replace(" ", "_"))
            searched_pams = sql_runner.list_pams()

            if pam in searched_pams:
                QtWidgets.QMessageBox.about(self, "Error",
                                            f"PAM {pam} has already been searched, choose a different pam")
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
        display_grna = DisplayGuideRNA(database_list=self.availible_databases)

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

            gene_mask_dictionary = {
                'genes': [items.replace("_", "") if "_" in items else items for items in user_options['genes']],
                'masks': [items.replace("_", "") if "_" in items else items for items in user_options['masks']]
            }
            # Strand is r for reverse
            worker = CrisprInterference_worker(database=database,
                                               mismatch=mismatch,
                                               strand='r',
                                               max_grna=max_grna,
                                               genes_masks=gene_mask_dictionary,
                                               max_primer_size=max_primer_len)

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

            shutil.rmtree(tempdir)

    def grna_functions(self, trigger):
        if trigger == "stats":
            self.grna_stats_runner()

        if trigger == "search":
            self.grna_search_runner()

        if trigger == "display":
            self.grna_display_runner()

    # old funcs start here
    ########################################################################################################################################3
    # UI metadata controls
    def availibleData_controls(self):
        """
        Metadata availible including database and predicted gRNA's
        """

        groupBox = QtWidgets.QGroupBox("Availible guide RNA databases")

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

    # ui functions: metadata behaviour
    def _update_metadata(self, database, user_pam, user_mismatch):
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
        self.metadata_mismatches_value = database.list_mismatches(pam=user_pam)
        self.metadata1_value = database.gene_count()
        self.metadata2_value = database.gRNA_count()
        self.metadata3_value = database.genes_with_gRNA()

        if user_mismatch not in self.metadata_mismatches_value:
            user_mismatch = self.metadata_mismatches_value[0]

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

        self.sql_runner = SQL(database=str(self.search_databases.itemText(value)).replace(" ", "_"))

        self._update_metadata(database=self.sql_runner, user_pam=pams_value,
                              user_mismatch=mismatch_value)

    def update_by_current_pams(self, value):
        """
        Update metadata based on pams
        """
        organism = str(self.current_databases.currentText())
        pams_value = self.availible_pams[value]
        mismatch_value = str(self.metadata_mismatches.currentText())

        self.sql_runner = SQL(database=organism.replace(" ", "_"))

        self._update_metadata(database=self.sql_runner, user_pam=pams_value,
                              user_mismatch=mismatch_value)

    def update_by_current_mismatch(self, value):
        """
        Update metadata based on mismatches
        """
        organism = str(self.current_databases.currentText())
        pams_value = str(self.current_pams.currentText())
        mismatch_value = self.metadata_mismatches_value[value]

        self.sql_runner = SQL(database=organism.replace(" ", "_"))

        self._update_metadata(database=self.sql_runner, user_pam=pams_value,
                              user_mismatch=mismatch_value)

    # UI new database
    def create_new_db_controls(self):
        """
        This essentially controls the making of the BSgenome and creating the first sqlite3 database
        """
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

        self.submit_newdb_btn = QtWidgets.QPushButton('submit', self)
        self.submit_newdb_btn.setToolTip('Create the database')
        self.submit_newdb_btn.clicked.connect(self._run_db_creation)

        self.submission_progress = QtWidgets.QProgressBar(self)
        self.submission_progress.setMaximum(100)

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
        grid.addWidget(self.submit_newdb_btn, 5, 0, 1, 2)
        grid.addWidget(self.submission_progress, 5, 2, 1, 7)

        groupBox.setLayout(grid)

        return groupBox

    # ui behaviour: create database functions
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
                self.genome_fasta_uploadbtn.setStyleSheet("background: green; color: white")
            else:
                self.genome_fasta_uploadbtn.setText("Failed")
                self.genome_fasta_uploadbtn.setStyleSheet("background: red; color: white")

        if bool(self.filepaths["multifasta"]):
            gene_filepath = self.filepaths["multifasta"]
            gene_file, gene_extention = os.path.splitext(gene_filepath)
            gene_file = gene_file.split("/")[-1]

            self.gene_fasta_uploaddescription.setText(f"{gene_file} is loaded")

            if gene_extention == ".fa" or gene_extention == ".fasta":
                self.gene_fasta_uploadbtn.setText("success")
                self.gene_fasta_uploadbtn.setStyleSheet("background: green; color: white")
            else:
                self.gene_fasta_uploadbtn.setText("Failed")
                self.gene_fasta_uploadbtn.setStyleSheet("background: red; color: white")

        if bool(self.filepaths["gff"]):
            gff_filepath = self.filepaths["gff"]
            gff_file, gff_extention = os.path.splitext(gff_filepath)
            gff_file = gff_file.split("/")[-1]

            self.gff_uploaddescription.setText(f"{gff_file} is loaded")

            if gff_extention == ".gff":
                self.gff_uploadbtn.setText("success")
                self.gff_uploadbtn.setStyleSheet("background: green; color: white")
            else:
                self.gff_uploadbtn.setText("Failed")
                self.gff_uploadbtn.setStyleSheet("background: red; color: white")

    # UI predict guide RNA
    def predict_gRNA_controls(self):
        """
        This runs the python wrapper for the predict_gRNA script, needs to update the database
        """
        groupBox = QtWidgets.QGroupBox("Predict guide RNA")

        self.PAM_label = QtWidgets.QLabel("PAM")
        self.PAM_input = QtWidgets.QLineEdit("AGAAG")
        self.PAM_input.setToolTip("Choose the PAM")
        self.PAM_input.textChanged.connect(self.get_pam)

        self.PAM_scoring_label = QtWidgets.QLabel("Scoring metircs")
        self.PAM_scoring_CFD = QtWidgets.QRadioButton("CFD score")
        self.PAM_scoring_CFD.setEnabled(False)
        self.PAM_scoring_CFD.setToolTip("Only availible for <b>NGG</b> PAM")
        self.PAM_scoring_HSU_Zhang = QtWidgets.QRadioButton("Hsu-Zhang")
        self.PAM_scoring_HSU_Zhang.setToolTip("Only availible for <b>NGG</b> PAM")
        self.PAM_scoring_HSU_Zhang.setEnabled(False)

        self.PAM_scoring_group = QtWidgets.QButtonGroup()
        self.PAM_scoring_group.addButton(self.PAM_scoring_CFD)
        self.PAM_scoring_group.addButton(self.PAM_scoring_HSU_Zhang)
        self.PAM_scoring_group.setExclusive(True)

        self.mismatch_label = QtWidgets.QLabel("Mismatches")
        self.mismatch = QtWidgets.QSpinBox()
        self.mismatch.setMinimum(0)
        self.mismatch.setMaximum(15)
        self.mismatch.setValue(4)
        self.mismatch.setToolTip(
            "Number of base mismathces in the off targets, use 0 to find guide RNA's with no off targets")

        self.circ_chromosome_label = QtWidgets.QLabel("Circular")
        self.circ_chromosome = QtWidgets.QLineEdit("TRUE")
        self.circ_chromosome.setToolTip("Is the crhomosome cirular, <b>TRUE</b> for bacteria")

        self.search_databases_label = QtWidgets.QLabel("Organism")
        self.search_databases = QtWidgets.QComboBox()
        self.search_databases.addItems(self.availible_databases)
        self.search_databases.setToolTip("Choose your organism")

        self.submit_gRNA_predict = QtWidgets.QPushButton('submit', self)
        self.submit_gRNA_predict.setToolTip("Predict all guide RNA's")
        self.submit_gRNA_predict.clicked.connect(self._run_gRNA_prediction)

        self.gRNA_progress = QtWidgets.QProgressBar(self)
        self.gRNA_progress.setMaximum(100)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.PAM_label, 0, 0)
        grid.addWidget(self.PAM_input, 0, 1)
        grid.addWidget(self.PAM_scoring_label, 0, 2)
        grid.addWidget(self.PAM_scoring_CFD, 0, 3)
        grid.addWidget(self.PAM_scoring_HSU_Zhang, 0, 4)
        grid.addWidget(self.circ_chromosome_label, 1, 0)
        grid.addWidget(self.circ_chromosome, 1, 1)
        grid.addWidget(self.mismatch_label, 1, 2)
        grid.addWidget(self.mismatch, 1, 3, 1, 3)
        grid.addWidget(self.search_databases_label, 2, 0)
        grid.addWidget(self.search_databases, 2, 1, 1, 2)
        grid.addWidget(self.submit_gRNA_predict, 2, 3, 1, 4)
        grid.addWidget(self.gRNA_progress, 3, 1, 1, 9)

        groupBox.setLayout(grid)

        return groupBox

    # ui behaviour: predict gRNA functions
    def get_pam(self, pam):
        """
        gets the PAM chosen by user and processes it
        """
        is_enabled = True if pam == "NGG" else False
        self.PAM_scoring_HSU_Zhang.setEnabled(is_enabled)
        self.PAM_scoring_CFD.setEnabled(is_enabled)

    # UI process guide RNA
    def process_gRNA_data(self):
        """
        The ui for extracting the correct guide RNA's
        """
        groupBox = QtWidgets.QGroupBox("Extract guide RNA's")

        self.processed_import_dialog_label = QtWidgets.QLabel("Target genes")
        self.processed_import_dialog = QtWidgets.QPushButton("Import genelist", self)
        self.processed_import_dialog.setToolTip("Upload genes and/or masks")
        self.processed_import_dialog.clicked.connect(self.processed_dialog_functions)

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

        self.processed_submit = QtWidgets.QPushButton("submit")
        self.processed_submit.clicked.connect(self.run_gRNA_processing)
        self.processed_progressbar = QtWidgets.QProgressBar(self)
        self.processed_progressbar.setMaximum(100)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.processed_import_dialog_label, 0, 0)
        grid.addWidget(self.processed_primerlen_label, 0, 1)
        grid.addWidget(self.processed_import_dialog, 1, 0)
        grid.addWidget(self.processed_primerlen, 1, 1)
        grid.addWidget(self.processed_addnucleotides_5_label, 3, 0)
        grid.addWidget(self.processed_addnucleotides_5, 4, 0)
        grid.addWidget(self.processed_addnucleotides_3_label, 3, 1)
        grid.addWidget(self.processed_addnucleotides_3, 4, 1)
        grid.addWidget(self.processed_organism_label, 5, 0)
        grid.addWidget(self.processed_maxMismatch_label, 5, 1)
        grid.addWidget(self.processed_organism, 6, 0)
        grid.addWidget(self.processed_maxMismatch, 6, 1)
        grid.addWidget(self.processed_maxPrimers_label, 7, 0)
        grid.addWidget(self.processed_maxPrimers, 8, 0)
        grid.addWidget(self.processed_submit, 8, 1)
        grid.addWidget(self.processed_progressbar, 9, 0, 1, 2)
        groupBox.setLayout(grid)

        return groupBox

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

    # ui behaviour
    # ui functions: process gRNA behaviour

    def processed_dialog_functions(self):
        """
        Processes the dialog click events to upload data
        """
        gene_importer = ImportGenes(self)
        if gene_importer.exec_():
            self.gene_mask_dict = gene_importer.out()
            if bool(self.gene_mask_dict['genes']):
                self.processed_import_dialog.setText("success")
                self.processed_import_dialog.setStyleSheet("background: green; color: white")

    def update_processed_mismatches(self, value):
        self.processed_maxMismatch.clear()
        self.processeed_mismatches_value.clear()
        sqlrunner = SQL(database=str(self.processed_organism.itemText(value).replace(" ", "_")))
        self.processeed_mismatches_value = sqlrunner.list_mismatches_processing()
        self.processed_maxMismatch.addItems(self.processeed_mismatches_value)

    # Submit button functions
    # Create database
    def _run_db_creation(self):
        """
        Build and create the database and R functions here
        """
        self.statusBar().showMessage("Creating genome metadata")

        os.makedirs(os.path.join(self.root, "temp", "fastafiles"), exist_ok=True)
        tempdir = os.path.join(self.root, "temp")

        self.progressbar_counter += 10
        self.submission_progress.setValue(self.progressbar_counter)

        version = str(date.today()).split("-")[0:2]
        version = "".join(version[0] + "." + version[1])

        organism = str(self.organism_input.text()).replace(" ", "_")

        try:
            genus, species, strain = organism.split("_")
        except:
            QtWidgets.QMessageBox.about(self, "Error",
                                        "The organism needs to have 3 fields: <b>genus</b> <b>species</b> <b>strain</b> separated by spaces")
            self.progressbar_counter = 0
            self.submission_progress.setValue(self.progressbar_counter)
            self.statusBar().showMessage("Ready")
            return None

        species = species.lower()

        organism = f"{genus}_{species}_{strain}"

        company = self.company_input.text()

        if " " in company:
            company = company.replace(" ", "_")

        common_name = self.common_input.text()
        if " " in common_name:
            common_name = common_name.replace(" ", "_")

        if not self.filepaths['genome']:
            QtWidgets.QMessageBox.about(self, "Error", "No genome fasta file detected")
            self.progressbar_counter = 0
            self.submission_progress.setValue(self.progressbar_counter)
            self.statusBar().showMessage("Upload genome fasta file")
            return None

        if not self.filepaths['multifasta']:
            QtWidgets.QMessageBox.about(self, "Error", "No multigene fasta file detected")
            self.progressbar_counter = 0
            self.submission_progress.setValue(self.progressbar_counter)
            self.statusBar().showMessage("Upload  file")
            return None

        if not self.filepaths['gff']:
            QtWidgets.QMessageBox.about(self, "Error", "No gff file detected")
            self.submission_progress.setValue(self.progressbar_counter)
            self.statusBar().showMessage("Upload GFF file")
            return None

        db = Database(database=organism)
        db.create_new_database(gff_file=self.filepaths['gff'],
                               genome=self.filepaths['genome'],
                               multifasta=self.filepaths['multifasta'])

        tempfa = pd.read_sql("SELECT * FROM genome",
                             sqlite3.connect(os.path.join(os.getcwd(), "databases", f"{organism}.db")))
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
            'Package': f'BSgenome.{common_name}.{company}.{strain}',
            'Title': f'Genome sequence of {organism}',
            'Description': f'BSgenome genome obtained from the genome sequence of {organism}',
            'Version': version,
            'Author': self.author_input.text(),
            'License': "GPLv2",
            'organism': organism,
            'common_name': common_name,
            'provider': self.author_input.text(),
            'provider_version': f'{strain}.{str(date.today())}',
            'release_date': str(date.today()),
            'release_name': f'{company} {self.author_input.text()} {strain}',
            'organism_biocview': f'{genus}_{species}',
            'BSgenomeObjname': common_name,
            'seqnames': seqnames,
            'seqs_srcdir': os.path.join(tempdir, "fastafiles")

        }

        with open(os.path.join(tempdir, "dcf_file.dcf"), 'w') as dcf:
            for k, v in dcf_file.items():
                dcf.write(f"{k}: {v} \n")

        bsgenome_worker = BSgenome_worker()
        self.threadingPool.start(bsgenome_worker)

        while self.threadingPool.activeThreadCount() == 1:
            QtWidgets.QApplication.processEvents()
            self.submit_newdb_btn.setEnabled(False)
            self.submit_gRNA_predict.setEnabled(False)
            self.processed_submit.setEnabled(False)
            self.statusBar().showMessage("buidling sql database...")
            if self.progressbar_counter < 99:
                self.progressbar_counter += 2
                self.submission_progress.setValue(self.progressbar_counter)
                time.sleep(1)

        if self.threadingPool.waitForDone():
            self.availible_pams.clear()
            self.current_databases.clear()
            self.search_databases.clear()
            self.availible_databases.clear()
            self.availible_bsgenome.clear()

            if ".db" == [items for items in os.listdir(os.path.join(self.root, "databases"))][0]:
                os.remove(os.path.join(self.root, "databases", ".db"))

            self.availible_databases = [str(os.path.splitext(items)[0]).replace("_", " ") for items in
                                        os.listdir(os.path.join(self.root, "databases")) if items != "test.db"]
            self.availible_bsgenome = [str(os.path.splitext(items)[0]) for items in
                                       os.listdir(os.path.join(self.root, "R", "local_packages"))]
            self.current_databases.addItems(self.availible_databases)
            self.search_databases.addItems(self.availible_databases)
            self.progressbar_counter = 100
            self.submission_progress.setValue(self.progressbar_counter)
            self.submit_newdb_btn.setEnabled(True)
            self.submit_gRNA_predict.setEnabled(True)
            self.processed_submit.setEnabled(True)

            self.filepaths = self.filepaths.fromkeys(self.filepaths, "")
            self.gff_uploaddescription.setText("no file loaded")
            self.gff_uploadbtn.setText("upload")
            self.gff_uploadbtn.setStyleSheet("background: light; color: black")
            self.gene_fasta_uploaddescription.setText("no file loaded")
            self.gene_fasta_uploadbtn.setText("upload")
            self.gene_fasta_uploadbtn.setStyleSheet("background: light; color: black")
            self.genome_fasta_uploaddescription.setText("no file loaded")
            self.genome_fasta_uploadbtn.setText("upload")
            self.genome_fasta_uploadbtn.setStyleSheet("background: light; color: black")

            self.statusBar().showMessage("Ready")
            self.progressbar_counter = 0
            self.submission_progress.setValue(self.progressbar_counter)
            if ".db" == [items for items in os.listdir(os.path.join(self.root, "databases"))][0]:
                os.remove(os.path.join(self.root, "databases", ".db"))

    # worker functions
    def _run_gRNA_prediction(self):
        """
        Populate sqlite database with gRNA candidates and offtargets, the database needs to exist
        """
        self.statusBar().showMessage("Creating guide RNA database...")
        self.gRNA_progressbar_counter += 1
        self.gRNA_progress.setValue(self.gRNA_progressbar_counter)

        if not "temp" in os.listdir(self.root):
            os.mkdir(os.path.join(self.root, "temp"))

        tempdir = os.path.join(self.root, "temp")

        if "global_gRNA" not in tempdir:
            os.mkdir(os.path.join(tempdir, "global_gRNA"))

        global_gRNA = os.path.join(tempdir, "global_gRNA")

        chosen_org = str(self.search_databases.currentText())

        chosen_org_modified = chosen_org.replace(" ", "_")

        strain = chosen_org.split(" ")[-1]

        # create fasta file for processing
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

        with open(os.path.join(tempdir, f"{chosen_org_modified}.fasta"), 'w') as input_fasta:
            for header, sequence in zip(fasta_header, fasta_sequence):
                input_fasta.write(">" + header + '\n')
                input_fasta.write(sequence + '\n')

        fasta_header.clear()
        fasta_sequence.clear()

        fasta_header = [header for header in tempgenes['header']]
        fasta_sequence = [sequence for sequence in tempgenes['sequence']]

        with open(os.path.join(tempdir, f"{chosen_org_modified}_genes.fasta"), 'w') as input_fasta:
            for header, sequence in zip(fasta_header, fasta_sequence):
                input_fasta.write(">" + header + '\n')
                input_fasta.write(sequence + '\n')

        fasta_header.clear()
        fasta_sequence.clear()

        is_circular = self.circ_chromosome.text()
        avail_options = ["TRUE", "FALSE"]

        if str(is_circular) not in avail_options:
            circ_msg = QtWidgets.QMessageBox.warning(self, "warning",
                                                     "Cirular has to be TRUE or FALSE in capital letters",
                                                     QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.Ok)

        for objects in self.availible_bsgenome:
            detected_package = objects.split("_")[0]
            detected_package = detected_package.split(".")[-1]
            if strain == detected_package:
                bsgenome_package = str(objects.split("_")[0])

        if self.PAM_scoring_HSU_Zhang.isChecked():
            scoring_method = "Hsu-Zhang"
        else:
            scoring_method = "CFDscore"

        genus, species, strain = str(self.search_databases.currentText()).split(" ")
        species = species.lower()

        config_file = {
            'organism': [f"{genus} {species} {strain}"],
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
            'max_mismatch_gRNA': [self.mismatch.value()],
            'PAM_sequence': [self.PAM_input.text()],
            'PAM_length': [len(str(self.PAM_input.text()))],
            'scoring_method': [scoring_method]
        }

        config_file = pd.DataFrame(config_file)
        config_file.to_csv(os.path.join(tempdir, "config.txt"), index=False, sep="\t")

        findgRNA_worker = FindgRNA_worker()
        database = Database(database=chosen_org_modified)

        self.threadingPool.start(findgRNA_worker)

        while self.threadingPool.activeThreadCount() == 1:
            self.statusBar().showMessage("predicting all guide RNA's...")
            self.submit_gRNA_predict.setEnabled(False)
            self.submit_newdb_btn.setEnabled(False)
            self.processed_submit.setEnabled(False)
            QtWidgets.QApplication.processEvents()
            if self.gRNA_progressbar_counter < 90:
                self.gRNA_progressbar_counter += 1
                self.gRNA_progress.setValue(self.gRNA_progressbar_counter)
                time.sleep(2)

        if self.threadingPool.waitForDone():
            self.statusBar().showMessage("buidling global gRNA database...")
            database.create_gRNA_database(summary=os.path.join(global_gRNA, "Summary.xls"),
                                          offtarget=os.path.join(global_gRNA, "OfftargetAnalysis.txt"),
                                          config_file=os.path.join(tempdir, "config.txt"))

            self.gRNA_progressbar_counter += 9
            self.gRNA_progress.setValue(self.gRNA_progressbar_counter)

            # widgets
            self.current_pams.clear()
            self.metadata_mismatches.clear()

            # lists
            self.availible_pams.clear()
            self.metadata_mismatches_value.clear()

            # update add metadata updates
            sqlrunner_predict = SQL(database=chosen_org_modified)

            self.availible_pams = sqlrunner_predict.list_pams()
            self.metadata_mismatches_value = sqlrunner_predict.list_mismatches(pam=self.availible_pams[0])
            self.current_pams.addItems(self.availible_pams)
            self.metadata_mismatches.addItems(self.metadata_mismatches_value)

            self.gRNA_progressbar_counter += 1
            self.gRNA_progress.setValue(self.gRNA_progressbar_counter)
            self.submit_gRNA_predict.setEnabled(True)
            self.submit_newdb_btn.setEnabled(True)
            self.processed_submit.setEnabled(True)
            self.gRNA_progressbar_counter = 0
            self.gRNA_progress.setValue(self.gRNA_progressbar_counter)
            shutil.rmtree(tempdir)
            self.statusBar().showMessage("Ready")

    def run_gRNA_processing(self):
        """
        Essentially the database gets querried here
        """
        if "temp" not in os.listdir(self.root):
            tempdir = os.path.join(self.root, "temp")
            os.mkdir(os.path.join(tempdir))

        progress_bar_val = 0

        self.processed_import_dialog.setText("Import genelist")
        self.processed_import_dialog.setStyleSheet("background: light; color: black")
        database = str(self.processed_organism.currentText()).replace(" ", "_")
        mismatch = self.processed_maxMismatch.currentText()
        max_grna = self.processed_maxPrimers.value()
        max_primer_len = self.processed_primerlen.value()

        self.gene_mask_dict['genes'] = [items.replace("_", "") if "_" in items else items for items in
                                        self.gene_mask_dict['genes']]
        self.gene_mask_dict['masks'] = [items.replace("_", "") if "_" in items else items for items in
                                        self.gene_mask_dict['masks']]

        # Strand is r for reverse
        worker = CrisprInterference_worker(database=database,
                                           mismatch=mismatch,
                                           strand='r',
                                           max_grna=max_grna,
                                           genes_masks=self.gene_mask_dict,
                                           max_primer_size=max_primer_len)

        self.threadingPool.start(worker)

        while self.threadingPool.activeThreadCount() == 1:
            self.statusBar().showMessage("Gathering guide RNA's...")
            QtWidgets.QApplication.processEvents()
            self.submit_gRNA_predict.setEnabled(False)
            self.submit_newdb_btn.setEnabled(False)
            self.processed_submit.setEnabled(False)

            if progress_bar_val < 90:
                progress_bar_val += 1
                time.sleep(0.8)
                self.processed_progressbar.setValue(progress_bar_val)

        if self.threadingPool.waitForDone():
            self.statusBar().showMessage("Gathering data ...")
            self.candidate_gRNA_df = pd.read_csv(filepath_or_buffer=os.path.join(self.root, "temp", "candidates.txt"),
                                                 sep=",")
            self.backup_gRNA_df = pd.read_csv(filepath_or_buffer=os.path.join(self.root, "temp", "backup.txt"), sep=",")
            self.dropped_gRNA_df = pd.read_csv(filepath_or_buffer=os.path.join(self.root, "temp", "dropped.txt"),
                                               sep=",")
            self.offtarget_df = pd.read_csv(filepath_or_buffer=os.path.join(self.root, "temp", "offtargets.txt"),
                                            sep=",")

            cand_model, backup_model, dropped_model, offtargets_model = map(PandasModel, [self.candidate_gRNA_df,
                                                                                          self.backup_gRNA_df,
                                                                                          self.dropped_gRNA_df,
                                                                                          self.offtarget_df])
            # model_candidates = PandasModel(candidates)

            while progress_bar_val < 100:
                progress_bar_val += 1
                time.sleep(0.01)
                self.statusBar().showMessage("Formatting for display...")
                self.processed_progressbar.setValue(progress_bar_val)

            self.display_candidates.setModel(cand_model)
            self.display_backup.setModel(backup_model)
            self.display_dropped.setModel(dropped_model)
            self.display_offtargets.setModel(offtargets_model)

            self.submit_gRNA_predict.setEnabled(True)
            self.submit_newdb_btn.setEnabled(True)
            self.processed_submit.setEnabled(True)
            self.database_querried = True
            self.statusBar().showMessage("Ready")
            shutil.rmtree(tempdir)
            self.processed_progressbar.setValue(0)

    # internal event modifiers here
    def closeEvent(self, event):
        user_response = QtWidgets.QMessageBox.question(self, 'Message', "Are you sure you would like to quit",
                                                       QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                       QtWidgets.QMessageBox.No)

        files_in_root = os.listdir(os.path.dirname(os.path.abspath(__name__)))
        if "temp" in files_in_root:
            shutil.rmtree(os.path.join(os.path.dirname(os.path.abspath(__name__)), "temp"))

        if "scores.RDS" in files_in_root:
            os.remove(os.path.join(os.path.dirname(os.path.abspath(__name__)), "scores.RDS"))

        if "offTargets.RDS" in files_in_root:
            os.remove(os.path.join(os.path.dirname(os.path.abspath(__name__)), "offTargets.RDS"))


def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    gui = MainWindow()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
