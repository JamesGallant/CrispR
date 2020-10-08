import os
import pandas as pd
from collections import Counter

from PyQt5 import QtCore
from src.py.build_reference import BSgenome
from src.py.predict_gRNA import Get_gRNA
from src.py.annotate_off_target import Annotater
from src.py.cripr import RefineCripri
from src.py.database_tools import SQL
from src.py.misc_functions import *

from src.py import utilities


# from src.py.misc_functions import possible_pams_ranked

class BSgenome_worker(QtCore.QRunnable):
    def __init__(self):
        super().__init__()
        self.bs = BSgenome()
        self.root = os.path.dirname(os.path.abspath(__name__))

    @QtCore.pyqtSlot()
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

    @QtCore.pyqtSlot()
    def run(self):
        """
		Creates a gRNA dataset from the CrisprSeek library in R, then adds it to the organism DB
		"""
        self.gRNA.create_gRNA_library(config_file=self.config_file,
                                      out_dir=self.global_gRNA)

        annotater = Annotater(config_file=pd.read_csv(self.config_file, sep="\t"),
                              offtarget_file=pd.read_csv(self.offtarget_file, index_col=False, sep="\t"))

        offtarget_file = annotater.annotate()

        offtarget_file.to_csv(path_or_buf=os.path.join(self.global_gRNA, "OfftargetAnalysis.txt"), index=False,
                              header=True, sep="\t")


class CrisprInterference_worker(QtCore.QRunnable):
    def __init__(self, database, mismatch, strand, max_grna, genes_masks, max_primer_size, cas9_organism):
        super().__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))
        self.tempdir = os.path.join(self.root, "temp")
        self.database = database
        self.mismatch = mismatch
        self.strand = strand
        self.max_grna = max_grna
        self.gene_mask_dict = genes_masks
        self.max_primer_size = max_primer_size
        self.cas9_organism = cas9_organism
        self.utils = utilities.CrispinatorUitls(cas9=self.cas9_organism)

    def get_targeted_data(self, dataframe, gene_mask_dict):
        """:parameter gene_mask_dict dictionary with genes on header 'genes'
        :parameter dataframe takes pandas dataframe
        :returns pandas dataframe for downstream processing
        """
        genes = [genes.split("_")[0] for genes in dataframe['names']]
        dataframe['genes'] = genes
        query = gene_mask_dict['genes']
        out = pd.DataFrame()
        for items in query:
            if items in genes:
                grad_idx = [idx for idx, val in dataframe.iterrows() if items in val['genes']]
                out = out.append(dataframe.loc[grad_idx, :], ignore_index=True)

        out.drop(['genes'], axis=1, inplace=True)
        return out

    def scan_maxmismatches(self, candidates, backup):
        """
        :parameter candidates dict
        :parameter backup dict
        :returns dict

		Function takes dictionaries and checks if the number of guide RNAs if more than the max it removes some candidates
		Needs to consider A/G and PAMS preferentially Does not consider PAMS atm
		"""

        candidates, backup = map(pd.DataFrame, [candidates, backup])
        number_of_genes = [genes.split("_")[0] for genes in candidates['names']]
        out = Counter(number_of_genes)

        for genes, gene_counts in out.items():
            if gene_counts > self.max_grna:
                get_gene_idx = [idx for idx, value in enumerate(candidates['names']) if genes in value]
                target_rows = candidates.loc[get_gene_idx, :]
                target_rows.reset_index(inplace=True, drop=True)
                tuple_rank_list = [(idx, int(float(item))) for idx, item in enumerate(target_rows['score'])]
                tuple_rank_list.sort(key=lambda x: x[1])
                target_idx_to_pop = [t[0] for t in tuple_rank_list[-(len(get_gene_idx) - self.max_grna):]]
                target_rows.drop(columns=["score"], inplace=True)
                target_rows.reset_index(inplace=True, drop=True)
                move_to_backup = target_rows.loc[target_idx_to_pop, :]

                primers = [target_rows['names'][idx] for idx in target_idx_to_pop]
                candidates = candidates[~candidates['names'].isin(primers)]
                # candidates.drop(target_idx_to_pop, inplace=True)
                candidates.sort_values(by=['genes'], inplace=True)
                candidates.reset_index(inplace=True, drop=True)

                backup = backup.append(move_to_backup, ignore_index=True)
                backup.sort_values(by=['genes'], inplace=True)
                backup.reset_index(inplace=True, drop=True)

        return candidates.to_dict('list'), backup.to_dict('list')

    def grab_targets(self, candidates, backup, dropped, gene_mask_dict):
        """
		gets the genes correpsonding to the targets as given by the user returns the altered dictionaries
		"""
        user_candidates = {
            'names': [],
            'gRNA': [],
            'PAM': [],
            'genes': [],
            'rank': [],
            'score': [],
            'notes': []
        }

        user_backup = {
            'names': [],
            'gRNA': [],
            'PAM': [],
            'genes': [],
            'rank': [],
            'score': [],
            'notes': []
        }

        user_dropped = {
            'names': [],
            'gRNA': [],
            'PAM': [],
            'genes': [],
            'rank': [],
            'score': [],
            'notes': []
        }

        for genes in gene_mask_dict['genes']:
            get_candidate_idx = [idx for idx, value in enumerate(candidates['names']) if genes in value]
            get_backup_idx = [idx for idx, value in enumerate(backup['names']) if genes in value]
            get_dropped_idx = [idx for idx, value in enumerate(dropped['names']) if genes in value]

            for items in get_candidate_idx:
                user_candidates['names'].append(candidates['names'][items])
                user_candidates['gRNA'].append(candidates['gRNA'][items])
                user_candidates['PAM'].append(candidates['PAM'][items])
                user_candidates['genes'].append(candidates['genes'][items])
                user_candidates['rank'].append(candidates['rank'][items])
                user_candidates['score'].append(candidates['score'][items])
                user_candidates['notes'].append(candidates['notes'][items])

            for items in get_backup_idx:
                user_backup['names'].append(backup['names'][items])
                user_backup['gRNA'].append(backup['gRNA'][items])
                user_backup['PAM'].append(backup['PAM'][items])
                user_backup['genes'].append(backup['genes'][items])
                user_backup['rank'].append(backup['rank'][items])
                user_backup['score'].append(backup['score'][items])
                user_backup['notes'].append(backup['notes'][items])

            for items in get_dropped_idx:
                user_dropped['names'].append(dropped['names'][items])
                user_dropped['gRNA'].append(dropped['gRNA'][items])
                user_dropped['PAM'].append(dropped['PAM'][items])
                user_dropped['genes'].append(dropped['genes'][items])
                user_dropped['rank'].append(dropped['rank'][items])
                user_dropped['score'].append(dropped['score'][items])
                user_dropped['notes'].append(dropped['notes'][items])

        return user_candidates, user_backup, user_dropped

    def grab_offtargets(self, query, offtargets, offtarget_ids):
        """
		query corresponds to the candidates/backup dicts and off targets is the global offtarget files. Will return a 
		pandas dataframe
		"""
        out = pd.DataFrame()
        offtargets = pd.DataFrame(offtargets)

        for items in query['names']:
            if items in offtarget_ids:
                out = out.append(offtargets[offtargets['name'].str.match(items)], ignore_index=True)

        return out.to_dict('list')

    def list_comparison(self, list1, list2):
        """
		check if list1 has a element present in list 2

		"""
        return bool(set(list1) & set(list2))

    def negate_pam_mismatch(self, grna_dataframe, offtarget_dataframe, target_ids):
        """:arg grna_dataframe: the dataframe containing target gRNA primers
           :arg offtarget_dataframe: the off targets dataframe
           If there are mismatches in the pam, this may negate the pam completely. If the PAM is modified and does
           not match the known pams then  we can ignore that off target and delete it from the list
        """
        grna_dataframe = pd.DataFrame(grna_dataframe)
        offtarget_dataframe = pd.DataFrame(offtarget_dataframe)
        grna_with_offtargets = [idx for idx, val in enumerate(grna_dataframe['names']) if val in target_ids]
        if not grna_with_offtargets:
            return grna_dataframe.to_dict()
        else:
            for index_location in grna_with_offtargets:
                grna_name = grna_dataframe['names'][index_location]
                grna_pam = grna_dataframe['PAM'][index_location]
                offtarget_row = offtarget_dataframe.loc[offtarget_dataframe['name'] == grna_name].index[0]
                offtarget_grna = offtarget_dataframe['OffTargetSequence'][offtarget_row]
                offtarget_pam = offtarget_grna[-len(grna_pam):]

                if grna_pam != offtarget_pam:
                    availible_pams = possible_pams_ranked(cas9=self.cas9_organism).keys()
                    print(f"oftarget: {offtarget_pam}, avail: {availible_pams}")
                    if offtarget_pam not in availible_pams:
                        grna_dataframe['notes'][index_location] = "PASS: offtarget pam is invalid"

            return grna_dataframe.to_dict('list')

    def move_grna_by_offtargets(self, grna_dataframe, dropped_dataframe, offtarget_dataframe, masks):
        """
		This function will do the masking, if a grna is present in grna_dataframe['name'] that is also present in gene_mask_dict['masks']
		and that can be found in the offtarget database based on the function grab_offtargets it will be moved to dropped
		Returns:
		new dropped dataframe, new grna dataframe
		use for both canididates and backups
		"""
        temp = {
            'names': [],
            'gRNA': [],
            'PAM': [],
            'rank': [],
            'score': [],
            'genes': [],
            'notes': []
        }

        # gets the indices where offtarget genes are located in masking list from user, map to grna names to cross reference later and dedup
        drop_offtarget_idx = [idx for idx, value in enumerate(offtarget_dataframe['annotation']) if value in masks]
        grna_to_move = [offtarget_dataframe['name'][i] for i in drop_offtarget_idx]
        grna_to_move = list(set(grna_to_move))

        # map grna to grna dataframe and find incices where they coincide
        grna_idx = [idx for idx, value in enumerate(grna_dataframe['names']) if value in grna_to_move]

        # append the information to dropped data
        for grna_items in grna_idx:
            dropped_dataframe['names'].append(grna_dataframe['names'][grna_items])
            dropped_dataframe['gRNA'].append(grna_dataframe['gRNA'][grna_items])
            dropped_dataframe['rank'].append(grna_dataframe['rank'][grna_items])
            dropped_dataframe['score'].append(grna_dataframe['score'][grna_items])
            dropped_dataframe['genes'].append(grna_dataframe['genes'][grna_items])
            dropped_dataframe['PAM'].append(grna_dataframe['PAM'][grna_items])
            dropped_dataframe['notes'].append("Has off target that was blocked")

        # checks if its neccesary to remake the candidates and does so
        if len(grna_dataframe['names']) > len(grna_idx):
            remake_idx = [idx for idx, value in enumerate(grna_dataframe['names']) if idx not in grna_idx]
            for remake_index in remake_idx:
                temp['names'].append(grna_dataframe['names'][remake_index])
                temp['gRNA'].append(grna_dataframe['gRNA'][remake_index])
                temp['rank'].append(grna_dataframe['rank'][remake_index])
                temp['score'].append(grna_dataframe['score'][remake_index])
                temp['genes'].append(grna_dataframe['genes'][remake_index])
                temp['PAM'].append(grna_dataframe['PAM'][remake_index])
                temp['notes'].append(grna_dataframe['notes'][remake_index])
        else:
            temp = dict.fromkeys(grna_dataframe, [])

        return temp, dropped_dataframe

    def force_max_grna_in_candidates(self, candidates, backup, max_grna):
        """
        :type: Dictionary
        :returns: Dictionary
		Assumes that all backup genes are valid and can be used, this must run AFTER masking algorithms
		This algorithm does not discriminate between offtargets and position 20 grna's in backup
		Takes in a max grna from user as well as candidates and backup in either dict of pandas dataframe form
		"""

        # This algorithm does not discriminate between offtargets and position 20 grna's in backup
        candidates, backup = map(pd.DataFrame, [candidates, backup])
        candidates_genes_count, backup_genes_count = map(Counter, [candidates['genes'],
                                                                   backup['genes']])
        backup_drop_idx = []
        for genes, counts in candidates_genes_count.items():
            if counts < max_grna:
                genes_in_backup = backup_genes_count.get(genes, 0)
                # only runs if there are genes to move
                if genes_in_backup > 0:
                    max_number_to_move = max_grna - counts

                    moving_idx = [idx for idx, value in enumerate(backup['genes']) if value in genes]

                    target_rows = backup.iloc[moving_idx, :]
                    target_rows.reset_index(inplace=True)
                    tuple_rank_list = [(idx, int(item)) for idx, item in enumerate(target_rows['rank'])]
                    tuple_rank_list.sort(key=lambda x: x[1])

                    if len(moving_idx) >= max_number_to_move:
                        tuple_rank_list = tuple_rank_list[:max_number_to_move]
                    else:
                        tuple_rank_list = tuple_rank_list[:counts]

                    backup_index_to_move = [t[0] for t in tuple_rank_list]
                    target_rows = target_rows.iloc[backup_index_to_move, :]
                    drop_idx = [idx for idx in target_rows['index']]

                    target_rows.drop(columns=['index'], inplace=True)

                    candidates = candidates.append(target_rows, ignore_index=True)
                    candidates.sort_values(by=['genes'], inplace=True)
                    candidates.reset_index(inplace=True, drop=True)

                    backup.drop(drop_idx, inplace=True)

        backup.reset_index(drop=True, inplace=True)
        return candidates.to_dict('list'), backup.to_dict('list')

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
                primer_wo_pam_stop = primer_wo_pam_start + 20  # 20 is the primer length, fixed value

                if primer_wo_pam_stop + max_primer_size < len(gene_sequence):
                    counter = 0
                    while counter <= max_primer_size:
                        target_base_location = primer_wo_pam_stop + counter
                        target_base = sequence_swapped[target_base_location]
                        if target_base == "a" or target_base == "g":
                            grna_out = sequence_swapped[loc_in_gene:target_base_location+1]
                            grna_out = grna_out[::-1]
                            dataframe['gRNA'][idx] = grna_out.upper()

                            if "position 20 is" in dataframe['notes'][idx]:
                                dataframe['notes'][idx] = "PASS"

                            break
                        else:
                            counter += 1
                else:
                    counter = 1
                    if primer_wo_pam_stop == len(gene_sequence):
                        grna_out = sequence_swapped[loc_in_gene:len(gene_sequence)]
                        grna_out = grna_out[::-1]
                        dataframe['gRNA'][idx] = grna_out.upper()
                    else:

                        while counter < len(gene_sequence):
                            target_base_location = primer_wo_pam_stop + counter
                            target_base = sequence_swapped[target_base_location]
                            if target_base == "a" or target_base == "g":
                                grna_out = sequence_swapped[loc_in_gene:target_base_location+1]
                                grna_out = grna_out[::-1]
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
            gc_content = str(round(((guanine + cytosine) / len(seqs)) * 100, 3))
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

    @QtCore.pyqtSlot()
    def run(self):
        sqlrunner = SQL(database=os.path.join(self.root, "databases", self.database))
        gRNA_db = sqlrunner.get_global_gRNA(mismatch=str(self.mismatch))

        if bool(self.gene_mask_dict['genes']):
            query_data = self.get_targeted_data(dataframe=gRNA_db, gene_mask_dict=self.gene_mask_dict)
        else:
            query_data = gRNA_db

        multifasta = sqlrunner.get_gene_multifasta()

        gRNA_runner = RefineCripri(grna_dataframe=query_data,
                                   strand=self.strand,
                                   fasta_dataframe=multifasta, cas9=self.cas9_organism)

        candidates, backup, dropped = gRNA_runner.cripr_interference()

        candidates, backup, dropped = map(self.utils.annotate_dataframe,
                                          [candidates, backup, dropped])

        candidates, backup = self.scan_maxmismatches(candidates=candidates, backup=backup)

        offtargets = sqlrunner.get_offtargets_by_mismatch(mismatch=self.mismatch)
        offtargets.dropna(subset=['annotation'], inplace=True)
        offtargets = offtargets[offtargets['strand'] != '+']
        offtargets['annotation'] = offtargets['annotation'].apply(
            lambda x: x.replace("_", "") if isinstance(x, str) else x)

        offtargets = offtargets.query("gene != annotation")
        offtargets.reset_index(drop=True, inplace=True)

        offtarget_ids = list(set(offtargets['name']))
        candidates_has_offtargets = self.list_comparison(list1=candidates['names'], list2=offtarget_ids)
        backup_has_offtargets = self.list_comparison(list1=backup['names'], list2=offtarget_ids)

        if candidates_has_offtargets:
            candidate_off_ids = list(set(candidates['names']) & set(offtarget_ids))
            candidates_offtargets = self.grab_offtargets(query=candidates,
                                                         offtargets=offtargets,
                                                         offtarget_ids=offtarget_ids)

            candidates = self.negate_pam_mismatch(grna_dataframe=candidates,
                                                  offtarget_dataframe=candidates_offtargets,
                                                  target_ids=candidate_off_ids)

            candidates, dropped = self.move_grna_by_offtargets(grna_dataframe=candidates,
                                                               dropped_dataframe=dropped,
                                                               offtarget_dataframe=candidates_offtargets,
                                                               masks=self.gene_mask_dict['masks'])

            candidates_offtargets = pd.DataFrame(candidates_offtargets)
        else:
            candidates_offtargets = dict.fromkeys(offtargets, [])
            candidates_offtargets = pd.DataFrame(candidates_offtargets)

        if backup_has_offtargets:
            backup_off_ids = list(set(backup['names']) & set(offtarget_ids))
            backup_offtargets = self.grab_offtargets(query=backup, offtargets=offtargets, offtarget_ids=offtarget_ids)

            backup = self.negate_pam_mismatch(grna_dataframe=backup,
                                              offtarget_dataframe=backup_offtargets,
                                              target_ids=backup_off_ids)

            backup, dropped = self.move_grna_by_offtargets(grna_dataframe=backup,
                                                           dropped_dataframe=dropped,
                                                           offtarget_dataframe=backup_offtargets,
                                                           masks=self.gene_mask_dict['masks'])

            backup_offtargets = pd.DataFrame(backup_offtargets)

        else:
            backup_offtargets = dict.fromkeys(offtargets, [])
            backup_offtargets = pd.DataFrame(backup_offtargets)

        ## add ranking to pam, move between dataframes if ranking is fucked

        candidates, backup = self.force_max_grna_in_candidates(candidates=candidates, backup=backup,
                                                               max_grna=self.max_grna)

        candidates = self.force_ag_base(dataframe=candidates, max_primer_size=self.max_primer_size)

        backup = self.force_ag_base(dataframe=backup, max_primer_size=self.max_primer_size)

        candidates, backup, dropped = map(self.calculate_primer_len, [candidates, backup, dropped])

        candidates, backup, dropped = map(self.calculate_gc_content, [candidates, backup, dropped])

        candidates, backup, dropped = map(pd.DataFrame,
                                          [candidates, backup, dropped])

        offtarget_empty = [candidates_offtargets.empty, backup_offtargets.empty]

        final_offtargets = pd.DataFrame()
        if not offtarget_empty[0] and not offtarget_empty[1]:
            final_offtargets = candidates_offtargets
            final_offtargets['from'] = "candidates"
            backup_offtargets['from'] = "backup"
            final_offtargets.append(backup_offtargets, ignore_index=True)

        else:
            if not offtarget_empty[0]:
                final_offtargets = candidates_offtargets
                final_offtargets['from'] = "candidates"

            if not offtarget_empty[1]:
                final_offtargets = backup_offtargets
                final_offtargets['from'] = "backup"

        if final_offtargets.empty:
            final_offtargets = pd.DataFrame(columns=offtargets.columns)

        candidates.to_csv(os.path.join(self.root, "temp", "candidates.txt"), header=True, index=False, sep=",")
        backup.to_csv(os.path.join(self.root, "temp", "backup.txt"), header=True, index=False, sep=",")
        dropped.to_csv(os.path.join(self.root, "temp", "dropped.txt"), header=True, index=False, sep=",")
        final_offtargets.to_csv(os.path.join(self.root, "temp", "offtargets.txt"), header=True, index=False, sep=",")


class CustomSQL_worker(QtCore.QRunnable):
    def __init__(self, database, sql_query):
        super().__init__()
        self.root = os.path.dirname(os.path.abspath(__name__))

        if not "temp" in os.listdir(self.root):
            os.mkdir(os.path.join(self.root, "temp"))

        self.database = database
        self.sql_query = sql_query

    @QtCore.pyqtSlot()
    def run(self):
        sqlrunner = SQL(database=self.database)
        out = sqlrunner.custom_sql(statement=self.sql_query)
        out.to_csv(os.path.join(self.root, "temp", "query.txt"), header=True, index=False, sep=",")
