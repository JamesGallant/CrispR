import pandas as pd

from src.py import utilities


class RefineCripri:
    """
    Dataframe is a pandas dataframe
    """

    def __init__(self, grna_dataframe, strand, fasta_dataframe, cas9):
        self.grna_dataframe = grna_dataframe
        self.strand = strand
        self.fasta = fasta_dataframe
        self.utils = utilities.CrispinatorUitls()
        self.cas9 = cas9
        self.ranked_data = self.initial_scoring()

    def initial_scoring(self):
        rank_df = self.utils.rank_pams(dictionary=self.grna_dataframe.to_dict('list'), cas9=self.cas9)
        rank_df = pd.DataFrame(rank_df)

        score = [row['rank'] + 0.5 if row['gRNAsPlusPAM'][0] == "C" or row['gRNAsPlusPAM'][0] == "T" else row['rank']
                 for idx, row in rank_df.iterrows()]

        rank_df['score'] = score

        rank_df.rename(columns={'gRNAsPlusPAM': 'gRNA'}, inplace=True)

        return rank_df.to_dict('list')

    def detected_genes(self):
        detected = self.ranked_data['names'].tolist()
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
            'notes': [],
            'PAM': []
        }

        all_searched_genes = [header for header in self.fasta['header']]
        not_detected = list(set(all_searched_genes) - set(self.detected_genes()))
        for genes in not_detected:
            out['name'].append(genes)
            out['gRNA'].append("NA")
            out['PAM'].append("NA")
            out['notes'].append(f"No gRNA detected based on PAM: {set(self.ranked_data['PAM'])}")

        return out

    def initial_filter(self):
        """
        first pass iteration to filter for possible matches. This function will return three dictionaries:
        canidates which have passed constraints, backup candidates that can be accepted in follow up rounds and those that need to be ignored.
        Contstraints are gRNA has to be on the reverse strand and the 20th position is either a G/A
        """

        ranked_df = pd.DataFrame(self.ranked_data)

        drop_strand_idx = [idx for idx, row in ranked_df.iterrows() if row['names'][-1] != self.strand]
        keep_strand_idx = [idx for idx, row in ranked_df.iterrows() if row['names'][-1] == self.strand]

        temp = ranked_df.loc[keep_strand_idx, :]
        temp.reset_index(inplace=True, drop=True)
        to_candidates = [idx for idx, row in temp.iterrows() if row['score'] <= 1.5]
        to_backup = [idx for idx, row in temp.iterrows() if row['score'] > 1.5]

        dropped = ranked_df.loc[drop_strand_idx, :]
        candidates = temp.loc[to_candidates, :]
        backup = temp.loc[to_backup, :]

        candidates.reset_index(inplace=True, drop=True)
        backup.reset_index(inplace=True, drop=True)
        dropped.reset_index(inplace=True, drop=True)

        drop_notes = ["Primer on the reverse strand"] * len(dropped.index)
        candidates_notes = ["position 20 is C/T" if row['score'] % 1 == 0.5 else "PASS" for _, row in
                            candidates.iterrows()]
        backup_notes = ["position 20 is C/T" if row['score'] % 1 == 0.5 else f"PASS: PAM Rank is {row['score']}" for
                        _, row in backup.iterrows()]

        candidates['notes'] = candidates_notes
        backup['notes'] = backup_notes
        dropped['notes'] = drop_notes

        candidates.drop(['top5OfftargetTotalScore', 'primer_length', 'mismatch_length'], axis=1, inplace=True)
        backup.drop(['top5OfftargetTotalScore', 'primer_length', 'mismatch_length'], axis=1, inplace=True)
        dropped.drop(['top5OfftargetTotalScore', 'primer_length', 'mismatch_length'], axis=1, inplace=True)

        return candidates.to_dict('list'), backup.to_dict('list'), dropped.to_dict('list')

    def correct_gRNA_dictionary(self, candidates, backup, dropped_gRNA):
        """
        returns a mulitple dictionary
        checks if gRNA's are in backup and not featured in the valids. In this case the guide RNA from a specific
        gene needs to be moved for further processing. If it not found in valids then it will be moved to valid
        which disregards previous filtering criteria to reduce gRNA loss.
        """

        backup2 = {
            'names': [],
            'gRNA': [],
            'PAM': [],
            'rank': [],
            'score': [],
            'notes': []
        }

        candidate_genes = [names.split("_")[0] for names in candidates['names']]
        backup_genes = [names.split("_")[0] for names in backup['names']]
        genes_to_move = set(backup_genes) - set(candidate_genes)

        for genes in genes_to_move:
            for idx, items in enumerate(backup['names']):
                if genes in items:
                    candidates['names'].append(backup['names'][idx])
                    candidates['gRNA'].append(backup["gRNA"][idx])
                    candidates['PAM'].append(backup["PAM"][idx])
                    candidates['rank'].append(backup["rank"][idx])
                    candidates['score'].append(backup["score"][idx])
                    candidates['notes'].append(backup["notes"][idx])

        candidate_genes_new = [names for names in candidates['names']]
        for idx, value in enumerate(backup['names']):
            if value not in candidate_genes_new:
                backup2['names'].append(backup['names'][idx])
                backup2['gRNA'].append(backup['gRNA'][idx])
                backup2['PAM'].append(backup['PAM'][idx])
                backup2['rank'].append(backup['rank'][idx])
                backup2['score'].append(backup['score'][idx])
                backup2['notes'].append(backup['notes'][idx])

        return candidates, backup2, dropped_gRNA

    def has_offtarget(self, candidates, backup, dropped_gRNA):
        """
        Takes the three input dictionaries, performs the function and returns them again. This function will check if there are
        off targets present, if not it will be moved to candidates. If off targets are present it will be moved to backup for
        follow up. Executes a correction before returning the dataframe
        """

        ranked_df = pd.DataFrame(self.ranked_data)
        # candidates = pd.DataFrame(candidates)
        for identifiers in candidates["names"]:
            cross_reference = ranked_df.loc[ranked_df['names'] == identifiers].to_dict()
            off_target_presence = str(list(cross_reference['top5OfftargetTotalScore'].values())[0])
            if off_target_presence != "nan":
                backup['names'].append(identifiers)
                backup['gRNA'].append(str(list(cross_reference['gRNA'].values())[0]))
                backup['PAM'].append(str(list(cross_reference['PAM'].values())[0]))
                backup['rank'].append(str(list(cross_reference['rank'].values())[0]))
                # backup['genes'].append(str(list(cross_reference['genes'].values())[0]))
                backup['score'].append(int(list(cross_reference['score'].values())[0]))
                backup['notes'].append("Has off target")

                candidates = pd.DataFrame(candidates)
                candidates = candidates[candidates['names'] != identifiers]
                candidates.reset_index(drop=True, inplace=True)
                candidates = candidates.to_dict('list')

        newscore_cand = [score + 0.4 if "Has off target" in candidates['notes'][idx] else score for idx, score in
                         enumerate(candidates['score'])]

        newscore_back = [score + 0.4 if "Has off target" in backup['notes'][idx] else score for idx, score in
                         enumerate(backup['score'])]

        candidates['score'] = newscore_cand
        backup['score'] = newscore_back

        candidates, backup, dropped_gRNA = self.correct_gRNA_dictionary(candidates=candidates, backup=backup,
                                                                        dropped_gRNA=dropped_gRNA)

        return candidates, backup, dropped_gRNA

    def cripr_interference(self):
        candidates, backup, dropped = self.has_offtarget(
            candidates=self.correct_gRNA_dictionary(candidates=self.initial_filter()[0],
                                                    backup=self.initial_filter()[1],
                                                    dropped_gRNA=self.initial_filter()[2])[0],
            backup=self.correct_gRNA_dictionary(candidates=self.initial_filter()[0],
                                                backup=self.initial_filter()[1],
                                                dropped_gRNA=self.initial_filter()[2])[1],
            dropped_gRNA=self.correct_gRNA_dictionary(candidates=self.initial_filter()[0],
                                                      backup=self.initial_filter()[1],
                                                      dropped_gRNA=self.initial_filter()[2])[2])

        return candidates, backup, dropped
