import unittest
import pandas as pd

from pandas.testing import assert_frame_equal
from src.py.workers import CrisprInterference_worker
from collections import Counter


class CrispriHelpers:
    def __init__(self, gc_content_string: str, cas9: str, max_grna: int):
        self.gc_content_string = gc_content_string.lower()
        self.cas9 = cas9
        self.max_grna = max_grna

    def scan_mismatch_test(self):
        CrisprI = CrisprInterference_worker(database=None, mismatch=None, strand=None, max_grna=self.max_grna,
                                            genes_masks=None, max_primer_size=None, cas9_organism=None)

        candidates = {
            'name': ["gene1_primer1", "gene1_primer2", "gene1_primer3", "gene2_primer1", "gene2_primer2"],
            'gRNA': ["NNNN", "NNNN", "NNNN", "NNNN", "NNNN"],
            'PAM': ["AGAAG", "AGAAT", "AGAAT", "AGAAG", "AGAAC"],
            'rank': [1, 2, 2, 1, 5],
            'notes': ["PASS", "PASS", "Position 20 is", "Position 20 is", "Position 20 is"],
            'genes': ["Gene1", "Gene1", "Gene1", "Gene2", "Gene2"]
        }

        backup = {
            'name': ["gene1_primer5", "gene1_primer6", "gene2_primer3"],
            'gRNA': ["NNNN", "NNNN", "NNNN"],
            'PAM': ["AGAAT", "AGAAG", "AGAAG"],
            'rank': [2, 1, 1],
            'notes': ["PASS", "PASS", "PASS"],
            'genes': ["Gene1", "Gene1", "Gene2"]
        }

        candidates, backup = CrisprI.scan_maxmismatches(candidates=candidates,
                                                        backup=backup)

        return pd.DataFrame(candidates), pd.DataFrame(backup)

    def scan_mismatch_result(self):
        def _one():
            candidates = {
                'name': ["gene1_primer1", "gene2_primer1"],
                'gRNA': ["NNNN", "NNNN"],
                'PAM': ["AGAAG", "AGAAG"],
                'rank': [1, 1],
                'notes': ["PASS", "Position 20 is"],
                'genes': ["Gene1", "Gene2"]
            }

            backup = {
                'name': ["gene1_primer5", "gene1_primer6", "gene1_primer2",
                         "gene1_primer3", "gene2_primer3", "gene2_primer2"],
                'gRNA': ["NNNN", "NNNN", "NNNN", "NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAT", "AGAAG", "AGAAT", "AGAAT", "AGAAG", "AGAAC"],
                'rank': [2, 1, 2, 2, 1, 5],
                'notes': ["PASS", "PASS", "PASS", "Position 20 is" "PASS", "Position 20 is"],
                'genes': ["Gene1", "Gene1", "Gene1", "Gene1", "Gene2", "Gene2"]
            }
            return candidates, backup

        def _two():
            candidates = {
                'name': ["gene1_primer1", "gene1_primer2", "gene2_primer1", "gene2_primer2"],
                'gRNA': ["NNNN", "NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAG", "AGAAT", "AGAAG", "AGAAC"],
                'rank': [1, 2, 1, 5],
                'notes': ["PASS", "PASS", "Position 20 is", "Position 20 is"],
                'genes': ["Gene1", "Gene1", "Gene2", "Gene2"]
            }

            backup = {
                'name': ["gene1_primer5", "gene1_primer6", "gene1_primer3", "gene2_primer3"],
                'gRNA': ["NNNN", "NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAT", "AGAAG", "AGAAT", "AGAAG"],
                'rank': [2, 1, 2, 1],
                'notes': ["PASS", "PASS", "Position 20 is", "PASS"],
                'genes': ["Gene1", "Gene1", "Gene1", "Gene2"]
            }
            return candidates, backup

        def _three():
            candidates = {
                'name': ["gene1_primer1", "gene1_primer2", "gene1_primer3", "gene2_primer1", "gene2_primer2"],
                'gRNA': ["NNNN", "NNNN", "NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAG", "AGAAT", "AGAAT", "AGAAG", "AGAAC"],
                'rank': [1, 2, 2, 1, 5],
                'notes': ["PASS", "PASS", "Position 20 is", "Position 20 is", "Position 20 is"],
                'genes': ["Gene1", "Gene1", "Gene1", "Gene2", "Gene2"]
            }

            backup = {
                'name': ["gene1_primer5", "gene1_primer6", "gene2_primer3"],
                'gRNA': ["NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAT", "AGAAG", "AGAAG"],
                'rank': [2, 1, 1],
                'notes': ["PASS", "PASS", "PASS"],
                'genes': ["Gene1", "Gene1", "Gene2"]
            }

            return candidates, backup

        switch = {
            1: _one(),
            2: _two(),
            3: _three()
        }

        candidates, backup = switch.get(self.max_grna, 3)

        return pd.DataFrame(candidates), pd.DataFrame(backup)

    def gc_content_test(self):
        gc_content = self.gc_content_string
        CriprI = CrisprInterference_worker(database=None, mismatch=None, strand=None, max_grna=None, genes_masks=None,
                                           max_primer_size=None, cas9_organism=None)
        dataframe = pd.DataFrame({'gRNA': [gc_content]})
        out = CriprI.calculate_gc_content(dataframe)
        return str(out['gc_content'][0])

    def gc_content_result(self):
        gc_content = self.gc_content_string
        basecount = Counter(gc_content)
        guanine = basecount['g']
        cytosine = basecount['c']
        gc_content = str(round(((guanine + cytosine) / len(gc_content)) * 100, 3))
        return gc_content

    def negate_pam_mismatch_test(self):
        CriprII = CrisprInterference_worker(database=None, mismatch=None, strand=None, max_grna=None, genes_masks=None,
                                            max_primer_size=None, cas9_organism=self.cas9)

        dataframe = {
            'name': ['grna_primer1', 'grna_primer2', 'grna_primer3', 'grna_primer4'],
            'gRNA': ["NNNNNNNAGAAG", "NNCNNNNAGAAG", "NNNNNNNAGAAT", "NNNNNNNNNNNN"],
            'PAM': ['AGAAG', 'AGAAG', 'AGAAG', 'AGAAG'],
            'notes': ["PASS", "Has off target", "Has off target", "Has off target"]
        }
        offtarget = {
            'name': ['grna_primer2', 'grna_primer3', 'grna_primer4'],
            'OffTargetSequence': ['NNCNNNNAGAAG', 'NNCNNNNAGAAT', 'NNNNNNNNNNNN']
        }

        out = CriprII.negate_pam_mismatch(grna_dataframe=dataframe,
                                          offtarget_dataframe=offtarget)

        return pd.DataFrame(out)

    def negate_pam_mismatch_result(self):
        dataframe = {
            'name': ['grna_primer1', 'grna_primer2', 'grna_primer3', 'grna_primer4'],
            'gRNA': ["NNNNNNNAGAAG", "NNCNNNNAGAAG", "NNNNNNNAGAAT", "NNNNNNNNNNNN"],
            'PAM': ['AGAAG', 'AGAAG', 'AGAAG', 'AGAAG'],
            'notes': ["PASS", "Has off target", "Has off target", "PASS"]
        }
        return pd.DataFrame(dataframe)

    def calculate_primerlenght_test(self):
        CriprI = CrisprInterference_worker(database=None, mismatch=None, strand=None, max_grna=None, genes_masks=None,
                                           max_primer_size=None, cas9_organism=None)

        out = CriprI.calculate_primer_len(dataframe={'gRNA': ["N", "NN", "NNN"]})
        return list(out['primer_length'])

    def calculate_primerlenght_result(self):
        return [1, 2, 3]

    def force_max_grna_test(self):
        CriprI = CrisprInterference_worker(database=None, mismatch=None, strand=None, max_grna=None, genes_masks=None,
                                           max_primer_size=None, cas9_organism=None)

        candidates = {
            'name': ["gene1_primer1", "gene1_primer2", "gene2_primer1"],
            'gRNA': ["NNNN", "NNNN", "NNNN"],
            'PAM': ["AGAAG", "AGAAG", "AGAAG"],
            'rank': [1, 1, 1],
            'notes': ["PASS", "PASS", "PASS"],
            'genes': ["Gene1", "Gene1", "Gene2"]
        }

        backup = {
            'name': ["gene1_primer3", "gene1_primer4", "gene2_primer2"],
            'gRNA': ["NNNN", "NNNN", "NNNN"],
            'PAM': ["AGAAT", "AGAAG", "AGAAG"],
            'rank': [2, 1, 1],
            'notes': ["PASS", "PASS", "PASS"],
            'genes': ["Gene1", "Gene1", "Gene2"]
        }

        candidates, backup = CriprI.force_max_grna_in_candidates(candidates=candidates,
                                                                 backup=backup,
                                                                 max_grna=self.max_grna)
        return pd.DataFrame(candidates), pd.DataFrame(backup)

    def force_max_grna_result(self):
        def _one():
            candidates = {
                'name': ["gene1_primer1", "gene1_primer2", "gene2_primer1"],
                'gRNA': ["NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAG", "AGAAG", "AGAAG"],
                'rank': [1, 1, 1],
                'notes': ["PASS", "PASS", "PASS"],
                'genes': ["Gene1", "Gene1", "Gene2"]
            }

            backup = {
                'name': ["gene1_primer3", "gene1_primer4", "gene2_primer1"],
                'gRNA': ["NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAT", "AGAAG", "AGAAG"],
                'rank': [2, 1, 1],
                'notes': ["PASS", "PASS", "PASS"],
                'genes': ["Gene1", "Gene1", "Gene2"]
            }
            return candidates, backup

        def _two():
            candidates = {
                'name': ["gene1_primer1", "gene1_primer2", "gene2_primer1", "gene2_primer2"],
                'gRNA': ["NNNN", "NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAG", "AGAAG", "AGAAG", "AGAAG"],
                'rank': [1, 1, 1, 1],
                'notes': ["PASS", "PASS", "PASS", "PASS"],
                'genes': ["Gene1", "Gene1", "Gene2", "Gene2"]
            }

            backup = {
                'name': ["gene1_primer3", "gene1_primer4"],
                'gRNA': ["NNNN", "NNNN"],
                'PAM': ["AGAAT", "AGAAG"],
                'rank': [2, 1],
                'notes': ["PASS", "PASS"],
                'genes': ["Gene1", "Gene1"]
            }
            return candidates, backup

        def _three():
            candidates = {
                'name': ["gene1_primer1", "gene1_primer2", "gene1_primer4", "gene2_primer1", "gene2_primer2"],
                'gRNA': ["NNNN", "NNNN", "NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAG", "AGAAG", "AGAAG", "AGAAG", "AGAAG"],
                'rank': [1, 1, 1, 1, 1],
                'notes': ["PASS", "PASS", "PASS", "PASS", "PASS"],
                'genes': ["Gene1", "Gene1", "Gene1", "Gene2", "Gene2"]
            }

            backup = {
                'name': ["gene1_primer3"],
                'gRNA': ["NNNN"],
                'PAM': ["AGAAT"],
                'rank': [2],
                'notes': ["PASS"],
                'genes': ["Gene1"]
            }
            return candidates, backup

        def _four():
            candidates = {
                'name': ["gene1_primer1", "gene1_primer2", "gene1_primer4", "gene1_primer3", "gene2_primer1",
                         "gene2_primer2"],
                'gRNA': ["NNNN", "NNNN", "NNNN", "NNNN", "NNNN", "NNNN"],
                'PAM': ["AGAAG", "AGAAG", "AGAAG", "AGAAT", "AGAAG", "AGAAG"],
                'rank': [1, 1, 1, 2, 1, 1],
                'notes': ["PASS", "PASS", "PASS", "PASS", "PASS", "PASS"],
                'genes': ["Gene1", "Gene1", "Gene1", "Gene1", "Gene2", "Gene2"]
            }

            backup = {
                'name': [],
                'gRNA': [],
                'PAM': [],
                'rank': [],
                'notes': [],
                'genes': []
            }
            return candidates, backup

        switch = {
            1: _one(),
            2: _two(),
            3: _three(),
            4: _four()
        }

        candidates, backup = switch.get(self.max_grna, 4)

        return pd.DataFrame(candidates), pd.DataFrame(backup)


CrispriRunner = CrispriHelpers(gc_content_string="AGTC", cas9="Streptococcus thermophilus", max_grna=3)


class CrisprInterferenceWorkerTestCase(unittest.TestCase):
    def test_max_mismatch_candidates(self):
        try:
            assert_frame_equal(CrispriRunner.scan_mismatch_result()[0],
                               CrispriRunner.scan_mismatch_test()[0])
        except AssertionError as e:
            raise self.failureException(e)

    def test_max_mismatch_backup(self):
        try:
            assert_frame_equal(CrispriRunner.scan_mismatch_result()[1],
                               CrispriRunner.scan_mismatch_test()[1])
        except AssertionError as e:
            raise self.failureException(e)

    def test_negate_pam_mismatch(self):
        try:
            assert_frame_equal(CrispriRunner.negate_pam_mismatch_result(),
                               CrispriRunner.negate_pam_mismatch_test())
        except AssertionError as e:
            raise self.failureException(e)

    def test_gc_content(self):
        self.assertEqual(CrispriRunner.gc_content_result(),
                         CrispriRunner.gc_content_test())

    def test_primerlength(self):
        self.assertEqual(CrispriRunner.calculate_primerlenght_result(),
                         CrispriRunner.calculate_primerlenght_test())

    def test_force_max_gRNA_candidates(self):
        try:
            assert_frame_equal(CrispriRunner.force_max_grna_result()[0],
                               CrispriRunner.force_max_grna_test()[0])
        except AssertionError as e:
            raise self.failureException(e)

    def test_force_max_gRNA_backup(self):
        try:
            assert_frame_equal(CrispriRunner.force_max_grna_result()[1],
                               CrispriRunner.force_max_grna_test()[1])
        except AssertionError as e:
            raise self.failureException(e)

if __name__ == '__main__':
    unittest.main()
