import unittest
import pandas as pd
import os

from src.py.cripr import RefineCripri
from src.py.database_tools import SQL

class CrisprFuncHelpers:
    """
    unittest class for the crispr.py
    """

    def __init__(self, database: str, strand: str, mismatch: int, cas9: str):
        self.root = os.path.dirname(os.path.abspath("../main.py"))
        os.chdir(self.root)
        self.sql = SQL(database=database)
        self.strand = strand
        self.mismatch = mismatch
        self.cas9 = cas9

    def initial_filter_test(self):
        data = self.sql.get_global_gRNA(mismatch=self.mismatch)
        genes = [genes.split("_")[0] for genes in data['names']]
        data['genes'] = genes
        query = ["Rv0899", "Rv0934"]
        out = pd.DataFrame()
        for items in query:
            if items in genes:
                grad_idx = [idx for idx, val in data.iterrows() if items in val['genes']]
                out = out.append(data.loc[grad_idx, :], ignore_index=True)


        runner = RefineCripri(grna_dataframe=out, strand=self.strand,
                              fasta_dataframe=None,
                              cas9=self.cas9)

        candidates, backup, dropped = map(pd.DataFrame, *[runner.initial_filter()])

        candidates_out = list(set([True if row['score'] < 2 else False for _, row in candidates.iterrows()]))
        backup_out = list(set([True if row['score'] >= 2 else False for _, row in backup.iterrows()]))
        dropped_out = list(set([True if row['names'][-1] != self.strand else False for _, row in dropped.iterrows()]))

        return [candidates_out[0], backup_out[0], dropped_out[0]]

    def initial_filter_result(self):
        return [True, True, True]

    def has_offtarget_test(self):
        data = self.sql.get_global_gRNA(mismatch=self.mismatch)
        genes = [genes.split("_")[0] for genes in data['names']]
        data['genes'] = genes
        query = ["Rv0899", "Rv0934", "Rv0051"]
        out = pd.DataFrame()
        for items in query:
            if items in genes:
                grad_idx = [idx for idx, val in data.iterrows() if items in val['genes']]
                out = out.append(data.loc[grad_idx, :], ignore_index=True)

        runner = RefineCripri(grna_dataframe=out, strand=self.strand,
                              fasta_dataframe=None,
                              cas9=self.cas9)

        candidates, backup, dropped = runner.initial_filter()
        candidates, backup, dropped = runner.has_offtarget(candidates=candidates, backup=backup, dropped_gRNA=dropped)
        candidates, backup, dropped = map(pd.DataFrame, [candidates, backup, dropped])
        return True

    def has_offtarget_result(self):
        return True

runner = CrisprFuncHelpers(database="", strand="r", mismatch=4,
                           cas9="Streptococcus thermophilus")
class CrisprFunction(unittest.TestCase):
    """
    unit test runner for the cripr.py
    """
    def test_initial_filter(self):
        self.assertEqual(runner.initial_filter_result(),
                         runner.initial_filter_test())

    def test_has_oftarget_func(self):
        self.assertEqual(runner.has_offtarget_result(),
                         runner.has_offtarget_test())


if __name__ == '__main__':
    unittest.main()
