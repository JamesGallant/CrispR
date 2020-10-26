import unittest
from src.py import find_guide_rna

class TestRunner:
    def __init__(self):
        self.runner = find_guide_rna.SmithWaterman()
        self.pairwise = find_guide_rna.Pairwise(reference="ttttagctattttagctattt",
                                                sequence="agcaa")

    def test(self):
        out = self.pairwise.match()

        #out = self.runner.smith_waterman(reference="ggggatggggag",
        #                                 query="atgaa")
        print("out")
        for alignments in sorted(out):
            print(f"mismatch = {alignments.score}")
            seq_idx = alignments.aligned
            start_ref = seq_idx[0][0][0]
            stop_ref = seq_idx[0][0][1]
            start_ref = seq_idx[0][0][0]
            stop_ref = seq_idx[0][0][1]
            print(str("agcaa")[start:stop])


    def result(self):
        pass


class FindgRNATest(unittest.TestCase):
    def test_smithwaterman(self):
        self.assertEqual(TestRunner().result(), TestRunner().test())


if __name__ == '__main__':
    unittest.main()
