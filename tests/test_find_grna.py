import unittest
from src.py import find_guide_rna

reference = "atgccgaacctcgacacgaacacccagaatcggaccaacttgttacgggcggcgctcatcgccgtcggcttgatcttcaccttcgccgtctaccccctcaccatcatctggccatccggctggagttggggtcacgggtcatcgcattacctcacgatgatcatcggcatctacgccacacttggcgtcttcctgctgatcgcggcacgagatccgctggcccaccggagcctgatctggttcacggtggtatcgagcgtcgttcacgccgcgatcatggccgcgcaggcgatcggcgacccacacgagcgcggccacctggccggcgatgtgcccgctctggtgatcgtcgccgttgccctgggtctgttgatgcggggcgccgaaaccccgcgacgagtccagacacccgggtag"
primer = "CTCCCGCCGGAAGCGCAGATAGAAG".lower()
class TestRunner:
    def __init__(self):
        self.runner = find_guide_rna.SmithWaterman()

        self.pairwise = find_guide_rna.Pairwise(reference="ttttttggacctttttggaccattt",
                                                query="ggactt",
                                                config={'mode': 'global',
                                                        'open_gap_score': -0.5,
                                                        'extend_gap_score': -0.1,
                                                        'target_end_gap_score': 0.0,
                                                        'query_end_gap_score': 0.0
                                                        })

    def test(self):
        alignment = self.pairwise.align()
        coords = self.pairwise.coordinates()
        ref = "ttttttggatcttttt"
        print("testing")
        print(coords)
        print(f"start: {coords[0][0]}, stop: {coords[0][1]}")
        print(f"sequence: {ref[coords[0][0]:coords[0][1]]}")

        for objects in alignment:
            print(f"score: {objects.score}")
            print(objects)

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
