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
    def result(self):
        pass


class FindgRNATest(unittest.TestCase):
    def test_smithwaterman(self):
        self.assertEqual(TestRunner().result(), TestRunner().test())


if __name__ == '__main__':
    unittest.main()
