import unittest
from src.py import find_guide_rna

class TestRunner:
    def __init__(self):
        self.runner = find_guide_rna.SmithWaterman()

    def test(self):
        out = self.runner.smith_waterman(reference="atgtcatcgggcaattcatctctgggaattatcgtcgggatcgacgattcaccggccgcacaggttgcggtgcggtgggcagctcgggatgcggagttgcgaaaaatccctctgacgctcgtgcacgcggtgtcgccggaagtagccacctggctggaggtgccactgccgccgggcgtgctgcgatggcagcaggatcacgggcgccacctgatcgacgacgcactcaaggtggttgaacaggcttcgctgcgcgctggtccccccacggtccacagtgaaatcgttccggcggcagccgttcccacattggtcgacatgtccaaagacgcagtgctgatggtcgtgggttgtctcggaagtgggcggtggccgggccggctgctcggttcggtcagttccggcctgctccgccacgcgcactgtccggtcgtgatcatccacgacgaagattcggtgatgccgcatccccagcaagcgccggtgctagttggcgttgacggctcgtcggcctccgagctggcgaccgcaatcgcattcgacgaagcgtcgcggcgaaacgtggacctggtggcgctgcacgcatggagcgacgtcgatgtgtcggagtggcccggaatcgattggccggcaactcagtcgatggccgagcaggtgctggccgagcggttggcgggttggcaggagcggtatcccaacgtagccataacccgcgtggtggtgcgcgatcagccggcccgccagctcgtccaacgctccgaggaagcccagctggtcgtggtcggcagccggggccgcggcggctacgccggaatgctggtggggtcggtaggcgaaaccgttgctcagctggcgcggacgccggtcatcgtggcacgcgagtcgctgacttag",
                                         query="atgtctcgggcacatctctggg")
        print(out)


    def result(self):
        pass


class FindgRNATest(unittest.TestCase):
    def test_smithwaterman(self):
        self.assertEqual(TestRunner().result(), TestRunner().test())


if __name__ == '__main__':
    unittest.main()
