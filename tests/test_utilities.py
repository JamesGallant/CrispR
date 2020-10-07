import unittest

from src.py import utilities

class UtilsRunner:
    def __init__(self, cas9: str):
        self.cas9 = cas9
        self.crispruitls = utilities.CrispinatorUitls(cas9=self.cas9)

    def test_pam_ranking(self):
        test = {
            'PAM': ["AGAAG", "AGAAC", "AGAAT"]
        }
        out = self.crispruitls.rank_pams(dictionary=test)

        return out.get('rank', None)

    def result_pam_ranking(self):
        return [1, 5, 2]

    def test_anotate_dataframe(self):
        test = {
            'names': ["Rv0001_gr15"]
        }
        out = self.crispruitls.annotate_dataframe(dictionary=test)
        return out.get('genes', None)

    def result_annotate_dataframe(self):
        return ["Rv0001"]

utils = UtilsRunner(cas9="Streptococcus thermophilus")
class UtilitiesTestCase(unittest.TestCase):
    def test_pam_ranking_func(self):
        self.assertEqual(utils.result_pam_ranking(),
                         utils.test_pam_ranking())

    def test_annotate_dataframe_func(self):
        self.assertEqual(utils.result_annotate_dataframe(),
                         utils.test_anotate_dataframe())

if __name__ == '__main__':
    unittest.main()
