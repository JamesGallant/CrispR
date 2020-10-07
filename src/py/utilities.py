import pandas as pd

from src.py.misc_functions import *

class CrispinatorUitls:
    def __init__(self, cas9: str):
        self.cas9 = cas9

    def rank_pams(self, dictionary):
        """:argument dictionary is a dict
        requires a PAM header
        :returns: dictionary
        """
        err = False if dictionary['PAM'] else True

        if err:
            raise AssertionError("Pams not present")

        df = pd.DataFrame(dictionary)
        ranking_dict = possible_pams_ranked(cas9=self.cas9)
        pam_rank = [ranking_dict.get(pams, None)[0] for pams in df['PAM']]
        df['rank'] = pam_rank
        return df.to_dict('list')

    def annotate_dataframe(self, dictionary=None):
        """:parameter dictionary
        Takes dicts or pandas dataframe. Requires a 'names' key
        :returns dictionary
        """
        dictionary['genes'] = [genes.split("_")[0] for genes in dictionary['names']]

        return dictionary
