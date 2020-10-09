import pandas as pd

class CrispinatorUitls:
    def __init__(self):

        self.complement_base = {
            'g': 'c',
            'G': 'C',
            'a': 't',
            'A': 'T',
            'c': 'g',
            'C': 'G',
            't': 'a',
            'T': 'A'
        }

    def possible_pams_ranked(self, cas9):
        """
        pam is the dictionary key, values at pos 0: rank, pos1: fold repression, pos2: scale pos3: wildcard bases
        :returns dict
        """
        if cas9 == "Streptococcus thermophilus":
            return {
                'AGAAG': [1, "high", 216.7, 2],
                'AGAAT': [2, "high", 216.2, 2],
                'AGAAA': [3, "high", 158.1, 2],
                'GGAAG': [4, "high", 145.2, 2],
                'AGAAC': [5, "high", 120.5, 2],
                'GGAAA': [6, "high", 110.5, 2],
                'AGCAT': [7, "high", 84.6, 2],
                'AGGAG': [8, "high", 82.2, 2],
                'AGGAT': [9, "high", 64.7, 2],
                'AGCAA': [10, "high", 53.4, 2],
                'GGAAC': [11, "high", 51.5, 2],
                'GGAAT': [12, "low", 47.3, 2],
                'AGCAG': [13, "low", 42.2, 2],
                'AGGAA': [14, "low", 38.5, 2],
                'AGGAC': [15, "low", 25.5, 2],
                'GGGAG': [16, "low", 24.7, 2],
                'GGGAT': [17, "low", 24.2, 2],
                'GGGAA': [18, "low", 12.3, 2],
                'AGCAC': [19, "low", 11.9, 2],
                'GGGAC': [20, "low", 7.9, 2],
                'GGCAT': [21, "low", 6.7, 2],
                'GGCAG': [22, "low", 4.0, 2],
                'GGCAA': [23, "low", 3.3, 2],
                'GGCAC': [24, "low", 2.7, 2]
            }

    def rank_pams(self, dictionary, cas9):
        """:argument dictionary is a dict
        requires a PAM header
        :returns: dictionary
        """
        err = False if dictionary['PAM'] else True

        if err:
            raise AssertionError("Pams not present")

        df = pd.DataFrame(dictionary)
        ranking_dict = self.possible_pams_ranked(cas9=cas9)
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

    def complement(self, sequence: str):
        """:parameter sequence DNA sequence to complement must be AGTC
        :returns string
        """
        valid_strings = ['A', 'a', 'G', 'g', 'C', 'c', 'T', 't']
        err = [False if string in valid_strings else True for string in sequence]
        if True in err:
            raise ValueError(f"None DNA characters in {sequence}")

        return "".join([self.complement_base.get(base, "N") for base in sequence])

    def reverse_compliment(self, sequence: str):
        """:parameter sequence DNA sequence to complement must be AGTC also reverses the sequence
        :returns string
        """
        valid_strings = ['A', 'a', 'G', 'g', 'C', 'c', 'T', 't']

        err = [False if string in valid_strings else True for string in sequence]
        if True in err:
            raise ValueError(f"None DNA characters in {sequence}")

        revcomp = "".join([self.complement_base.get(base, "N") for base in sequence])
        return revcomp[:: -1]
