import numpy as np
import itertools
from Bio import pairwise2
from Bio import Align


from Bio import Align

class Pairwise:
    """
    Performs paiwise alignment using smith waterman type algos.
    param: reference: String reference sequence to match against
    param: query: String sequence to match, the query sequence
    param: config: A dictionary containing the following keys {mode, extend_gap_score, target_end_gap_score,
    query_end_gap_score}
    defaults from https://biopython.org/docs/1.75/api/Bio.Align.html
    """
    
    def __init__(self, reference: str, query: str, config: dict):
        self.aligner = Align.PairwiseAligner()
        self.reference = reference
        self.query = query
        self.aligner.mode = config.get('mode', 'global')
        self.aligner.open_gap_score = config.get('open_gap_score', -0.5)
        self.aligner.extend_gap_score = config.get('extend_gap_score', -0.1)
        self.aligner.target_end_gap_score = config.get('target_end_gap_score', 0.0)
        self.aligner.query_end_gap_score = config.get('query_end_gap_score', 0.0)


    def __init__(self):
        pass

    def sequence_matrix(self, reference: str, query: str, score=2, gap_cost=2):
        """"
        :parameter reference: the reference sequence, i.e. longer strings
        :parameter query": The query sequence, i.e. smaller string within longer string
        :parameter score: default is three form wiki. If there is a gap = -3, match = 3
        :parameter gap_cost: default is two from wiki. If there is a gap on each adjacent cell add -2 penalty
        :returns a numpy 2D matrix
             A  T  G  C
         [0][0][0][0][0]
       A [0][3][0][0][0]
       C [0][0][0][0][3]
       G [0][0][0][3][0]
       C [0][0][0][0][3]
        """
        matrix = np.zeros((len(reference) + 1, len(query) + 1), np.int)

        for reference_idx, query_idx in itertools.product(range(1, matrix.shape[0]), range(1, matrix.shape[1])):
            match = matrix[reference_idx - 1, query_idx - 1] + (
                score if reference[reference_idx - 1] == query[query_idx - 1] else -score)
            delete = matrix[reference_idx - 1, query_idx] - gap_cost
            insert = matrix[reference_idx, query_idx - 1] - gap_cost
            matrix[reference_idx, query_idx] = max(match, delete, insert, 0)


        self.reference = reference
        self.query = query

    def align(self):
        """
        returns alignment object
        """
        return self.aligner.align(self.reference, self.query)


    def coordinates(self):
        """
        returns the start and stop coordinates as a list of tuples (start, stop) from the reference
        """
        return [objects.aligned[0][0] for objects in sorted(self.align())]

    def smith_waterman(self, reference, query, match_score=3, gap_cost=2):
        reference, query = reference.upper(), query.upper()
        matrix = self.sequence_matrix(reference, query, match_score, gap_cost)
        b_, pos = self.traceback(matrix, query)
        print(b_)
        return pos, pos + len(b_)


class Pairwise:
    def __init__(self, reference: str, sequence: str):
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.open_gap_score = -0.5
        self.aligner.extend_gap_score = -0.1
        self.aligner.target_end_gap_score = 0.0
        self.aligner.query_end_gap_score = 0.0

        self.reference = reference
        self.sequence = sequence

    def match(self):
        return self.aligner.align(self.reference, self.sequence)

