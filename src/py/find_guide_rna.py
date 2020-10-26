import numpy as np
import itertools
from Bio import pairwise2
from Bio import Align

class SmithWaterman:
    """
    Python and numpy implementation of the smith waterman algorithm.
    The algorithm has three main parts: Scoring, backtracking and calculating the start and end index
    """

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

        return matrix

    def traceback(self, matrix, b, b_='', old_i=0):
        # flip H to get index of **last** occurrence of H.max() with np.argmax()
        H_flip = np.flip(np.flip(matrix, 0), 1)
        i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
        i, j = np.subtract(matrix.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
        if matrix[i, j] == 0:
            return b_, j

        b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
        return self.traceback(matrix[0:i, 0:j], b, b_, i)

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


