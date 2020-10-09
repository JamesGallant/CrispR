import  numpy as np
class SmithWaterman:
    """
    Python and numpy implementation of the smith waterman algorithm.
    The algorithm has three main parts: Scoring, backtracking and calculating the start and end index
    """
    def __init__(self):
        pass

    def sequence_matrix(self, reference: str, query: str, score=3, gap_cost=2):
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
