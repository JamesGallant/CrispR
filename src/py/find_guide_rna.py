
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
