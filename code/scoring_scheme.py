import numpy as np


class ScoringScheme:
    def __init__(self, gap_extend=-1, gap_start=0.0, match_score=1.0, mismatch_penalty=-1.0):
        self.index_by_symbol = None
        self.scoring_matrix = None
        self.gap_extend = gap_extend
        self.gap_start = gap_start
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty

    def get_symbols(self):
        """
        Gets all the symbols (excluding gaps) that are in the scoring system.

        :return: a list of symbols
        """
        if self.index_by_symbol is None:
            return None
        else:
            symbols = [None for _ in self.index_by_symbol]
            for symbol, index in self.index_by_symbol.items():
                symbols[index] = symbol
            return symbols

    def score(self, sequence1, sequence2):
        """
        Scores the alignment of two sequences based on their symbol-by-symbol scores going left-to-right.
        The method expects the sequences to be the same length, but otherwise will quietly cut off the extra symbols
        from the longer alignment.

        :param sequence1: the first sequence in the alignment to score
        :param sequence2: the second sequence in the alignment to score
        :return: the score of the alignment
        """
        score = 0.0
        seq1_gap = False
        seq2_gap = False
        for symbol1, symbol2 in zip(sequence1, sequence2):
            assert symbol1 != "-" or symbol2 != "-", "Encountered alignment between two gaps"

            if symbol1 == "-":
                score += self.gap_extend if seq1_gap else self.gap_start + self.gap_extend
                seq1_gap, seq2_gap = True, False

            elif symbol2 == "-":
                score += self.gap_extend if seq2_gap else self.gap_start + self.gap_extend
                seq1_gap, seq2_gap = False, True

            elif self.scoring_matrix is not None:
                assert symbol1 in self.index_by_symbol, "Encountered symbol not in the current scoring matrix"
                assert symbol2 in self.index_by_symbol, "Encountered symbol not in the current scoring matrix"
                row = self.index_by_symbol[symbol1]
                col = self.index_by_symbol[symbol2]
                score += self.scoring_matrix[row][col]
            else:
                score += self.match_score if symbol1 == symbol2 else self.mismatch_penalty
        return score

    def load_matrix(self, filename):
        """
        Sets the scoring matrix of the scoring system based on the scoring matrix of a file.

        :param filename: the file containing the new matrix
        :return: None
        """
        with open(filename, mode='r') as f:
            # reads the first line of the file, disregarding the first comment lines
            headline = f.readline()
            while headline[0] == "#":  # iterates past all lines with comments
                headline = f.readline()
            self.index_by_symbol = {symbol.strip(): i for i, symbol in enumerate(headline.split())}

            # fill the scoring matrix
            n_symbols = len(self.index_by_symbol)
            self.scoring_matrix = np.zeros((n_symbols, n_symbols))
            for line in f:
                row = line.split()
                symbol = row.pop(0)
                i = self.index_by_symbol[symbol]
                for j, val in enumerate(row):
                    self.scoring_matrix[i][j] = float(val)
