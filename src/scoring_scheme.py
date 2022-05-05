import numpy as np


class ScoringScheme:
    def __init__(self, match=1.0, mismatch=-1.0, gap_extend=-1.0, gap_open=0.0, semi_global=False):
        self.match = match
        self.mismatch = mismatch
        self.gap_extend = gap_extend
        self.gap_open = gap_open
        self.semi_global = semi_global
        self.symbol_to_index = None
        self.scoring_matrix = None

    def get_symbols(self):
        """
        Creates a list of all symbols in the scoring matrix in the order in which they appear.
        Returns none if no scoring matrix has been loaded.

        :return: a list of symbols
        """
        if self.scoring_matrix is None:
            return None
        else:
            symbols = np.empty(len(self.symbol_to_index), dtype=str)
            for symbol, index in self.symbol_to_index.items():
                symbols[index] = symbol
            return list(symbols)

    def score(self, alignment1: str, alignment2: str):
        """
        Scores the alignment of two sequences based on their symbol-by-symbol scores going left-to-right.
        The method expects alignments to be the same length.
        It will quietly stop evaluating the score at the end of the shorter alignment.

        :param alignment1: the first sequence in the alignment to score
        :param alignment2: the second sequence in the alignment to score
        :return: the score of the alignment
        """
        if self.semi_global:
            alignment1, alignment2 = alignment1.lstrip("-"), alignment2.lstrip("-")  # disregard gaps at the beginning
            alignment1, alignment2 = alignment1[-min(len(alignment1), len(alignment2)):], \
                                     alignment2[-min(len(alignment1), len(alignment2)):]
            alignment1, alignment2 = alignment1.rstrip("-"), alignment2.rstrip("-")  # disregard gaps at the end
            alignment1, alignment2 = alignment1[:min(len(alignment1), len(alignment2))], \
                                     alignment2[:min(len(alignment1), len(alignment2))]

        score = 0.0
        for i, (symbol1, symbol2) in enumerate(zip(alignment1, alignment2)):

            if symbol1 == "-" and symbol2 == "-":
                raise ValueError("Encountered alignment between two gaps")

            elif symbol1 == "-":
                score += self.gap_extend
                score += self.gap_open if (i > 0) and alignment1[i - 1] != "-" else 0.0

            elif symbol2 == "-":
                score += self.gap_extend
                score += self.gap_open if (i > 0) and alignment2[i - 1] != "-" else 0.0

            elif self.scoring_matrix is not None:
                if (symbol1 not in self.symbol_to_index) or (symbol2 not in self.symbol_to_index):
                    raise ValueError("Encountered a symbol not in the scoring matrix")
                row = self.symbol_to_index[symbol1]
                col = self.symbol_to_index[symbol2]
                score += self.scoring_matrix[row][col]
            else:
                score += self.match if symbol1 == symbol2 else self.mismatch

        return score

    def __call__(self, alignment1, alignment2):
        """
        Calling the scoring scheme will score the two sequences passed in.
        You can find the behavior of the scoring in the "score" method.

        :param alignment1: the first sequence in the alignment to score
        :param alignment2: the second sequence in the alignment to score
        :return: the score of the alignment
        """
        return self.score(alignment1, alignment2)

    def load_matrix(self, filename):
        """
        Sets the scoring matrix of the scoring system based on the scoring matrix of a file.

        :param filename: the file containing the new matrix
        :return: None
        """
        with open(filename, mode='r') as f:
            headline = f.readline()
            while headline[0] == "#":  # iterates past all lines with comments
                headline = f.readline()
            self.symbol_to_index = {symbol.strip(): i for i, symbol in enumerate(headline.split())}

            # fill the scoring matrix
            n_symbols = len(self.symbol_to_index)
            self.scoring_matrix = np.zeros((n_symbols, n_symbols))
            for line in f:
                row = line.split()
                symbol = row.pop(0)
                i = self.symbol_to_index[symbol]
                for j, score in enumerate(row):
                    self.scoring_matrix[i, j] = float(score)

    def __str__(self):
        """
        Creates a string representation of the scoring scheme.
        If a scoring matrix is loaded, then it uses a standard PAM or BLOSUM style matrix.

        :return: a string of the scoring scheme
        """
        if self.scoring_matrix is None:
            s = f"Match Score: {self.match}\n" + \
                f"Mismatch Penalty: {self.mismatch}\n" + \
                f"Gap Open Penalty: {self.gap_open}\n" + \
                f"Gap Extension Penalty: {self.gap_extend}"
        else:
            s = f"# Gap Open Penalty: {self.gap_open}\n" + \
                f"# Gap Extension Penalty: {self.gap_extend}\n"

            symbols = self.get_symbols()
            s += "   " + "  ".join(symbols)
            for i in range(self.scoring_matrix.shape[0]):
                s += "\n" + symbols[i]
                for j in range(self.scoring_matrix.shape[1]):
                    score = self.scoring_matrix[i, j]
                    score = int(score) if score.is_integer() else score
                    s += f" {score}" if score < 0.0 else f"  {score}"

        return s
