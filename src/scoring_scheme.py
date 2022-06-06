import numpy as np


class ScoringScheme:
    """
    This class contains the template for a symbol and alignment scoring scheme
    (note that the "-" symbol is reserved for gaps).
    If the scoring scheme does not recognize a symbol, or if a scoring matrix has not been specified, then the scheme
    will use the default match/mismatch/gap scores respectively.
    A scoring matrix can be automatically loaded from a file.
    """

    def __init__(self, match_score: float = 1.0, mismatch_score: float = -1.0, gap_score: float = -1.0) -> None:
        self.symbol_to_index = None
        self.scoring_matrix = None
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score

    def get_symbols(self) -> list[str]:
        """
        Creates a list of all the symbols (excluding gaps) that are in the scoring matrix.
        The list is ordered based on the order that the symbols appear in the scoring matrix.

        :return: a list of symbols
        """
        symbols = ["" for _ in range(len(self.symbol_to_index))]
        for symbol, i in self.symbol_to_index.items():
            symbols[i] = symbol
        return symbols

    def score_symbols(self, x: str, y: str) -> float:
        """
        Scores the alignment of two symbols. "-" is reserved for gaps.

        :param x: the first symbol
        :param y: the second symbol
        :return: the score of the alignment
        """
        if (x == "-") and (y == "-"):
            return 0.0

        if (x == "-") or (y == "-"):
            return self.gap_score

        if (self.symbol_to_index is not None) and (x in self.symbol_to_index) and (y in self.symbol_to_index):
            r = self.symbol_to_index[x]
            c = self.symbol_to_index[y]
            return self.scoring_matrix[r][c]

        return self.match_score if x == y else self.mismatch_score

    def score_alignment(self, x: str, y: str) -> float:
        """
        Scores the alignment of two sequences based on their symbol-by-symbol scores
        (does not change based off the ordering of the symbols alignments).
        If x and y do not have the same length, then the score will treat the missing symbols as gaps.

        :param x: the first alignment sequence
        :param y: the second alignment sequence
        :return: the score the alignment of the two sequences
        """
        score = 0.0
        for i in range(min(len(x), len(y))):
            score += self.score_symbols(x[i], y[i])
        score += abs(len(x) - len(y)) * self.gap_score
        return score

    def load_matrix(self, filename: str) -> None:
        """
        Sets the scoring matrix of the scoring system based on the scoring matrix of a file.
        Files should be of the form:\n
        X A B C\n
        A 1 2 3\n
        B 4 5 6\n
        C 7 8 9\n

        :param filename: the file containing the new matrix
        :return: None
        """
        with open(filename, mode="r") as f:
            file_matrix = [line.strip().split() for line in f.readlines()]

            # takes the symbols from the first row
            self.symbol_to_index = {s: i for i, s in enumerate(file_matrix[0][1:])}

            # takes the score values from the interior of the matrix
            n_symbols = len(self.symbol_to_index)
            self.scoring_matrix = np.array((n_symbols, n_symbols))
            for r in range(n_symbols):
                for c in range(n_symbols):
                    self.scoring_matrix[r][c] = float(file_matrix[r + 1][c + 1])

    def get_lambda(self, prior, precision=0.001) -> float:
        """
        Calculates the value of the normalizing lambda value that is the unique solution to

        sum_ij(p_i * p_j * exp(s_ij * lambda))

        :param prior: the underlying distribution of the symbols
        :param precision: the precision of the lambda value
        :return: value of lambda for the scoring matrix
        """
        prior = np.reshape(np.array(prior), (1, -1))
        min_bound = 0.0

        max_bound = 1.0
        found_max = False
        while not found_max:
            total = prior.T @ np.exp(self.scoring_matrix * max_bound) @ prior
            if total > 1:
                found_max = True
            else:
                found_max *= 2

        while max_bound - min_bound >= precision:
            mid = (max_bound - min_bound) / 2
            total = prior.T @ np.exp(self.scoring_matrix * mid) @ prior
            if total > 1.0:
                max_bound = mid
            if total < 1.0:
                min_bound = mid

        return (max_bound - min_bound) / 2

    def get_transition_matrix(self, lamb: float, prior: np.ndarray) -> np.ndarray:
        """
        Get the transition probability matrix using the provided lambda

        :param lamb: the lambda value use to calculate probability values
        :param prior: and ndarray
        :return: the transition probability matrix`
        """
        prior = np.array(prior)
        probs = np.exp(np.copy(self.scoring_matrix) * lamb)
        probs = probs.T * prior
        probs = probs.T * prior
        for i in range(probs.shape[0]):
            probs[i] /= sum(probs[i])  # normalizes so that all rows sum to 1
        return probs
