import argparse
import sys
import matplotlib.pyplot as plt
from scoring_scheme import ScoringScheme
import numpy as np

inf = float("inf")


def gapped_extension(db_seq, q_seq, scoring_scheme, X):
    """
    Uses a variation of Gotoh's algorithm to perform local alignment with a fixed starting point.
    This variation implements a drop-off metric X and utilizes a dictionary/map in place of a dynamic programming
    matrix. Rather than iterating over and calculating the scores of all cells, it successively calculates the scores
    of cells adjacent to previously-scores cells. When the score of a cell drops more than X below the previously-seen
    maximum, does not calculate the scores of the cells adjacent to that cell. As a result, some cells will not have
    scores, and the dictionary prevents those cells from taking space in memory (unlike an array).

    :param db_seq: the database sequence which we are aligning
    :param q_seq: the query sequence which we are aligning
    :param scoring_scheme: the scoring system used in the alignment
    :param X: the drop-off threshold
    :return: a max alignment score and a list of alignments that produce that score
    """
    # initialize constants
    gap_open = scoring_scheme.gap_open
    gap_extend = scoring_scheme.gap_extend
    max_i = len(q_seq)
    max_j = len(db_seq)

    # initialize the scoring matrices
    M = dict()
    M[0, 0] = 0.0

    I = dict()
    I[0, 0] = -inf

    D = dict()
    D[0, 0] = -inf

    B = dict()
    B[0, 0] = 0.0

    # populate the dynamic programming matrix
    frontier = {(1, 0), (0, 1)}
    best_score, best_ends = 0.0, []
    while len(frontier) > 0:
        next_frontier = set()
        for i, j in frontier:

            if i == 0:
                M[i, j], I[i, j], D[i, j] = -inf, -inf, gap_open + (j * gap_extend)

            elif j == 0:
                M[i, j], I[i, j], D[i, j] = -inf, gap_open + (i * gap_extend), -inf

            else:
                M[i, j] = -inf if (i - 1, j - 1) not in M else \
                    max(M[i - 1, j - 1], I[i - 1, j - 1], D[i - 1, j - 1]) + scoring_scheme(db_seq[j - 1], q_seq[i - 1])

                I[i, j] = -inf if (i - 1, j) not in M else \
                    max(M[i - 1, j] + gap_open, I[i - 1, j], D[i - 1, j] + gap_open) + gap_extend

                D[i, j] = -inf if (i, j - 1) not in M else \
                    max(M[i, j - 1] + gap_open, I[i, j - 1] + gap_open, D[i, j - 1]) + gap_extend

            curr_score = max(M[i, j], D[i, j], I[i, j])

            B[i, j] = max(B[i - 1, j - 1] if (i - 1, j - 1) in B else -inf,
                          B[i, j - 1] if (i, j - 1) in B else -inf,
                          B[i - 1, j] if (i - 1, j) in B else -inf,
                          curr_score)

            # only extend a cell if its score has not dropped away from the maximum
            if (B[i, j] - curr_score < X) and (i + 1 <= max_i):
                next_frontier.add((i + 1, j))
            if (B[i, j] - curr_score < X) and (j + 1 <= max_j):
                next_frontier.add((i, j + 1))

            if curr_score > best_score:
                best_score = curr_score
                best_ends = []
            if curr_score == best_score:
                best_ends.append((i, j))

        frontier = next_frontier

    # backtrace to find the alignments that created the max alignment score
    alignments = []
    backtraces = [("", "", "M", i, j) for i, j in best_ends]
    while len(backtraces) > 0:
        next_backtraces = []
        for q_suffix, d_suffix, matrix, i, j in backtraces:
            if i == 0:
                alignments.append(((j * "-") + q_suffix, d_suffix))
            elif j == 0:
                alignments.append((q_suffix, (i * "-") + d_suffix))
            else:
                if matrix == "M":
                    prev_cell = (i - 1, j - 1)
                    q_symbol, d_symbol = q_seq[i - 1], db_seq[j - 1]
                    q_suffix, d_suffix = q_symbol + q_suffix, d_symbol + d_suffix
                    prev_score = M[i, j] - scoring_scheme(q_symbol, d_symbol)
                    if (prev_cell in M) and (M[prev_cell] == prev_score):
                        next_backtraces.append((q_suffix, d_suffix, "M", i - 1, j - 1))
                    if (prev_cell in I) and (I[prev_cell] == prev_score):
                        next_backtraces.append((q_suffix, d_suffix, "I", i - 1, j - 1))
                    if (prev_cell in D) and (D[prev_cell] == prev_score):
                        next_backtraces.append((q_suffix, d_suffix, "D", i - 1, j - 1))

                if matrix == "I":
                    prev_cell = (i - 1, j)
                    q_symbol, d_symbol = q_seq[i - 1], "-"
                    q_suffix, d_suffix = q_symbol + q_suffix, d_symbol + d_suffix
                    if prev_cell in M and I[i, j] == M[prev_cell] + gap_open + gap_extend:
                        next_backtraces.append((q_suffix, d_suffix, "M", i - 1, j))
                    if prev_cell in I and I[i, j] == I[prev_cell] + gap_extend:
                        next_backtraces.append((q_suffix, d_suffix, "I", i - 1, j))
                    if prev_cell in D and I[i, j] == D[prev_cell] + gap_open + gap_extend:
                        next_backtraces.append((q_suffix, d_suffix, "D", i - 1, j))

                if matrix == "D":
                    prev_cell = (i, j - 1)
                    q_symbol, d_symbol = "-", db_seq[j - 1]
                    q_suffix, d_suffix = q_symbol + q_suffix, d_symbol + d_suffix
                    if prev_cell in M and D[i, j] == M[prev_cell] + gap_open + gap_extend:
                        next_backtraces.append((q_suffix, d_suffix, "M", i, j - 1))
                    if prev_cell in I and D[i, j] == I[prev_cell] + gap_open + gap_extend:
                        next_backtraces.append((q_suffix, d_suffix, "I", i, j - 1))
                    if prev_cell in D and D[i, j] == D[prev_cell] + gap_extend:
                        next_backtraces.append((q_suffix, d_suffix, "D", i, j - 1))

        backtraces = next_backtraces

    return best_score, alignments, {"M": M, "I": I, "D": D, "B": B}


def main():
    description = "Takes each hit or hit group within a file and performs gapped extension to the left and right." \
                  "If the extension score drops more than X below its previous maximum, then the extension ceases and" \
                  "the end is set at the previous maximum. After extension, all the alignments scoring over S are" \
                  "printed."
    parser = argparse.ArgumentParser()
    parser.add_argument("database_file", type=str,
                        help="This file contains the database sequence in which the hits were found.")
    parser.add_argument("query_file", type=str,
                        help="This file contains the query sequence in which the hits were found.")
    parser.add_argument("hits_file", type=str,
                        help="This file contains the hits found between the database sequence and the query sequence."
                             "The hits are either individual (one per line), or in groups (separated by commas).")
    parser.add_argument("X", type=float, help="The dropoff threshold value.")
    parser.add_argument("--match", "-m", type=float, default=1.0,
                        help="The score of two matching symbols (used if scoring matrix is not provided).")
    parser.add_argument("--mismatch", "-e", type=float, default=-1.0,
                        help="The score of two mismatching symbols (used if scoring matrix is not provided).")
    parser.add_argument("--gap-open", "-o", type=float, default=-2.0,
                        help="The penalty for beginning an insertion or deletion")
    parser.add_argument("--gap-extend", "-g", type=float, default=-1.0,
                        help="The penalty for extending an insertion or deletion")
    parser.add_argument("--matrix", "-s", type=str, help="The scoring matrix for alignments")
    parser.add_argument("--display", "-d", action="store_true",
                        help="If this flag is set, then the extensions will be displayed using a matplot plot with"
                             "the original hits for each sequence displayed as well.")
    args = parser.parse_args(sys.argv[1:])

    db_seq = ""
    with open(args.database_file) as f:
        for line in f:
            if line[0] != "#":
                db_seq += line.strip()

    q_seq = ""
    with open(args.query_file) as f:
        for line in f:
            if line[0] != "#":
                q_seq += line.strip()

    seeds = []
    with open(args.hits_file, mode="r") as f:
        for line in f:
            if line[0] != "#":
                group = line.split(",")
                start_i, start_j, start_k = group[0].split()
                end_i, end_j, end_k = group[-1].split()
                i = int(start_i)
                j = int(start_j)
                k = (int(end_i) - int(start_i)) + int(end_k)
                seeds.append((i, j, k))

    scoring_scheme = ScoringScheme(args.match, args.mismatch, gap_open=args.gap_open, gap_extend=args.gap_extend)
    if args.matrix is not None:
        scoring_scheme.load_matrix(args.matrix)

    score = None
    extensions = set()
    for i, j, k in seeds:
        l_offset, r_offset = (k // 2) - 1, k // 2

        # in order to perform leftward extension, we flip the sequences and perform rightward extension
        # this saves us from defining a near-identical left-extending function
        db_prefix = db_seq[i + l_offset::-1]
        q_prefix = q_seq[j + l_offset::-1]
        l_score, l_extensions, l_matrices = gapped_extension(db_prefix, q_prefix, scoring_scheme, args.X)

        db_suffix = db_seq[i + r_offset:]
        q_suffix = q_seq[j + r_offset:]
        r_score, r_extensions, r_matrices = gapped_extension(db_suffix, q_suffix, scoring_scheme, args.X)

        score = l_score + r_score
        for l_extension in l_extensions:
            for r_extension in r_extensions:
                db_extension = l_extension[0][::-1] + r_extension[0]
                q_extension = l_extension[1][::-1] + r_extension[1]
                extensions.add((db_extension, q_extension))

        if args.display:
            display_matrix = np.zeros((len(q_seq), len(db_seq)))
            for r, c in l_matrices["B"].keys():
                prev_best = l_matrices["B"][r, c]
                curr_score = max(l_matrices["M"][r, c], l_matrices["I"][r, c], l_matrices["D"][r, c])
                display_matrix[j + l_offset - (r - 1), i + l_offset - (c - 1)] = \
                    1 - ((prev_best - curr_score) / (2 * args.X))
            for r, c in r_matrices["B"].keys():
                prev_best = r_matrices["B"][r, c]
                curr_score = max(r_matrices["M"][r, c], r_matrices["I"][r, c], r_matrices["D"][r, c])
                display_matrix[j + r_offset + (r - 1), i + r_offset + (c - 1)] = \
                    1 - ((prev_best - curr_score) / (2 * args.X))

            plt.imshow(display_matrix, cmap="Greys")
            plt.show()

    print(f"# Database Sequence: {db_seq}")
    print(f"# Query Sequence: {q_seq}")
    print(f"# Dropoff (X): {args.X}")
    print(f"# {len(extensions)} Extensions with the optimal score {score}: ")
    for db_extension, q_extension in extensions:
        print()
        print(db_extension)
        print(q_extension)


if __name__ == "__main__":
    main()
