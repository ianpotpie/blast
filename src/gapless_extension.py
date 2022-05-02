import argparse
import sys
from scoring_scheme import ScoringScheme
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt


def filter_extensions_by_score(extensions, db_seq, query_seq, scoring_scheme, S1):
    """
    Filters out all extensions from a list of extensions based on their scores

    :param extensions: a list of extensions
    :param db_seq: the database sequence form which the extensions were taken
    :param query_seq: the query sequence from which the extensions were taken
    :param scoring_scheme: the scoring scheme used to score the extensions
    :param S1: the threshold score used to discriminate between extensions
    :return: the list of extensions that surpass the threshold score
    """
    filtered_extensions = []
    for i, j, k in extensions:
        if scoring_scheme(db_seq[i:i + k], query_seq[j:j + k]) > S1:
            filtered_extensions.append((i, j, k))
    return filtered_extensions


def gapless_extension(db_seq, db_start, query_seq, query_start, k, scoring_scheme, X):
    """
    Takes in the position of a seed of size k and returns the largest viable extension based on a scoring scheme and the
    falloff limit for the random walk score. The extension is provided in the form of its positions in both the
    search sequence and the query along with the length of the extension.

    :param db_seq: the sequence in which the seed was matched
    :param db_start: the index of the seed location in the database sequence
    :param query_seq: the query from which the seed was derived
    :param query_start: the index of the seed location in the query sequence
    :param k: the length of the seed
    :param scoring_scheme: a function for scoring symbols
    :param X: the limit for the score falloff when performing the extension
    :return: the left and right bounds of the extension of the seed
    """
    db_kmer = db_seq[db_start: db_start + k]
    query_kmer = query_seq[query_start: query_start + k]
    seed_score = scoring_scheme(db_kmer, query_kmer)

    # extend the query/database sequence to the right
    max_score = cur_score = seed_score
    for i, (db_symbol, query_symbol) in enumerate(zip(db_seq[db_start + k:], query_seq[query_start + k:])):
        if max_score - cur_score > X:
            break
        cur_score += scoring_scheme(db_symbol, query_symbol)
        if cur_score >= max_score:
            max_score = cur_score
            k = k + (i + 1)

    # extend the query/database sequence to the left
    max_score = cur_score = seed_score
    for i, (db_symbol, query_symbol) in enumerate(zip(db_seq[db_start - 1::-1], query_seq[db_start - 1::-1])):
        if max_score - cur_score > X:
            break
        cur_score += scoring_scheme(db_symbol, query_symbol)
        if cur_score >= max_score:
            max_score = cur_score
            db_start, query_start, k = db_start - (i + 1), query_start - (i + 1), k + (i + 1)

    return db_start, query_start, k


# This is not really used in any papers that I found.
# def join_hits(hits, A):
#     """
#     Combines all hits within a distance of A from one another.
#
#     :param hits: a list of seed-hit tuples (db_index, query_index, length)
#     :param A: a distance threshold for merging hits
#     :return: a new dictionary of hits to their
#     """
#     hits_by_diagonal = {}
#     for db_index, query_index, k in hits:
#         diagonal = db_index - query_index  # this defines the diagonal
#         if diagonal not in hits_by_diagonal:
#             hits_by_diagonal[diagonal] = []
#         hits_by_diagonal[diagonal].append((db_index, query_index, k))
#
#     joined_hits = []
#     for aligned_hits in hits_by_diagonal.values():
#         prev_i, prev_j = float("-inf"), float("-inf")
#         for i, j, k in sorted(aligned_hits):
#             if i - prev_i > A:
#                 joined_hits.append((i, j, k))
#             else:
#                 start_i, start_j, combined_k = joined_hits[-1]
#                 joined_hits[-1] = (start_i, start_j, max(combined_k, (i + k) - start_i))
#             prev_i, prev_j = i, j
#
#     return joined_hits


def display_gapless_extensions(extensions, db_seq=None, query_seq=None, hits=None):
    """

    :param extensions: a sequence of hits to be plotted (db_index, query_index, length)
    :param db_seq: the database sequence from which the hits were taken (optional)
    :param query_seq: the query sequence from which the hits were taken (optional)
    :param hits: the hits from which the extensions were created (optional)
    :return: None
    """
    fig, ax = plt.subplots()
    fig.suptitle("Gapless Extension\n")
    ax.set_xlabel("Database Sequence")
    ax.set_ylabel("Query Sequence")
    alignments = [[(i, j), (i + k - 1, j + k - 1)] for i, j, k in extensions]
    lc = LineCollection(alignments, linewidths=1, zorder=0)
    ax.add_collection(lc)
    if db_seq is not None:
        ax.set_xticks([i for i in range(len(db_seq))])
        ax.set_xticklabels([s for s in db_seq])
    if query_seq is not None:
        ax.set_yticks([j for j in range(len(query_seq))])
        ax.set_yticklabels([s for s in query_seq])
    if hits is not None:
        db_indices = [i for i, j, k in extensions]
        query_indices = [j for i, j, k in extensions]
        ax.scatterplot(db_indices, query_indices, c="black", s=3)
    ax.tick_params(labelsize=4)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    plt.show()


def extract_hits(hit_file):
    """
    Takes a file containing a list of seeds and converts the data in the file to a list of seeds (index, kmer tuples).
    In the file, lines containing comments will begin with the "#" symbol.
    Each line will contain at most one seed, which will have the index and kmer of the seed separated by whitespace.

    :param hit_file: file containing hit data
    :return: the list of hits in the file
    """
    hits = []
    with open(hit_file) as f:
        for line in f:
            if line[0] != "#":
                db_index, query_index, k = line.split()
                hits.append((int(db_index), int(query_index), int(k)))
    return hits


def main():
    parser = argparse.ArgumentParser("Combine and/or extend all hits.")
    parser.add_argument("db_seq")
    parser.add_argument("query_seq")
    parser.add_argument("hit_file")
    parser.add_argument("X")
    parser.add_argument("S1")
    parser.add_argument("--match", "-m", type=float, default=1)
    parser.add_argument("--mismatch", "-n", type=float, default=-1)
    parser.add_argument("--matrix", "-s", type=str)
    args = parser.parse_args(sys.argv[1:])

    hits = extract_hits(args.hit_file)

    scoring_scheme = ScoringScheme(match=args.match, mismatch=args.mismatch)
    if args.matrix is not None:
        scoring_scheme.load_matrix(args.matrix)

    extensions = []
    for i, j, k in hits:
        extension = gapless_extension(args.db_seq, i, args.query_seq, j, k, scoring_scheme, args.X)
        extensions.append(extension)

    print(f"# Database Sequence: {args.db_seq}")
    print(f"# Query Sequence: {args.query_seq}")
    print(f"# Hits: {hits}")
    print(f"# Dropoff Threshold (X): {args.X}")
    print(f"# Extension Threshold (S1): {args.S1}")
    print(f"# Extensions:")
    for i, j, k in extensions:
        print(f"{i} {j} {k}")

    display_gapless_extensions(extensions, args.db_seq, args.query_seq, hits)


if __name__ == "__main__":
    main()