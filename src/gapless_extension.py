import argparse
import sys
from scoring_scheme import ScoringScheme
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def gapless_extension(db_seq, q_seq, seed, scoring_scheme, X):
    """
    Takes in the position of a seed of size k and returns the largest viable extension based on a scoring scheme and the
    falloff limit for the random walk score. The extension is provided in the form of its positions in both the
    search sequence and the query along with the length of the extension.

    :param db_seq: the sequence in which the seed was matched
    :param q_seq: the query from which the seed was derived
    :param seed: the location and length of an initial high-matching region
    :param scoring_scheme: a function for scoring symbols
    :param X: the limit for the score falloff when performing the extension
    :return: the left and right bounds of the extension of the seed
    """
    i, j, k = seed
    db_kmer = db_seq[i: i + k]
    q_kmer = q_seq[j: j + k]
    seed_score = scoring_scheme(db_kmer, q_kmer)

    # extend the query/database sequence to the right
    if (i + k) < len(db_seq) and (j + k) < len(q_seq):
        max_score = cur_score = seed_score
        best_index = -1
        for index, (db_symbol, q_symbol) in enumerate(zip(db_seq[i + k:], q_seq[j + k:])):
            if max_score - cur_score > X:
                break
            cur_score += scoring_scheme(db_symbol, q_symbol)
            if cur_score >= max_score:
                max_score = cur_score
                best_index = index
        k = k + best_index + 1

    # extend the query/database sequence to the left
    if i > 0 and j > 0:
        max_score = cur_score = seed_score
        best_index = -1
        for index, (db_symbol, q_symbol) in enumerate(zip(db_seq[i - 1::-1], q_seq[j - 1::-1])):
            if max_score - cur_score > X:
                break
            cur_score += scoring_scheme(db_symbol, q_symbol)
            if cur_score >= max_score:
                max_score = cur_score
                best_index = index
        i, j, k = i - best_index + 1, j - best_index + 1, k + best_index + 1

    return i, j, k


def display_gapless_extensions(extensions, seeds=None):
    """
    Displays a list of extensions using matplotlib.

    :param extensions: a sequence of hits to be plotted (db_index, query_index, length)
    :param seeds: the hits from which the extensions were created (optional)
    :return: None
    """
    fig, ax = plt.subplots()
    fig.suptitle("Gapless Extension", y=0.06)
    ax.set_xlabel("Database")
    ax.set_ylabel("Query")

    segments = []
    for i, j, k in extensions:
        segments.append(((i, j), (i + k, j + k)))
    ax.add_collection(LineCollection(segments, zorder=0))

    if seeds is not None:
        db_indices = [i + (k / 2) for i, _, k in seeds]
        query_indices = [j + (k / 2) for _, j, k in seeds]
        ax.scatter(db_indices, query_indices, s=1, c="r", zorder=1)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    plt.show()


def main():
    description = "Takes each hit or hit group within a file and performs gapless extension to the left and right." \
                  "If the extension score drops more than X below its previous maximum, then the extension ceases and" \
                  "the end is set at the previous maximum. After extension, all the sequences scoring over S are" \
                  "printed."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("database_file", type=str,
                        help="This file contains the database sequence in which the hits were found.")
    parser.add_argument("query_file", type=str,
                        help="This file contains the query sequence in which the hits were found.")
    parser.add_argument("hits_file", type=str,
                        help="This file contains the hits found between the database sequence and the query sequence."
                             "The hits are either individual (one per line), or in groups (separated by commas).")
    parser.add_argument("X", type=float, help="The dropoff threshold value.")
    parser.add_argument("S", type=float,
                        help="The threshold score for considering an extension to be statistically significant.")
    parser.add_argument("--display", "-d", action="store_true",
                        help="If this flag is set, then the extensions will be displayed using a matplot plot with"
                             "the original hits for each sequence displayed as well.")
    parser.add_argument("--match", "-m", type=float, default=1,
                        help="The score of two matching symbols (used if scoring matrix is not provided).")
    parser.add_argument("--mismatch", "-e", type=float, default=-1,
                        help="The score of two mismatching symbols (used if scoring matrix is not provided).")
    parser.add_argument("--matrix", "-s", type=str,
                        help="The file containing the scoring matrix to use for evaluating extensions.")
    args = parser.parse_args(sys.argv[1:])

    db_seq = ""
    with open(args.database_file, mode="r") as f:
        for line in f:
            db_seq += line.strip()

    q_seq = ""
    with open(args.query_file, mode="r") as f:
        for line in f:
            q_seq += line.strip()

    seeds = []
    with open(args.hits_file, mode="r") as f:
        for line in f:
            if line[0] != "#":
                group = line.split(",")
                first_i, first_j, first_k = group[0].split()
                last_i, last_j, last_k = group[-1].split()
                k = int(last_k) + int(last_i) - int(first_i)
                seeds.append((int(first_i), int(first_j), k))

    scoring_scheme = ScoringScheme(match=args.match, mismatch=args.mismatch)
    if args.matrix is not None:
        scoring_scheme.load_matrix(args.matrix)

    extensions = set()
    for seed in seeds:
        extension = gapless_extension(db_seq, q_seq, seed, scoring_scheme, args.X)
        extensions.add(extension)

    filtered_extensions = []
    for i, j, k in extensions:
        if scoring_scheme(db_seq[i:i + k], q_seq[j:j + k]) > args.S:
            filtered_extensions.append((i, j, k))

    print(f"# Database Sequence: {db_seq}")
    print(f"# Query Sequence: {q_seq}")
    print(f"# Seeds: {seeds}")
    print(f"# Dropoff Threshold (X): {args.X}")
    print(f"# Extension Threshold (S): {args.S}")
    print(f"# {len(filtered_extensions)} Extensions:")
    for i, j, k in filtered_extensions:
        print(f"{i} {j} {k}")

    if args.display:
        display_gapless_extensions(extensions, seeds)
        display_gapless_extensions(filtered_extensions, seeds)


if __name__ == "__main__":
    main()
