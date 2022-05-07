import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from knuth_morris_pratt import kmp_search


def find_hits(db_seq, seeds):
    """
    Finds instances of the seeds in the database sequence.

    :param db_seq: the database sequence
    :param seeds: a list of seeds (query_index, kmer)
    :return: a list of hits (db_index, query_index, length)
    """
    kmer_to_indices = {}  # dictionary prevents searching for the same kmer twice
    for q_index, kmer in seeds:
        if kmer not in kmer_to_indices:
            kmer_to_indices[kmer] = []
        kmer_to_indices[kmer].append(q_index)

    hits = []
    for seed, q_indices in kmer_to_indices.items():
        k = len(seed)
        db_indices = kmp_search(db_seq, seed)
        for db_index in db_indices:
            for q_index in q_indices:
                hits.append((db_index, q_index, k))

    return hits


def extract_seeds(seeds_file):
    """
    Takes a file containing a list of seeds and converts the data in the file to a list of seeds (index, kmer tuples).
    In the file, lines containing comments will begin with the "#" symbol.
    Each line will contain at most one seed, which will have the index and kmer of the seed separated by whitespace.

    :param seeds_file:
    :return:
    """
    seeds = []
    with open(seeds_file) as f:
        for line in f:
            if line[0] != "#":
                index, kmer = line.split()
                seeds.append((int(index), kmer))
    return seeds


def display_hits(hits, db_seq=None, q_seq=None):
    """
    Displays a dot-plot of the hits resulting from an alignment using matplotlib.

    :param hits: a sequence of hits to be plotted (db_index, query_index, length)
    :param db_seq: the database sequence from which the hits were taken
    :param q_seq: the query sequence from which the hits were taken
    :return:
    """
    fig, ax = plt.subplots()
    fig.suptitle("BLAST Hit Identification\n")
    ax.set_xlabel("Database Sequence")
    ax.set_ylabel("Query Sequence")
    alignments = [[(i, j), (i + k - 1, j + k - 1)] for i, j, k in hits]
    lc = LineCollection(alignments, linewidths=1, zorder=0)
    ax.add_collection(lc)
    db_indices = [i for i, _, _ in hits]
    q_indices = [j for _, j, _ in hits]
    ax.scatter(db_indices, q_indices, s=10, c="black", zorder=1)
    if db_seq is not None:
        ax.set_xticks([i for i in range(len(db_seq))])
        ax.set_xticklabels([s for s in db_seq])
    if q_seq is not None:
        ax.set_yticks([j for j in range(len(q_seq))])
        ax.set_yticklabels([s for s in q_seq])
    ax.tick_params(labelsize=4)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    plt.show()


def main():
    description = "Given a database sequence and a file containing a list of seeds, this script finds all instances " \
                  "of each seed in the database sequence and returns them as a list of \'hits\'. If the display flag " \
                  "is present, then a plot of the alignment will be shown. The original query sequence can be " \
                  "provided as well to improve the plot."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("database_sequence", type=str)
    parser.add_argument("seeds_file", type=str)
    parser.add_argument("--query-sequence", "-q", type=str)
    parser.add_argument("--display", "-d", type=bool, action="store_true")
    args = parser.parse_args(sys.argv[1:])

    seeds = extract_seeds(args.seeds_file)

    hits = find_hits(args.database_sequence, seeds)

    print(f"# Database Sequence: {args.database_sequence}")
    if args.query_sequence is not None:
        print(f"# Query Sequence: {args.query_sequence}")
    print(f"# Seeds: {seeds}")
    print(f"# Hits:")
    for db_index, q_index, k in hits:
        print(f"{db_index} {q_index} {k}")

    if args.display:
        display_hits(hits, args.database_sequence, args.query_sequence)


if __name__ == "__main__":
    main()
