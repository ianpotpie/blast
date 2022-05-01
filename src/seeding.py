import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from knuth_morris_pratt import kmp_search


def find_hits(db_seq: str, seeds: list):
    """
    Finds instances of the seeds in the database sequence.

    :param db_seq: the database sequence
    :param seeds: a list of seeds (query_index, kmer)
    :return: a list of hits (db_index, query_index, length)
    """
    kmer_to_indices = {}  # dictionary prevents searching for the same kmer twice
    for query_index, kmer in seeds:
        if kmer not in kmer_to_indices:
            kmer_to_indices[kmer] = []
        kmer_to_indices[kmer].append(query_index)

    hits = []
    for seed, query_indices in kmer_to_indices.items():
        k = len(seed)
        db_indices = kmp_search(db_seq, seed)
        for db_index in db_indices:
            for query_index in query_indices:
                hits.append((db_index, query_index, k))

    return hits


def extract_seeds(seed_file):
    """

    :param seed_file:
    :return:
    """
    seeds = []
    with open(seed_file) as f:
        for line in f:
            if line[0] != "#":
                index, kmer = line.split()
                seeds.append((int(index), kmer))
    return seeds


def display_hits(hits, db_seq=None, query_seq=None):
    """

    :param hits: a sequence of hits to be plotted (db_index, query_index, length)
    :param db_seq: the database sequence from which the hits were taken
    :param query_seq: the query sequence from which the hits were taken
    :return:
    """
    fig, ax = plt.subplots()
    fig.suptitle("BLAST Seeding\n")
    ax.set_xlabel("Database Sequence")
    ax.set_ylabel("Query Sequence")
    alignments = [[(i, j), (i + k - 1, j + k - 1)] for i, j, k in hits]
    lc = LineCollection(alignments, linewidths=1, zorder=0)
    ax.add_collection(lc)
    db_indices = [i for i, _, _ in hits]
    query_indices = [j for _, j, _ in hits]
    ax.scatter(db_indices, query_indices, s=10, c="black", zorder=1)
    if db_seq is not None:
        ax.set_xticks([i for i in range(len(db_seq))])
        ax.set_xticklabels([s for s in db_seq])
    if query_seq is not None:
        ax.set_yticks([j for j in range(len(query_seq))])
        ax.set_yticklabels([s for s in query_seq])
    ax.tick_params(labelsize=4)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Locate HSPs in a database sequence.")
    parser.add_argument("db_seq", type=str)
    parser.add_argument("seed_file", type=str)
    parser.add_argument("--query-seq", type=str)
    args = parser.parse_args(sys.argv[1:])

    seeds = extract_seeds(args.seed_file)
    hits = find_hits(args.db_seq, seeds)

    print(f"# Database Sequence: {args.db_seq}")
    if args.query_seq is not None:
        print(f"# Query Sequence: {args.query_seq}")
    print(f"# Seeds: {seeds}")
    print(f"# Hits:")
    for db_index, query_index, k in hits:
        print(f"{db_index} {query_index} {k}")

    display_hits(hits, args.db_seq, args.query_seq)


if __name__ == "__main__":
    main()
