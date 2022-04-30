import sys
import argparse
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import collections as mc
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


# def display_hits(hits):
#     """
#
#     :param hits:
#     :return:
#     """
#     alignments = [[(i, j), (i + k, j + k)] for i, j, k in hits]
#     lc = mc.LineCollection(alignments, linewidths=1)
#     fig, ax = pl.subplots()
#     pl.show()


def main():
    parser = argparse.ArgumentParser(description="Locate HSPs in a database sequence.")
    parser.add_argument("db_seq", type=str)
    parser.add_argument("seed_file", type=str)
    args = parser.parse_args(sys.argv[1:])

    seeds = extract_seeds(args.seed_file)
    hits = find_hits(args.db_seq, seeds)

    print(f"# Database Sequence: {args.db_seq}")
    print(f"# Seeds: {seeds}")
    print(f"# Hits:")
    for db_index, query_index, k in hits:
        print(f"{db_index} {query_index} {k}")


if __name__ == "__main__":
    main()
