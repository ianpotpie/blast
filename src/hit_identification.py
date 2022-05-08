import sys
import argparse
import matplotlib.pyplot as plt
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


def display_hits(hits):
    """
    Displays a dot-plot of the hits resulting from an alignment using matplotlib.

    :param hits: a sequence of hits to be plotted (db_index, query_index, length)
    :return: None
    """
    fig, ax = plt.subplots()
    fig.suptitle("Hit Identification", y=0.06)
    ax.set_xlabel("Database")
    ax.set_ylabel("Query")

    db_indices = [i for i, _, _ in hits]
    q_indices = [j for _, j, _ in hits]
    ax.scatter(db_indices, q_indices, s=1, c="b")
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
    parser.add_argument("database_file", type=str,
                        help="This file contains the database sequence in which we are looking for hits."
                             "Characters on multiple lines are concatenated together (interpreted as one sequence).")
    parser.add_argument("seeds_file", type=str,
                        help="This file contains the seeds for the hits.")
    parser.add_argument("--display", "-d", action="store_true",
                        help="If this tag is included, the hits will be displayed in a matlab plot.")
    args = parser.parse_args(sys.argv[1:])

    db_seq = ""
    with open(args.database_file) as f:
        for line in f:
            db_seq += line.strip()

    seeds = []
    with open(args.seeds_file) as f:
        for line in f:
            if line[0] != "#":
                index, kmer = line.split()
                seeds.append((int(index), kmer))

    hits = find_hits(db_seq, seeds)

    print(f"# Database Sequence: {db_seq}")
    print(f"# Seeds: {seeds}")
    print(f"# {len(hits)} Hits:")
    for db_index, q_index, k in hits:
        print(f"{db_index} {q_index} {k}")

    if args.display:
        display_hits(hits)


if __name__ == "__main__":
    main()
