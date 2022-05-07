import argparse
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import sys

inf = float("inf")


def group_hits(hits, A, N):
    """
    Finds all groups of N non-overlapping hits whose indices are within a distance of A from one another.

    :param hits: a list of seed-hit tuples (db_index, query_index, length)
    :param A: a distance threshold for grouping hits
    :param N: the number of hits in each group
    :return: a new dictionary of hits to their
    """
    hits_by_diagonal = {}
    for db_index, q_index, k in hits:
        diagonal = db_index - q_index  # this defines the diagonal
        if diagonal not in hits_by_diagonal:
            hits_by_diagonal[diagonal] = []
        hits_by_diagonal[diagonal].append((db_index, q_index, k))

    hit_groups = []
    for aligned_hits in hits_by_diagonal.values():
        aligned_hits = sorted(aligned_hits)
        start = 0
        while start <= len(aligned_hits) - N:
            hit_group = [aligned_hits[start]]
            valid_group = True
            for curr in range(start + 1, start + N):
                prev_i, prev_j, prev_k = aligned_hits[curr - 1]
                curr_i, curr_j, curr_k = aligned_hits[curr]
                if curr_i - prev_i > A:
                    valid_group = False
                    start = curr
                    break
                else:
                    hit_group.append((curr_i, curr_j, curr_k))
            if valid_group:
                hit_groups.append(hit_group)

    return hit_groups


def extract_hits(hits_file):
    """
    Extracts hits from a file and returns them as a list.

    :param hits_file: a file containing the hits
    :return: a list of hits tuples (i, j, k)
    """
    hits = []
    with open(hits_file, mode="r") as f:
        for line in f:
            if line[0] != "#":
                i, j, k = line.split()
                hits.append((int(i), int(j), int(k)))

    return hits


def display_hit_groups(hit_groups, db_seq=None, q_seq=None):
    """

    :param hit_groups:
    :param db_seq:
    :param q_seq:
    :return:
    """
    fig, ax = plt.subplots()
    fig.suptitle("BLAST Hit Groups\n")
    ax.set_xlabel("Database Sequence")
    ax.set_ylabel("Query Sequence")

    group_subsequences = []
    for group in hit_groups:
        first_hit = group[0]
        last_hit = group[-1]
        start_i, start_j, start_k = first_hit
        end_i, end_j, end_k = last_hit
        group_subsequences.append((start_i, start_j, (end_i - start_i) + end_k))

    group_alignments = [(i, j, i + k, j + k) for i, j, k in group_subsequences]
    lc = LineCollection(group_alignments, linewidths=1, zorder=0)
    ax.add_collection(lc)

    db_indices = []
    q_indices = []
    for group in hit_groups:
        for hit in group:
            i, j, k = hit
            db_indices.append((i))
            q_indices.append((j))

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
    description = "Finds all groups of N non-overlapping hits whose pairwise distances are less than A." \
                  "If an N is not provided then the default value is 2 (it looks for pairs)"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("hits_file", type=str)
    parser.add_argument("A", type=int)
    parser.add_argument("N", type=int, default=2)
    parser.add_argument("--display", "-d", action="store_true")
    args = parser.parse_args(sys.argv[1:])

    hits = extract_hits(args.hits_file)

    hit_groups = group_hits(hits, args.A, args.N)

    print(f"# Distance Limit (A): {args.A}")
    print(f"# Group Size (N): {args.N}")
    print(f"# {len(hit_groups)} Hit Groups:")
    for hit_group in hit_groups:
        hit_group = [" ".join(hit) for hit in hit_group]
        print(",".join(hit_group))

    if args.display:
        display_hit_groups(hit_groups)


if __name__ == "__main__":
    main()
