import argparse
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
    for db_index, query_index, k in hits:
        diagonal = db_index - query_index  # this defines the diagonal
        if diagonal not in hits_by_diagonal:
            hits_by_diagonal[diagonal] = []
        hits_by_diagonal[diagonal].append((db_index, query_index, k))

    hit_groups = []
    for aligned_hits in hits_by_diagonal.values():
        aligned_hits = sorted(aligned_hits)
        index = 0
        while index <= len(aligned_hits) - N:
            hit_group = [aligned_hits[index]]
            next_index, valid_group = index + 1, True
            for ith_hit, hit in enumerate(aligned_hits[index + 1:index + N]):
                prev_i, prev_j, prev_k = hit_group[-1]
                curr_i, curr_j, curr_k = hit
                if curr_i - prev_i > A:
                    valid_group = index + ith_hit + 1, False
                    break
            if valid_group:
                hit_groups.append(hit_group)

    return hit_groups


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hits_file", type=str)
    parser.add_argument("A", type=int)
    parser.add_argument("--group-size", "-n", type=int, default=2)
    args = parser.parse_args(sys.argv[1:])

    hits = []
    with open(args.hits_file, mode="w") as f:
        for line in f:
            if line[0] != "#":
                i, j, k = line.split()
                hits.append((int(i), int(j), int(k)))

    hit_groups = group_hits(hits, args.A, args.group_size)

    print(f"# Distance Limit (A): {args.A}")
    print(f"# Group Size (N): {args.group_size}")
    print(f"# {len(hit_groups)} Hit Groups:")
    for hit_group in hit_groups:
        print(" ".join(hit_group) + "\n")


if __name__ == "__main__":
    main()
