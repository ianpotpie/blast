import argparse
import matplotlib.pyplot as plt
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
        diagonal = db_index - q_index  # the index difference is unique to a diagonal
        if diagonal not in hits_by_diagonal:
            hits_by_diagonal[diagonal] = []
        hits_by_diagonal[diagonal].append((db_index, q_index, k))

    hit_groups = []
    for aligned_hits in hits_by_diagonal.values():
        aligned_hits = sorted(aligned_hits)
        start = 0  # tracking the start index allows us to avoid re-evaluating subsequences of invalid groups
        while start <= len(aligned_hits) - N:
            hit_group = [aligned_hits[start]]
            valid_group = True
            for curr in range(start + 1, start + N):
                prev_i, prev_j, prev_k = aligned_hits[curr - 1]
                curr_i, curr_j, curr_k = aligned_hits[curr]
                if curr_i - prev_i <= A:
                    hit_group.append((curr_i, curr_j, curr_k))
                else:
                    valid_group = False
                    start = curr  # if a gap exceeds A, then we can move the start index over the gap
                    break

            if valid_group:
                hit_groups.append(hit_group)
                start += 1

    return hit_groups


def display_hit_groups(hit_groups, original_hits=None):
    """
    Displays the hit groups (and original hits if provided) on a matplot scatter-plot.

    :param hit_groups: a list of lists of hits (tuples)
    :param original_hits: a list of tuples
    :return: None
    """
    fig, ax = plt.subplots()
    fig.suptitle("Hit Groups", y=0.06)
    ax.set_xlabel("Database")
    ax.set_ylabel("Query")

    if original_hits is not None:
        db_indices = [i + (k / 2) for i, _, k in original_hits]
        q_indices = [j + (k / 2) for _, j, k in original_hits]
        ax.scatter(db_indices, q_indices, s=1, c="b")

    grouped_hits = set()
    for group in hit_groups:
        for hit in group:
            grouped_hits.add(hit)  # using a set is faster that displaying overlapping hits
    db_indices = []
    q_indices = []
    for i, j, k in grouped_hits:
        db_indices.append(i + (k / 2))
        q_indices.append(j + (k / 2))
    ax.scatter(db_indices, q_indices, s=1, c="r")

    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    plt.show()


def main():
    description = "Finds all groups of N non-overlapping hits whose pairwise distances are less than A." \
                  "If an N is not provided then the default value is 2 (it looks for pairs)"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("hits_file", type=str, help="The file containing the hits in which we will find groups.")
    parser.add_argument("A", type=int,
                        help="The maximum possible distance between non-overlapping hits to consider them part of the"
                             "same group.")
    parser.add_argument("N", type=int, default=2, help="The number of hits to consider a group viable.")
    parser.add_argument("--display", "-d", action="store_true",
                        help="Displays the hit groups that are found in red and the original hits in blue.")
    args = parser.parse_args(sys.argv[1:])

    hits = []
    with open(args.hits_file, mode="r") as f:
        for line in f:
            if line[0] != "#":
                i, j, k = line.split()
                hits.append((int(i), int(j), int(k)))

    hit_groups = group_hits(hits, args.A, args.N)

    print(f"# Distance Limit (A): {args.A}")
    print(f"# Group Size (N): {args.N}")
    print(f"# {len(hit_groups)} Hit Groups:")
    for hit_group in hit_groups:
        hit_group = [f"{i} {j} {k} " for i, j, k in hit_group]
        print(",".join(hit_group))

    if args.display:
        display_hit_groups(hit_groups, hits)


if __name__ == "__main__":
    main()
