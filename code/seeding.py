import sys
import itertools
from scoring_scheme import ScoringScheme
from knuth_morris_pratt import kmp_search


def get_neighborhoods(sequence, scoring_scheme, k, T):
    """
    Gets all the k-mers that have a high-scoring subsequence within the sequence.

    :param sequence: the query string
    :param scoring_scheme: the system for scoring alignments
    :param k: the size of the k-mers
    :param T: the threshold score for seeds
    :return: a dictionary of high-scoring kmers to their query indices
    """

    symbols = scoring_scheme.get_symbols()
    indices_by_kmer = {}
    for kmer in itertools.product(symbols, repeat=k):
        kmer = "".join(kmer)
        indices_by_kmer[kmer] = []
        for i in range(len(sequence) - k):
            subsequence = sequence[i: i + k]
            align_score = scoring_scheme.score(kmer, subsequence)
            if align_score > T:
                indices_by_kmer[kmer].append(i)
        if len(indices_by_kmer[kmer]) == 0:
            del indices_by_kmer[kmer]
    return indices_by_kmer


def seeding(database, query, scoring_scheme, k, T):
    """
    Given a database and a set of seeds, finds all exact matches of those seeds in the database.

    :param database: a database of strings
    :param query: a query sequence
    :param scoring_scheme: an alignment scorer
    :param k: the length of the kmers of the query
    :param T: the minimum score of an alignment in the query to be considered a hit
    :return: a dictionary of queries to their matches (a tuple of the sequence index and position in that sequence)
    """

    indices_by_kmer = get_neighborhoods(query, scoring_scheme, k, T)

    seeds = []
    for nth_seq, sequence in enumerate(database):
        for kmer, query_indices in indices_by_kmer.items():
            hits = kmp_search(sequence, kmer)
            for seq_index in hits:
                for query_index in query_indices:
                    seeds.append((nth_seq, seq_index, query_index))
    seeds.sort()

    return seeds


def main():
    database_file = sys.argv[1]
    with open(database_file) as f:
        database = [line.strip() for line in f.readlines()]

    query_file = sys.argv[2]
    with open(query_file) as f:
        query = f.readline().strip()

    matrix_file = sys.argv[3]
    scoring_scheme = ScoringScheme()
    scoring_scheme.load_matrix(matrix_file)

    k = int(sys.argv[4])

    T = float(sys.argv[5])

    seeds = seeding(database, query, scoring_scheme, k, t)

    print(query)
    print(k)
    print(int(T))
    print(len(seeds))
    for nth_seq, i_pos, query_index in seeds:
        print(f"Sequence {i_seq} Position {i_pos} Q-index {i_query}")


if __name__ == "__main__":
    main()
