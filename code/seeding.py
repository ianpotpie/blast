import sys
import itertools
import heapq
from scoring_scheme import ScoringScheme
from knuth_morris_pratt import kmp_search


def get_seeds(seq: str, k: int):
    """
    Retrieves all k-mers from a query sequence and returns them as a list of (index, kmer) tuples.

    :param seq: the sequence from which seeds are extracted
    :param k: the size of the seeds/kmers
    :return: a list of seeds, (kmer,index) tuples
    """
    seeds = []
    for i in range((len(seq) + 1) - k):
        kmer = seq[i:i + k]
        seeds.append((i, kmer))
    return seeds


def get_neighborhoods(seq: str, k: int, scoring_scheme: ScoringScheme, T: float):
    """
    Gets all the k-mers that have a high-scoring subsequence within the sequence and returns them as a list of
    (index, kmer) tuples.

    :param seq: the sequence from which seeds are extracted
    :param k: the size of the k-mers/seeds
    :param scoring_scheme: the system for scoring alignments
    :param T: the threshold score for seeds
    :return: a list of seeds, (kmer, index) tuples
    """
    seeds = []
    symbols = scoring_scheme.get_symbols()
    for i in range((len(seq) + 1) - k):
        subseq = seq[i: i + k]
        for kmer in itertools.product(symbols, repeat=k):
            kmer = "".join(kmer)
            if scoring_scheme(kmer, subseq) > T:
                seeds.append((i, kmer))
    return seeds


def find_hits(db_seq: str, seeds: list):
    """
    Finds instances of the seeds in the database sequence.

    :param db_seq: the database sequence
    :param seeds: a list of seeds (query_index, kmer)
    :return: a list of hits (db_index, query_index, length)
    """
    kmer_to_indices = {kmer: [] for i, kmer in seeds}  # dictionary prevents searching for the same kmer twice
    for query_index, kmer in seeds:
        kmer_to_indices[kmer].append(query_index)

    hits = []
    for seed, query_indices in kmer_to_indices.items():
        k = len(seed)
        db_indices = kmp_search(db_seq, seed)
        for db_index in db_indices:
            for query_index in query_indices:
                hits.append((db_index, query_index, k))

    return hits


def combine_hits(hits, A):
    """
    Combines all hits with ends that either overlap or have ends within a distance of A on the diagonal.

    :param hits: a list of seed-hit tuples (db_index, query_index, length)
    :param A: a distance threshold for merging hits
    :return: a new dictionary of hits to their
    """
    hits_by_diagonal = {}
    for db_index, query_index, k in hits:
        diagonal = db_index - query_index  # this defines the diagonal
        if diagonal not in hits_by_diagonal:
            hits_by_diagonal[diagonal] = []
        heapq.heappush(hits_by_diagonal[diagonal], (db_index, query_index, k))  # hits will be sequential

    combined_hits = []
    for diagonal, aligned_hits in hits_by_diagonal.items():
        start_i, start_j, combined_k = aligned_hits[0]
        prev_i, prev_j = start_i, start_j
        for i, j, k in aligned_hits:
            if i - prev_i <= A:
                combined_k = (i - start_i) + k
            else:
                combined_hits.append((start_i, start_j, combined_k))
                start_i, start_j, combined_k = i, j, k
            prev_i, prev_j = i, j

    return combined_hits


def extend_hits(hits, X, S):
    pass


def main():
    # the database file should contain sequences seperated by "\n"
    db_seq = sys.argv[1]

    # the query file should contain a single query string
    query_seq = sys.argv[2]

    # the matrix file contains the scoring matrix
    matrix_file = sys.argv[3]
    scoring_scheme = ScoringScheme()
    scoring_scheme.load_matrix(matrix_file)

    # this is the size of the seeds/kmers
    k = int(sys.argv[4])

    # the minimum alignment score (we are guaranteed alignments over this threshold)
    T = float(sys.argv[5])

    seeds = seeding(database, query, scoring_scheme, k, T)

    print(query)
    print(k)
    print(int(T))
    print(len(seeds))
    for nth_seq, seq_index, query_index in seeds:
        print(f"Sequence {nth_seq} Position {seq_index} Q-index {query_index}")


if __name__ == "__main__":
    main()
