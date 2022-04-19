import sys
from scoring_scheme import ScoringScheme


def extend_seed(sequence, query, sequence_index, query_index, k, scoring_scheme, falloff_limit):
    """
    Takes in the position of a seed of size k and returns the largest viable extension based on a scoring scheme and the
    falloff limit for the random walk score. The extension is provided in the form of its positions in both the
    search sequence and the query along with the length of the extension.

    :param sequence: the sequence in which the seed was matched
    :param query: the query from which the seed was derived
    :param sequence_index: the index of the seed in the sequence
    :param query_index: the index of the seed in the query
    :param k: the size of the seed match
    :param scoring_scheme: a function for scoring symbols
    :param falloff_limit: the limit for the score falloff when performing the extension
    :return: the left and right bounds of the extension of the seed
    """

    sequence_kmer = sequence[sequence_index: sequence_index + k]
    query_kmer = query[query_index: query_index + k]
    seed_score = scoring_scheme(sequence_kmer, query_kmer)

    # extend the query/database sequence to the right
    score, max_score = seed_score, seed_score
    right_bound = k
    max_index = min(len(sequence) - sequence_index, len(query) - query_index)
    for index in range(k, max_index):
        if max_score - score > falloff_limit:
            break
        next_sequence_symbol = sequence[sequence_index + index]
        next_query_symbol = query[query_index + index]
        score += scoring_scheme(next_sequence_symbol, next_query_symbol)
        if score > max_score:
            max_score = score
            right_bound = index + 1

    # extend the query/database sequence to the left
    score, max_score = seed_score, seed_score
    left_bound = 0
    max_index = min(sequence_index, query_index)
    for index in range(1, max_index + 1):
        if max_score - score > falloff_limit:
            break
        next_sequence_symbol = sequence[sequence_index - index]
        next_query_symbol = query[query_index - index]
        score += scoring_scheme(next_sequence_symbol, next_query_symbol)
        if score > max_score:
            max_score = score
            left_bound = -index

    sequence_index += left_bound
    query_index += left_bound
    extension_size = right_bound - left_bound
    return sequence_index, query_index, extension_size


def get_high_scoring_pairs(database, query, scoring_scheme, k, seeds, falloff_limit, cutoff_limit):
    """
    This takes in a database, a query sequence, and the locations of seeds found between the two. It extends the seeds
    and finds any high-scoring alignment pairs between the database and the query sequences.

    :param database: a list of sequences (strings)
    :param query: a sequence which we are trying to locate in the database (or a high-scoring match)
    :param scoring_scheme: the function used to score matches of the query in the database
    :param k: the size of the seeds whose matches we are provided
    :param seeds: locations of regions similar to k-sized subsequences of our query
    :param falloff_limit: the most the score of an extension can fall before termination
    :param cutoff_limit: the minimum score of a successful extension
    :return: a list of high-scoring alignment pairs
    """

    high_scoring_pairs = set()
    for seed in seeds:
        database_index, sequence_index, query_index = seed
        sequence = database[database_index]
        sequence_index, query_index, extension_size = extend_seed(sequence, query, sequence_index, query_index,
                                                                  k, scoring_scheme, falloff_limit)

        sequence_alignment = sequence[sequence_index: sequence_index + extension_size]
        query_alignment = query[query_index: query_index + extension_size]
        alignment_score = scoring_scheme(sequence_alignment, query_alignment)

        # the formatting of each tuple is necessary for sorting
        if alignment_score > cutoff_limit:
            high_scoring_pairs.add((alignment_score,
                                    len(query_alignment),
                                    -database_index,
                                    -sequence_index,
                                    -query_index,
                                    sequence_alignment,
                                    query_alignment))

    # we sort by score, then size, then indices
    high_scoring_pairs = list(high_scoring_pairs)
    high_scoring_pairs.sort(reverse=True)

    # after sorting, we reformat each tuple
    n_high_scoring_pairs = len(high_scoring_pairs)
    for i in range(n_high_scoring_pairs):
        high_scoring_pair = high_scoring_pairs[i]
        alignment_score = high_scoring_pair[0]
        sequence_alignment = high_scoring_pair[5]
        query_alignment = high_scoring_pair[6]
        database_index = -1 * high_scoring_pair[2]
        sequence_index = -1 * high_scoring_pair[3]
        query_index = -1 * high_scoring_pair[4]
        high_scoring_pairs[i] = (database_index, sequence_index, query_index,
                                 query_alignment, sequence_alignment, alignment_score)

    return high_scoring_pairs


def main():
    db_file = sys.argv[1]
    with open(db_file) as f:
        database = [line.strip() for line in f.readlines()]

    matrix_file = sys.argv[2]
    scoring_scheme = ScoringScheme()
    scoring_scheme.load_matrix(matrix_file)

    seeds = []
    seeds_file = sys.argv[3]
    with open(seeds_file) as f:
        query = f.readline().strip()
        k = int(f.readline())
        T = float(f.readline())
        n_seeds = int(f.readline())
        for line in f.readlines():
            line = line.split()
            seed = (int(line[1]), int(line[3]), int(line[5]))
            seeds.append(seed)

    falloff_limit = float(sys.argv[4])

    cutoff_limit = float(sys.argv[5])

    high_scoring_pairs = get_high_scoring_pairs(database, query, scoring_scheme, k, seeds, falloff_limit, cutoff_limit)

    print(query)
    print(k)
    print(int(T))
    print(len(high_scoring_pairs))
    print(35 * "-")
    for i_seq, i_pos, i_query, query_alignment, seq_alignment, score in high_scoring_pairs:
        print(f"Sequence {i_seq} Position {i_pos} Q-index {i_query}")
        print(query_alignment)
        print(seq_alignment)
        print(int(score))
        print(35 * "-")


if __name__ == "__main__":
    main()
