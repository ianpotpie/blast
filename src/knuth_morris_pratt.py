import sys


def build_failure_function(sequence):
    """
    Takes in a sequence and returns an array where the ith index of the array contains the size of the largest proper
    prefix of the sequence that is also a proper suffix of the substring ending at index i.

    :param sequence: a sequence of values
    :return: the failure function of the sequence
    """

    # no proper prefix or suffix for size 1
    failure_function = [0]

    prefix_end = 0
    subsequence_end = 1
    while subsequence_end < len(sequence):
        prefix_symbol = sequence[prefix_end]
        suffix_symbol = sequence[subsequence_end]

        # if the values match, then we increment the failure function
        if prefix_symbol == suffix_symbol:
            prefix_end += 1
            failure_function.insert(subsequence_end, prefix_end)
            subsequence_end += 1

        # on a mismatch, decrement the prefix end until it hits 1
        else:
            if prefix_end != 0:
                prefix_end = failure_function[prefix_end - 1]
            else:
                failure_function.insert(subsequence_end, prefix_end)
                subsequence_end += 1

    return failure_function


def kmp_search(database, query):
    """
    Uses the KMP algorithm in order to search a database sequence for instance of the query sequence.

    :param database: a sequence
    :param query: a sequence
    :return: a list of indices indicating locations of query sequences in the database sequence
    """
    hits = []

    failure_function = build_failure_function(query)
    db_index = 0
    q_index = 0
    while db_index < len(database):

        # if the characters match then we can advance the indices
        if query[q_index] == database[db_index]:
            q_index += 1
            db_index += 1

        # if we reach the end of the query, then we have a match, and we revert the q_index with the failure function
        if q_index == len(query):
            hits.append(db_index - q_index)
            q_index = failure_function[q_index - 1]

        # mismatch after q_index matches
        elif db_index < len(database) and query[q_index] != database[db_index]:
            if q_index != 0:
                q_index = failure_function[q_index - 1]
            else:
                db_index += 1

    return hits


def main():
    database = sys.argv[1]

    query = sys.argv[2]

    hits = kmp_search(database, query)

    for hit in hits:
        print(hit)


if __name__ == "__main__":
    main()
