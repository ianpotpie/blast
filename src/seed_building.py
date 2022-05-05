from scoring_scheme import ScoringScheme
from itertools import product
import argparse
import sys


def get_kmers(seq, k):
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


def get_neighborhood(seq, scoring_scheme, T):
    """
    Gets all the kmers that have a match score of more than T with the sequence. I call this the "prefix pruning"
    algorithm, since it builds up kmers from common prefixes whilst discarding inviable candidates.

    :param seq: the sequence whose neighborhood we are finding
    :param scoring_scheme: the scoring scheme use for sequence comparisons
    :param T: the threshold score for potential matches.
    :return: a list of sequences with an alignment score greater than T with the provided sequence
    """
    symbols = scoring_scheme.get_symbols()
    prefixes = [("", 0.0)]
    for i in range(len(seq)):
        suffix = seq[i + 1:]
        best_suffix_score = scoring_scheme(suffix, suffix)
        next_prefixes = []
        for prefix, prefix_score in prefixes:
            for symbol in symbols:
                symbol_score = scoring_scheme(seq[i], symbol)
                if prefix_score + symbol_score + best_suffix_score > T:  # this prevents extending unusable prefixes
                    next_prefixes.append((prefix + symbol, prefix_score + symbol_score))
        prefixes = next_prefixes
    return [seq for seq, score in prefixes]


def get_seeds(seq: str, k: int, scoring_scheme: ScoringScheme, T: float):
    """
    Gets all the k-mers that have a high-scoring subsequence within the sequence and returns them as a list of
    (index, kmer) tuples.

    :param seq: the sequence from which seeds are extracted
    :param k: the size of the k-mers/seeds
    :param scoring_scheme: the system for scoring alignments
    :param T: the threshold score for seeds
    :return: a list of seeds, (index, kmer) tuples
    """
    kmer_to_indices = {}  # prevents the need for checking a kmer multiple times
    for i in range((len(seq) + 1) - k):
        kmer = seq[i: i + k]
        if kmer not in kmer_to_indices:
            kmer_to_indices[kmer] = []
        kmer_to_indices[kmer].append(i)

    seeds = []
    for subseq, indices in kmer_to_indices.items():
        neighborhood = get_neighborhood(subseq, scoring_scheme, T)
        for index in indices:
            for kmer in neighborhood:
                seeds.append((index, kmer))
    return seeds


def main():
    parser_description = "Provided a source sequence and an integer k, this script finds all subsequences of the " \
                         "source sequence with length k. If a threshold score 'T' is provided, then it will find all " \
                         "kmers with gapless alignments scoring over T against some subsequence in the source " \
                         "sequence. The script prints a list of all such kmers and their positions in the source " \
                         "sequence. The user may provide their own scoring criteria by setting the match/mismatch " \
                         "scores or loading a scoring matrix file."
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument("seq", type=str)
    parser.add_argument("k", type=int)
    parser.add_argument("--threshold", "-t", type=float)
    parser.add_argument("--match", "-m", type=float, default=1)
    parser.add_argument("--mismatch", "-e", type=float, default=-1)
    parser.add_argument("--matrix", "-s", type=str)
    args = parser.parse_args(sys.argv[1:])

    if args.threshold is None:
        seeds = get_kmers(args.seq, args.k)
    else:
        scoring_scheme = ScoringScheme(match=args.match, mismatch=args.mismatch)
        if args.matrix:
            scoring_scheme.load_matrix(args.matrix)
        seeds = get_seeds(args.seq, args.k, scoring_scheme, args.threshold)

    print(f"# Sequence: {args.seq}")
    print(f"# Length (k): {args.k}")
    if args.threshold is not None:
        print(f"# Threshold (T): {args.threshold}")
    print(f"# {len(seeds)} Seeds:")
    for index, kmer in seeds:
        print(f"{index} {kmer}")


if __name__ == "__main__":
    main()
