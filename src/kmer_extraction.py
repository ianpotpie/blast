from scoring_scheme import ScoringScheme
import argparse
import sys

DNA_SYMBOLS = ["A", "C", "G", "T"]
RNA_SYMBOLS = ["A", "C", "G", "T", "U"]
PROTEIN_SYMBOLS = ["A", "G", "I", "L", "P", "V", "F", "W", "Y", "D",
                   "E", "R", "H", "K", "S", "T", "C", "M", "N", "Q"]
ALPHABET = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
            "V", "W", "X", "Y", "Z"]


def get_kmers(seq, k):
    """
    Retrieves all k-mers from a query sequence and returns them as a list of (index, kmer) tuples.

    :param seq: the sequence from which seeds are extracted
    :param k: the size of the seeds/kmers
    :return: a list of seeds, (kmer,index) tuples
    """
    kmers = []
    for i in range((len(seq) + 1) - k):
        kmer = seq[i:i + k]
        kmers.append((i, kmer))
    return kmers


def build_neighborhood(seq, symbols, scoring_scheme, T):
    """
    Gets all the kmers that have a match score of more than T with the sequence. I call this the "prefix pruning"
    algorithm, since it builds up kmers from common prefixes whilst discarding inviable candidates.

    :param seq: the sequence whose neighborhood we are finding
    :param symbols: the symbols available to construct the sequence
    :param scoring_scheme: the scoring scheme use for sequence comparisons
    :param T: the threshold score for potential matches.
    :return: a list of sequences with an alignment score greater than T with the provided sequence
    """
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
    return [kmer for kmer, score in prefixes]


def get_neighborhoods(seq, symbols, k, scoring_scheme, T):
    """
    Gets all the k-mers that have a high-scoring subsequence within the sequence and returns them as a list of
    (index, kmer) tuples.

    :param seq: the sequence from which seeds are extracted
    :param symbols: the symbols available for the sequence
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

    kmers = []
    for subseq, indices in kmer_to_indices.items():
        neighborhood = build_neighborhood(subseq, symbols, scoring_scheme, T)
        for index in indices:
            for kmer in neighborhood:
                kmers.append((index, kmer))
    return kmers


def main():
    description = "Provided a source sequence and an integer k, this script finds all subsequences of the source " \
                  "sequence with length k. If a threshold score 'T' is provided, then it will find all kmers with " \
                  "gapless alignments scoring over T against some subsequence in the source sequence. The script " \
                  "prints a list of all such kmers and their positions in the source sequence. The user may provide " \
                  "their own scoring criteria by setting the match/mismatch scores or loading a scoring matrix file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("sequence", type=str,
                        help="The file containing the sequence from which kmers will be extracted."
                             "Characters on multiple lines are concatenated together (assumed to be the same sequence.")
    parser.add_argument("k", type=int,
                        help="The length of kmers/seeds\' to be extracted.")
    parser.add_argument("--threshold", "-t", type=float,
                        help="The threshold score for each kmer/seed to be added to the neighborhood of a subsequence.")
    parser.add_argument("--match", "-m", type=float, default=1,
                        help="The score of a match within an alignment (used if matrix is not provided).")
    parser.add_argument("--mismatch", "-e", type=float, default=-1,
                        help="The score of a mismatch within an alignment (used if matrix is not provided).")
    parser.add_argument("--matrix", "-s", type=str,
                        help="The file containing a scoring matrix used to score alignments.")
    parser.add_argument("--type", type=str, choices=["DNA", "RNA", "PROTEIN", "ALPHABET"], default="DNA",
                        help="The type of sequence from which we are extracting kmers (used to determine the possible"
                             "replacement symbols available).")
    args = parser.parse_args(sys.argv[1:])

    with open(args.sequence, mode="r") as f:
        sequence = ""
        for line in f:
            sequence += line.strip()

    if args.type == "RNA":
        symbols = RNA_SYMBOLS
    elif args.type == "PROTEIN":
        symbols = PROTEIN_SYMBOLS
    elif args.type == "ALPHABET":
        symbols = ALPHABET
    else:
        symbols = DNA_SYMBOLS

    if args.threshold is None:
        kmers = get_kmers(sequence, args.k)
    else:
        scoring_scheme = ScoringScheme(match=args.match, mismatch=args.mismatch)
        if args.matrix:
            scoring_scheme.load_matrix(args.matrix)
        kmers = get_neighborhoods(sequence, symbols, args.k, scoring_scheme, args.threshold)

    print(f"# Sequence: {sequence}")
    print(f"# Length (k): {args.k}")
    if args.threshold is not None:
        print(f"# Threshold (T): {args.threshold}")
    print(f"# {len(kmers)} Kmers:")
    for index, kmer in kmers:
        print(f"{index} {kmer}")


if __name__ == "__main__":
    main()
