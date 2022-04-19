import sys
from scoring_scheme import ScoringScheme


def main():
    database_file = sys.argv[1]
    with open(database_file) as f:
        database = []
        for line in f:
            database.append(line.strip())

    query_file = sys.argv[2]
    with open(query_file) as f:
        queries = []
        for line in f:
            queries.append(line.strip())

    score_matrix_file = sys.argv[3]
    scoring_scheme = ScoringScheme()
    scoring_scheme.load_matrix(score_matrix_file)

    k = int(sys.argv[4])

    T = float(sys.argv[5])


if __name__ == "__main__":
    main()
