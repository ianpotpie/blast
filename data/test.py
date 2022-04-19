import sys

with open("./test.txt") as f:
    print(f.readline())
    for line in f:
        print(line)