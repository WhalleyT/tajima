#!/usr/bin/env python3

import os
import subprocess
import argparse

from math import sqrt
from itertools import combinations

"""
This is a more data-agnostic implementation of Tajima's
D. It doesn't need any file format conversions (i.e. it
works on FASTA). There is a function to call MUSCLE, but
by default it accepts a raw FASTA
"""

def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-F", dest="fasta_file", type=str,
                        help="FASTA file of interest. Sequences must \
                        be of the same length.", required=True)
    parser.add_argument("--muscle", "-M", dest="muscle", action="store_true",
                        help="Flag. Run MUSCLE and use aligned output \
                        in calculation.")
    
    args = parser.parse_args()
    return(args)


def _run_msa(file):
    """Run clustal from executable."""
    subprocess.call(["muscle", "-in", file, "-out", "MUSCLE.out"])

def _calculate_pairwise(sequences):
    """Calculate pi, number of pairwise differences."""
    for seq in sequences:
        if len(seq) != len(sequences[0]):
            raise("All sequences must have the same length.")

    numseqs = len(sequences)

    num = float(numseqs * (numseqs - 1)) / float(2)

    combos = combinations(sequences, 2)
    counts = []
    for pair in combos:
        seqA = pair[0]
        seqB = pair[1]
        count = sum(1 for a, b in zip(seqA, seqB) if a != b)
        counts.append(count)

    return(float(sum(counts)) / float(num))


def _calculate_segregating_sites(sequences):
    """Calculate S, number of segregation sites)."""
    # Assume if we're in here seqs have already been checked
    combos = combinations(sequences, 2)
    indexes = []
    for pair in combos:
        seqA = pair[0]
        seqB = pair[1]
        for idx, (i, j) in enumerate(zip(seqA, seqB)):
            if i != j:
                indexes.append(idx)

    indexes = list(set(indexes))

    S = len(indexes)
    n = len(sequences)

    denom = 0
    for i in range(1, n):
        denom += (float(1) / float(i))

    return float(S / denom)


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def _read_sequences(file):
    sequences = []
    with open(file) as f:
        for name, seq in read_fasta(f):
            sequences.append(seq)
    return sequences


def _D(l, pi, s):
    a1 = sum([1.0/i for i in range(1, l)])
    a2 = sum([1.0/(i**2) for i in range(1, l)])
    
    b1 = float(l+1)/(3*(l-1))
    b2 = float(2 * ((l**2) + l + 3)) / (9*l*(l-1))
    
    c1 = b1 - 1.0/a1
    c2 = b2 - float(l+2)/(a1 * l) + float(a2)/(a1 ** 2)
    
    e1 = float(c1) / a1
    e2 = float(c2) / ( (a1**2) + a2 )
    
    D = (float(pi - (float(s)/a1)) /
         sqrt((e1 * s)+
         ((e2 * s) * (s - 1))))
    
    return D


def tajima():
    args = _parse_args()
    
    if args.muscle:
        _run_msa(args.fasta_file)
        sequences = _read_sequences("MUSCLE.out")
        os.remove("MUSCLE.out")
    else:
        sequences = _read_sequences(args.fasta_file)

    pi = _calculate_pairwise(sequences)
    S = _calculate_segregating_sites(sequences)

    """
    Now we have pi (pairwise differences) and s (number
    of segregating sites). This gives us 'little d', so
    now we need to divide it by sqrt of variance.
    """
    l = len(sequences)
    D = _D(l, pi, S)
    print("D = %f" %D)

if __name__ == "__main__":
    tajima()
