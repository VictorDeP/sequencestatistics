"""
DNA and protein probability calculators
Author: Victor De Pillecyn
Date: 21-01-2022
Last update: 21-01-2022
A series of functions to calculate DNA and protein probabilities
"""
from Bio.Align import substitution_matrices
from math import comb
from numpy import arange
import re

BLOSUM62 = substitution_matrices.load("BLOSUM62")
DNA_ALPHABET = ['A', 'C', 'G', 'T']
RNA_ALPHABET = ['A', 'C', 'G', 'U']
PROTEIN_ALPHABET = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']


def GCcontent(DNAseq):
    """
    Calculates GC content of DNA sequence
    Args:
        DNAseq (string): DNA sequence
    Returns:
        float: GC content in decimal notation (f.e. 0.67)
    """
    return len(re.findall('[GC]', DNAseq))/len(DNAseq)

# GCcontent('ACGT')

def transcribeDNA(DNAseq):
    """
    Transcribes DNA sequence to RNA sequence (T -> U)
    Args:
        DNAseq (string): DNA sequence
    Returns:
        string: RNA sequence
    """
    return DNAseq.replace('T','U')

#transcribeDNA('ACTG')

def residucount(seq):
    """
    Counts number of residus in DNA or PROTEIN sequence
    Args:
        seq (string): DNA or PROTEIN sequence
    Returns:
        dict: Counts of residus
    """
    alphabet = set(seq)
    counts = {}
    for r in seq:
        if r in counts:
            counts[r] += 1
        else:
            counts[r] = 1
    return counts

#residucount('AAACT')

def residufreq(seq):
    """
    Counts relative number (frequencies) of residus in DNA or PROTEIN sequence
    Args:
        seq (string): DNA or PROTEIN sequence
    Returns:
        dict: Frequencies of residus (residu count/len(sequence))
    """
    alphabet = set(seq)
    counts = {}
    frequencies= {}
    for r in seq:
        if r in counts:
            counts[r] += 1
        else:
            counts[r] = 1
    for key in counts.keys():
        frequencies[key] = counts[key]/len(seq)
    return frequencies

#residufreq('AAACT')


def seqprobDNA(DNAseq, similarity = 1, equiprob=dict(zip(DNA_ALPHABET, [1/len(DNA_ALPHABET)]*len(DNA_ALPHABET)))):
    """
    Calculates probability of two DNA sequences sharing certain similarity by chance, with certain equiprobability of bases
    Args:
        DNAseq (string): DNA sequence
        similarity (integer): Decimal similarity of query and sequence found by chance
        equiprob (dict): Dictionary with residu frequencies (can be obtained with residufreq(DNAseq))
    Returns:
        float: Probability of finding sequence with similarity by chance
    """
    if not re.search('[^ACTG]', DNAseq):
        totalprobability = comb(len(DNAseq), round(len(DNAseq)*similarity)) # Number of combinations to place dissimilar residues (=1 if similarity of 1)
        for key in equiprob.keys(): # For each residu
            residuprobability = equiprob[key] # Probability of residu occuring
            # DNAseq.count(key)*similarity = number of residu that have to match, vs. DNAseq.count(key)*1-similarity = number of residu that can not match
            totalprobability *= (residuprobability**(DNAseq.count(key)*similarity))*((1-residuprobability)**(DNAseq.count(key)*(1-similarity)))
    else:
        totalprobability = 'Not a DNA strand'
    return totalprobability

seqprobDNA('TGATTGACTATGCTTTACCG', 8/20, residufreq('CCCCCCGATAT'))
seqprobDNA('ATCTCT', similarity = 0.5)

 def minseqprobDNA(DNAseq, minsimilarity = 0.5, equiprob=dict(zip(DNA_ALPHABET, [1/len(DNA_ALPHABET)]*len(DNA_ALPHABET)))):
    """
    Calculates probability of two DNA sequences sharing AT LEAST certain similarity by chance, with certain equiprobability of bases
    Args:
        DNAseq (string): DNA sequence
        minsimilarity (integer): Decimal minimal similarity of two sequences
        equiprob (dict): Dictionary with residu frequencies (can be obtained with residufreq(DNAseq))
    Returns:
        float: Probability of finding sequence with minsimilarity or more by chance
    """
    totalprobability = 0
    if not re.search('[^ACTG]', DNAseq):
        min_residue_match = round(len(DNAseq)*minsimilarity)
        print(min_residue_match)
        for similarity in range(minsimilarity, 1, 0.05): #NOT WORKING
            totalprobability += seqprobDNA(DNAseq, similarity=1, equiprob)
    else:
        totalprobability = 'Not a DNA strand'
    return totalprobability

minseqprobDNA('ATCTCT')

# Identical functions but for proteins

def seqprobPROTEIN(PROTEINseq, similarity = 1, equiprob=dict(zip(PROTEIN_ALPHABET, [1/len(PROTEIN_ALPHABET)]*len(PROTEIN_ALPHABET)))):
    """
    Calculates probability of two PROTEIN sequences sharing certain similarity by chance, with certain equiprobability of bases
    Args:
        PROTEINseq (string): PROTEIN sequence
        similarity (integer): Decimal similarity of two sequences
        equiprob (dict): Dictionary with residu frequencies (can be obtained with residufreq(PROTEINseq))
    Returns:
        float: Probability of finding sequence with similarity by chance
    """
    if not re.search('[^ARNDCQEGHILKMFPSTWYV]', PROTEINseq):
        totalprobability = comb(len(PROTEINseq), round(len(PROTEINseq)*similarity)) # Number of combinations to place dissimilar residus (=1 if similarity of 1)
        for key in equiprob.keys(): # For each residu
            residuprobability = equiprob[key] # Probability of residu occuring
            # PROTEINseq.count(key)*similarity = number of residu that have to match, vs. PROTEINseq.count(key)*1-similarity = number of residu that can not match
            totalprobability *= (residuprobability**(PROTEINseq.count(key)*similarity))*((1-residuprobability)**(PROTEINseq.count(key)*(1-similarity)))
    else:
        totalprobability = 'Not a PROTEIN strand'
    return totalprobability

seqprobPROTEIN('MQAMRCTLVP', 0.5)

def minseqprobPROTEIN(PROTEINseq, minsimilarity = 0.5, equiprob=dict(zip(PROTEIN_ALPHABET, [1/len(PROTEIN_ALPHABET)]*len(PROTEIN_ALPHABET)))):
    """
    Calculates probability of two PROTEIN sequences sharing AT LEAST certain similarity by chance, with certain equiprobability of bases
    Args:
        PROTEINseq (string): PROTEIN sequence
        similarity (integer): Decimal similarity of two sequences
        equiprob (dict): Dictionary with residu frequencies (can be obtained with residufreq(PROTEINseq))
    Returns:
        float: Probability of finding sequence with minsimilarity or more by chance
    """
    totalprobability = 0
    if not re.search('[^ARNDCQEGHILKMFPSTWYV]', PROTEINseq):
        for similarity in arange(minsimilarity, 1, 0.05):
            totalprobability += seqprobPROTEIN(PROTEINseq, similarity, equiprob)
    else:
        totalprobability = 'Not a PROTEIN strand'
    return totalprobability

minseqprobPROTEIN('MQAMRCTLVP', 0.5)
