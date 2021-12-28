"""
This code is provided by instructors from Rice University
by online Algorithmic Thinking course through online learning platform Coursera.
"""

import math
import random
from urllib.request import urlopen

# URLs for data files
PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
HUMAN_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"


###############################################
# provided code

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = urlopen(filename)
    ykeys = scoring_file.readline().decode('utf-8')
    ykeychars = ykeys.split()
    lines = [i.decode("utf-8") for i in scoring_file.readlines()]
    for line in lines:
        vals = line.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][ykey] = int(val)
    return scoring_dict


def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = urlopen(filename)
    protein_seq = protein_file.read().decode('utf-8')
    protein_seq = protein_seq.rstrip()
    return protein_seq


def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = urlopen(filename)

    # read in files as string
    words = word_file.read()

    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print
    "Loaded a dictionary with", len(word_list), "words"
    return word_list

def delete_dashes(string):
    """ Delete all dashes from a given string"""
    dashless_str = ""
    for char in string:
        if char == "-":
            continue
        else:
            dashless_str += char
    return dashless_str


def calculate_seqs_agreement(seq1, seq2):
    """
    Given two same length sequences, this function calculates the percentage of the agreeing element over whole sequences
    """
    assert len(seq1) == len(seq2)
    seq_length = len(seq1)
    agr_number = 0
    percentage = 0

    for key_i in range(seq_length):
        if seq1[key_i] == seq2[key_i]:
            agr_number += 1
    if agr_number:
        percentage = agr_number/float(seq_length)*100
    return percentage