import random
import itertools
from openpyxl import Workbook
from pathlib import Path


def random_seq_creator(length):

    bases = ["A", "C", "G", "T"]
    sequence = []

    for _ in range(length):
        base = random.choices(bases, weights=(29, 21, 21, 29), k=1)
        sequence = sequence + (base)

    sequence = ''.join(sequence)
    return sequence


def consecutive_g_count(sequence):

    consecutive_G = 0

    counter = ([[k, len(list(g))] for k, g in itertools.groupby(sequence)])
    #G_counter = [char for char in counter if counter [0]=="G"]
    for char in counter:
        if char[0] == "G":
            G_count = char[1]
            if G_count >= consecutive_G:
                consecutive_G = G_count

    return consecutive_G


def consecutive_c_count(sequence):

    consecutive_C = 0

    counter = ([[k, len(list(g))] for k, g in itertools.groupby(sequence)])
    #C_counter = [char for char in counter if counter [0]=="G"]
    for char in counter:
        if char[0] == "C":
            C_count = char[1]
            if C_count >= consecutive_C:
                consecutive_C = C_count

    return consecutive_C


def gc_content(sequence, length):

    GC_count = 0

    for base in sequence:
        if base == 'G' or base == 'C':
            GC_count = GC_count + 1

    GC_percentage = GC_count/length * 100

    return GC_percentage


def sequence_creator(length):

    should_restart = True

    while should_restart == True:

        should_restart = False
        sequence = random_seq_creator(length)
        consecutive_G = consecutive_g_count(sequence)
        consecutive_C = consecutive_c_count(sequence)
        gc_percentage = gc_content(sequence, length)

        if consecutive_C > 4 or consecutive_G > 4 or gc_percentage > 44:
            should_restart = True
            continue

    #print ('Finished,', "Consecutive C = ", consecutive_C, "Consecutive G = ", consecutive_G)

    return sequence, gc_percentage
