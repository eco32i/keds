__author__ = 'Mandel & Shamovsky'
"""Takes the following arguments: The barcode length, the number of barcodes, 
the number of errors to correct (only 1 or 2), minimum GC content(int, ie, 40),
max GC content(also int)"""

import random
import numpy as np
import argparse
import time


# Initialize lists
nucl_list = ['A', 'C', 'G', 'T']
barcode_list = []
tested = []
tests = 0
# We need to initialize this list to keep track of tested barcodes
# including duplicates in order to terminate program if too much time elapses


def gc_cont(barcode, length):
    """Returns the GC content of a barcode"""
    gc = 0.0
    for base in range(length):
        if barcode[base] == 'C' or barcode[base] == 'G':
            gc += 1
        else:
            gc += 0
    cont = gc / length
    return cont


def make_barcode():
    """Generates a random barcode from nucl_list"""
    barcode = ''
    while barcode == '':
        for i in range(length):
            barcode += random.choice(nucl_list)
        if maxgc >= gc_cont(barcode, length) >= mingc:
            bar_code = barcode
        else:
            barcode = ''
    return bar_code


def SLdistance(s1, s2):
    """Calculates the hamming distance between s1 and s2"""

    # Initiate np array
    matrix = np.zeros((len(s1) + 1, len(s2) + 1), dtype=np.int)
    matrix[:, 0] = np.array([i for i in xrange(len(s1) + 1)])
    matrix[0, :] = np.array([i for i in xrange(len(s2) + 1)])

    # // Classical Levenshtein part
    for i in xrange(1, len(s1) + 1):
        for j in xrange(1, len(s2) + 1):
            cost = 0
            if s1[i - 1] != s2[j - 1]:
                cost = 1
            matrix[i, j] = min(matrix[(i - 1), (j - 1)] + cost,
                               matrix[i, (j - 1)] + 1,
                               matrix[(i - 1), j] + 1)
    min_distance = matrix[len(s1)][len(s2)]

    # New Sequence-Levenshtein part

    # Truncating
    for i in xrange(0, len(s1) + 1):
        min_distance = min(min_distance, matrix[i, len(s2)])
    # Elongating
    for j in xrange(0, len(s2) + 1):
        min_distance = min(min_distance, matrix[len(s1), j])
    return min_distance


def compare_distances(new_barcode, num_errors):
    """Compares the sequence-Levenshtein distance between
    new barcode and old barcodes
    Uses the S-L distance depending on # errors
    to correct (2 * k + 1) k = errors"""
    # Count number of barcodes with 'bad' distances
    count = 0
    global barcode_list
    if num_errors == 1:
        for barcode in barcode_list:
            if SLdistance(new_barcode, barcode) < 3:
                count += 1
    elif num_errors == 2:
        for barcode in barcode_list:
            if SLdistance(new_barcode, barcode) < 5:
                count += 1
    return count


def complement(barcode):
    """returns the complement of the barcode"""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = ''
    for base in barcode:
        complement += complement_dict[base]
    return complement


def compare_complements(new_barcode):
    """Returns a count > 0 if generated barcode is a complement of
    any in current list"""
    complement_count = 0
    global barcode_list
    for barcode in barcode_list:
        if complement(barcode) == complement(new_barcode):
            complement_count += 1
    return complement_count


def compare_repeat(barcode, length):
    """Returns a count > 0 if 2 consecutive bases in a barcode are the same"""
    count = 0
    for i in range(length - 1):
        if barcode[i] == barcode[i + 1]:
            count += 1
    return count


def compare_barcodes(num_errors):
    """Main, monster function...
    Compares a barcode list, which can correct up to 'num_errors'.
    Also does an ongoing comparison of new generated barcodes and checks for:
    1. The desired S-L distance between each barcode
    2. Excludes self-complements
    3. Excludes any barcodes that contain 2 duplicate consecutive bases"""
    global tests
    new_barcode = make_barcode()
    tests += 1
    if new_barcode not in tested:
        tested.append(new_barcode)
    if new_barcode not in barcode_list:
        distance_count = compare_distances(new_barcode, num_errors)
        complement_count = compare_complements(new_barcode)
        repeat_count = compare_repeat(new_barcode, length)
        if distance_count > 0 or complement_count > 0 or repeat_count > 0:
            pass
        else:
            barcode_list.append(new_barcode)
    else:
        pass

class BarcodeGenerator(object):
    length = 6
    number_barcodes = 10
    errors = 1
    mingc = 45.
    maxgc = 55.

    def __init__(self, kwargs):
# Not sure this is a right way to do it
        for k,v in kwargs.items:
            self[k] = v


    def generate(self):
        
        now = time.time()
        future = now + 10
        while len(barcode_list) < self.number_barcodes:
            if time.time() >= future:
                print "This is taking too much time...statistics show that if you \
    run the program again,\nyou may get quicker results."
                print "Didn't create the desired number of barcodes... \
    only generated {0}. Better luck next time!".format(len(barcode_list))
                break
            compare_barcodes(errors)
            if tests >= 200000:
                print "Generated ONLY {0} barcodes before termination".format(len(barcode_list))
                break
        barcode_list.sort()
        if len(barcode_list) == 10:
            print "Correcting up to {0} error(s)".format(errors)
            print "SUCCESS"
            if 1 == errors:
                min_dist = 3
            elif 2 == errors:
                min_dist = 5
            print "Created {0} barcodes of length {1}, with a S-L distance of at least {2} " \
                "and a gc content range between {3}% and {4}%.\n{5}" \
                .format(len(barcode_list), length, min_dist, args.mingc, args.maxgc, barcode_list)



def main():
    # generate(.....)
    

    parser = argparse.ArgumentParser()
    parser.add_argument('length', nargs='?', default=6, type=int,
                    help='takes the number nucleotides for each barcode')
    parser.add_argument('number_barcodes', nargs='?', default=10, type=int,
                    help='takes the number of barcodes to generate')
    parser.add_argument('errors', nargs='?', default=1, type=int, choices=[1, 2],
                    help='gives the number of mismatches')
    parser.add_argument('mingc', nargs='?', default=45.0, type=float,
                    help='enter the minimum gc count')
    parser.add_argument('maxgc', nargs='?', default=55.0, type=float,
                    help='enter the maximum gc count')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    mingc = args.mingc / 100
    maxgc = args.maxgc / 100
    number_barcodes = args.number_barcodes
    errors = args.errors
    length = args.length
    if args.verbose:
        print 'Took the following arguments:\nlength: {0}\nnumber_barcodes: {1}\nerrors: {2}\n\
            mingc: {3}\nmaxgc: {4}'.format(args.length, args.number_barcodes, args.errors, args.mingc,
                               args.maxgc)



    generator = BarcodeGenerator(args)
    generator.generate()


if __name__ == '__main__':
    main()




