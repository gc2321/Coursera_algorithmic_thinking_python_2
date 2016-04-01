"""
Provide code and solution for Application 4
"""

DESKTOP = True

import math
import random
import urllib2
import matplotlib.pyplot as plt
import numpy as np
import string

if DESKTOP:
    import matplotlib.pyplot as plt
    import project4 as student
else:
    import simpleplot
    import userXX_XXXXXXX as student
    

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
    scoring_file = urllib2.urlopen(filename)
    ykeys = scoring_file.readline()
    ykeychars = ykeys.split()
    for line in scoring_file.readlines():
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
    protein_file = urllib2.urlopen(filename)
    protein_seq = protein_file.read()
    protein_seq = protein_seq.rstrip()
    return protein_seq


# Q1
#compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix)
seq_fly = read_protein(FRUITFLY_EYELESS_URL)
seq_human = read_protein(HUMAN_EYELESS_URL)
score_matrix = read_scoring_matrix(PAM50_URL)
alignment_matrix = student.compute_alignment_matrix(seq_human, seq_fly, score_matrix, False)

result = student.compute_local_alignment(seq_human, seq_fly, score_matrix, alignment_matrix)
#print result[0]
#human
#print result[1]
#fly
#print result[2]


# Q2
seq_pax = read_protein(CONSENSUS_PAX_URL)

#fly and pax domain
# alignment_matrix_global = student.compute_alignment_matrix(result[2], seq_pax, score_matrix, True)
# result2 = student.compute_global_alignment(result[2], seq_pax, score_matrix, alignment_matrix_global)
#print result2[0]
#print result2[1]
#print result2[2]

#human and pax domain
# seq_human_dashless =''
# for each in result[1]:
#     if each !='-':
#         seq_human_dashless +=each

#print seq_human_dashless
# alignment_matrix_global2 = student.compute_alignment_matrix(seq_human_dashless, seq_pax, score_matrix, True)
# result2b = student.compute_global_alignment(seq_human_dashless, seq_pax, score_matrix, alignment_matrix_global2)
#print result2b[0]
#print result2b[1]
#print result2b[2]


# Q4
def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    """
    generate null distribution of amino acid at specific position
    :param seq_x: seq_x
    :param seq_y: seq_y
    :param scoring_matrix:  scoring matrix
    :param num_trials: number of trials
    :return: a dictionary of scoring_distribution
    """
    scoring_distr= {}
    for i in xrange(1, num_trials+1):
        # random seq from seq_y
        rand_y = ''.join(random.sample(seq_y, len(seq_y)))

        alignment_matrix = student.compute_alignment_matrix(seq_x, rand_y, scoring_matrix, False)
        result = student.compute_local_alignment(seq_x, rand_y, scoring_matrix, alignment_matrix)

        scoring_distr[i]= result[0]

    return scoring_distr

#distr = generate_null_distribution(seq_human, seq_fly, score_matrix, 1000)


#plot graph
# x = []
# for each in distr:
#     x.append(distr[each])

# fig = plt.figure()
# ax = fig.add_subplot(111)
# n, bins, rectangles = ax.hist(x, 1000, normed=True)
# fig.canvas.draw()
# plt.title('Distribution of score of local alignment of\nHumanEyeless and random sequence of FlyEyeless')
# plt.ylabel('Distribution')
# plt.xlabel('Score')
# plt.show()


# Q5
# mean = np.mean(x)
# sd = np.std(x)
# z_score = (result[0]-mean)/sd
# print "mean: "+ str(mean)
# print "standard deviation: "+str(sd)
# print "z-score: "+str(z_score)

# Q8

def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = urllib2.urlopen(filename)

    # read in files as string
    words = word_file.read()

    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print "Loaded a dictionary with", len(word_list), "words"
    return word_list

words = read_words(WORD_LIST_URL)

def check_spelling(check_word, dist, word_list):
    """
    check spelling of check_word
    :param check_word: word to check
    :param dist: edit distance
    :param word_list: list of wrod (dictionary)
    :return: set of words from word_list that has the distance of 'dist' from check_word
    """
    result =[]
    alphabet = list(string.ascii_lowercase)
    score_matrix = student.build_scoring_matrix(alphabet, 2, 1, 0)

    for each in word_list:
        alignment_matrix = student.compute_alignment_matrix(each, check_word, score_matrix, True)
        global_align = student.compute_global_alignment(each, check_word, score_matrix, alignment_matrix)
        distance = len(each)+len(check_word)-global_align[0]
        if distance <= dist:
            result.append(each)

    return result

result = check_spelling('humble', 1, words)

for each in result:
    print each

print len(result)