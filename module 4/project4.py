"""
project 4
"""

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    build a scoring matrix for sequence alignments
    :param alphabet: a set of characters
    :param diag_score: scores A-X, x!=A, etc...
    :param off_diag_score: scores for A-A, T-T, etc...
    :param dash_score: score for any entry indexed by one or more dashes
    :return: dictionary of dictionaries, eg.'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2},
    """
    alist = list(alphabet)
    alist.append('-')

    dict1 = {}
    for val1 in alist:
        dict2 = {}
        for val2 in alist:
            if (val1 is '-') or (val2 is '-'):
                dict2[val2] = dash_score
            elif val1 is val2:
                dict2[val2] = diag_score
            else:
                dict2[val2] = off_diag_score

        dict1[val1]= dict2

    return dict1

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    compute alignment matrix of two sequences
    :param seq_x: one sequence
    :param seq_y: another sequence
    :param scoring_matrix: dictionary of scoring matrix
    :param global_flag: either global alignment or local alignment
    :return: alignment matrix
    """

    matrix = [[0]*(len(seq_y)+1) for _ in range(len(seq_x)+1)]

    def check_global(value):
        """
        check if global_flag is true
        """
        if value < 0 and global_flag==False:
            return 0
        else:
            return value

    # global alignment
    for idx in range(1, len(seq_x)+1):
        matrix [idx][0] += matrix [idx-1][0]
        matrix [idx][0] += scoring_matrix[seq_x[idx-1]]['-']
        matrix [idx][0] = check_global(matrix[idx][0])

    for idx in range(1, len(seq_y)+1):
        matrix [0][idx] += matrix [0][idx-1]
        matrix [0][idx] += scoring_matrix['-'][seq_y[idx-1]]
        matrix [0][idx] = check_global(matrix [0][idx])

    for idx1 in range(1, len(seq_x)+1):
        for idx2 in range(1, len(seq_y)+1):
            value1 = matrix[idx1-1][idx2-1]+scoring_matrix[seq_x[idx1-1]][seq_y[idx2-1]]
            value2 = matrix[idx1-1][idx2]+scoring_matrix[seq_x[idx1-1]]['-']
            value3 = matrix[idx1][idx2-1]+scoring_matrix['-'][seq_y[idx2-1]]

            #print idx1, idx2, seq_x[idx1-1], seq_y[idx2-1]
            matrix[idx1][idx2] += max(value1, value2, value3)
            matrix[idx1][idx2] = check_global(matrix[idx1][idx2])

    return matrix

#def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
#score_matrix = build_scoring_matrix(set(['A', 'C', 'T', 'G']), 5, 2, -2)
#print score_matrix
#matrix1 = compute_alignment_matrix('AC', 'TAG', score_matrix, True)

# Q12 homework
#score_matrix = build_scoring_matrix(set(['A', 'C', 'T', 'G']), 10, 4, -6)
#matrix1 = compute_alignment_matrix('AA', 'TAAT', score_matrix, False)

#for row in matrix1:
#    print row

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    computer alignment of two sequences
    :param seq_x: one sequence
    :param seq_y: another sequence
    :param scoring_matrix: scoring matrix
    :param alignment_matrix: alignment matrix
    :return:
    """
    idx_i = len(seq_x)
    idx_j = len(seq_y)

    if idx_i <1 and idx_j <1:
        return (0, seq_x, seq_y)

    align_x = ""
    align_y = ""

    while idx_i !=0 and idx_j !=0:
        if alignment_matrix[idx_i][idx_j]==alignment_matrix[idx_i-1][idx_j-1]+scoring_matrix[seq_x[idx_i-1]][seq_y[idx_j-1]]:
            align_x = seq_x[idx_i-1]+align_x
            align_y = seq_y[idx_j-1]+align_y
            idx_i -=1
            idx_j -=1
        else:
            if alignment_matrix[idx_i][idx_j]==alignment_matrix[idx_i-1][idx_j]+scoring_matrix[seq_x[idx_i-1]]['-']:
                align_x = seq_x[idx_i-1]+align_x
                align_y = '-'+align_y
                idx_i -=1
            else:
                align_x = '-'+align_x
                align_y = seq_y[idx_j-1]+align_y
                idx_j -=1

    while idx_i !=0:
        align_x = seq_x[idx_i-1]+align_x
        align_y = '-'+align_y
        idx_i -=1

    while idx_j !=0:
        align_x = '-'+align_x
        align_y = seq_y[idx_j-1]+align_y
        idx_j -=1

    align_matrix_final = compute_alignment_matrix(align_x, align_y, scoring_matrix, True)
    score = align_matrix_final[len(align_x)][len(align_y)]

    return (score, align_x, align_y)

#print compute_global_alignment('AA', 'TAAT', score_matrix, matrix1)

def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    computer alignment of two sequences
    :param seq_x: one sequence
    :param seq_y: another sequence
    :param scoring_matrix: scoring matrix
    :param alignment_matrix: alignment matrix
    :return:
    """
    if len(seq_x) <1 and len(seq_y) <1:
        return (0, seq_x, seq_y)

    align_x = ""
    align_y = ""

    #obtain max
    max_value = 0
    for idx1 in reversed(range (len(alignment_matrix))):
        for idx2 in reversed(range(len(alignment_matrix[idx1]))):
            #print idx1, idx2
            if alignment_matrix[idx1][idx2] > max_value:
                max_value = alignment_matrix[idx1][idx2]
                location = (idx1, idx2)

    idx_i = location[0]
    idx_j = location[1]

    while idx_i !=0 and idx_j !=0:
        if alignment_matrix[idx_i][idx_j]==alignment_matrix[idx_i-1][idx_j-1]+scoring_matrix[seq_x[idx_i-1]][seq_y[idx_j-1]]:
            align_x = seq_x[idx_i-1]+align_x
            align_y = seq_y[idx_j-1]+align_y
            idx_i -=1
            idx_j -=1
        else:
            if alignment_matrix[idx_i][idx_j]==alignment_matrix[idx_i-1][idx_j]+scoring_matrix[seq_x[idx_i-1]]['-']:
                align_x = seq_x[idx_i-1]+align_x
                align_y = '-'+align_y
                idx_i -=1
            else:
                align_x = '-'+align_x
                align_y = seq_y[idx_j-1]+align_y
                idx_j -=1

    #print align_x, align_y

    ax_list = list(align_x)
    ay_list = list(align_y)

    while ax_list[0]=='-' or ay_list[0]=="-":
        ax_list.pop(0)
        ay_list.pop(0)

    while ax_list[-1]=='-' or ay_list[-1]=="-":
        ax_list.pop()
        ay_list.pop()

    align_x = ''.join(ax_list)
    align_y = ''.join(ay_list)

    align_matrix_final = compute_alignment_matrix(align_x, align_y, scoring_matrix, False)

    return (align_matrix_final[len(align_x)][len(align_y)], align_x, align_y)

