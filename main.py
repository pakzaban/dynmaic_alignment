# coding=utf-8
"""
Dynamic Programming project for DNA sequence alignment
"""


def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Builds a scoring matrix for comparing DNA nucleotides
    :param alphabet_string: string, list, tuple or set representing the DNA alphabet (nucleotides)
    :param diag_score: score if the nucleotides match
    :param off_diag_score: score if the nucleotides mismatch
    :param dash_score: score if one is a dash
    :return: scoring matrix as dictionary of dictionaries
    """

    # turn alphabet from whatever iterable data type into a string
    alphabet_string = ""
    for item in alphabet:
        alphabet_string += item
    alphabet_string += "-"

    #build matrix as dictionary of dictionaries
    matrix_size = len(alphabet_string)
    scoring_matrix = {}
    for row in range(matrix_size):
        scoring_matrix[alphabet_string[row]] = {}
        for col in range(matrix_size):
            if alphabet_string[row] == alphabet_string[col]:
                scoring_matrix[alphabet_string[row]][alphabet_string[row]] = diag_score
            else:
                scoring_matrix[alphabet_string[row]][alphabet_string[col]] = off_diag_score
            if alphabet_string[row] == "-" or alphabet_string[col] == "-":
                scoring_matrix[alphabet_string[row]][alphabet_string[col]] = dash_score
    return scoring_matrix

def compute_alignment_matrix(seq_x,seq_y,scoring_matrix,global_flag):
    """
    :param seq_x:
    :param seq_y:
    :param scoring_matrix:
    :param global_flag:
    :return:
    """
    length_x = len(seq_x) + 1
    length_y = len(seq_y) + 1
    s_matrix = [[0 for dummy_col in range(length_y)] for dummy_row in range(length_x)]
    if global_flag:
        for i_index in range(1, length_x):
            s_matrix[i_index][0] = s_matrix[i_index-1][0] + scoring_matrix[seq_x[i_index-1]]["-"]
        for j_index in range(1, length_y):
            s_matrix[0][j_index] = s_matrix[0][j_index-1] + scoring_matrix["-"][seq_y[j_index-1]]
        for i_index in range(1, length_x):
            for j_index in range(1, length_y):
                first = s_matrix[i_index-1][j_index-1] + scoring_matrix[seq_x[i_index-1]][seq_y[j_index-1]]
                second = s_matrix[i_index-1][j_index] + scoring_matrix[seq_x[i_index-1]]["-"]
                third = s_matrix[i_index][j_index-1] + scoring_matrix["-"][seq_y[j_index-1]]
                s_matrix[i_index][j_index] = max(first,second,third)
    else:
        for i_index in range(1, length_x):
            s_matrix[i_index][0] = max(s_matrix[i_index - 1][0] + scoring_matrix[seq_x[i_index - 1]]["-"],0)
        for j_index in range(1, length_y):
            s_matrix[0][j_index] = max(s_matrix[0][j_index - 1] + scoring_matrix["-"][seq_y[j_index - 1]],0)
        for i_index in range(1, length_x):
            for j_index in range(1, length_y):
                first = max(s_matrix[i_index - 1][j_index - 1] + scoring_matrix[seq_x[i_index - 1]][seq_y[j_index - 1]],0)
                second = max(s_matrix[i_index - 1][j_index] + scoring_matrix[seq_x[i_index - 1]]["-"],0)
                third = max(s_matrix[i_index][j_index - 1] + scoring_matrix["-"][seq_y[j_index - 1]],0)
                s_matrix[i_index][j_index] = max(first, second, third)

    return s_matrix


def compute_global_alignment(seq_x,seq_y,scoring_matrix,alignment_matrix):
    """
    This function computes a global alignment of seq_x and seq_y using the global alignment matrix alignment_matrix.The
    function returns a tuple of the form (score,align_x,align_y) where score is the score of the global alignment align_x
    and align_y. Note that align_x and align_y should have the same length and may include the padding character ’-’.
    :param seq_x: string
    :param seq_y: string
    :param scoring_matrix: dict of dict from build_scoring_matrix()
    :param alignment_matrix: list of lists from compute_alignment_matrix()
    :return: tuple (score,align_x,align_y)
    """
    i_index = len(seq_x)
    j_index = len(seq_y)
    score = alignment_matrix[i_index][j_index]
    align_x = ""
    align_y = ""
    while i_index != 0 and j_index != 0:
        if alignment_matrix[i_index][j_index] == alignment_matrix[i_index-1][j_index-1] + scoring_matrix[seq_x[i_index-1]][seq_y[j_index-1]]:
            align_x = seq_x[i_index-1] + align_x
            align_y = seq_y[j_index-1] + align_y
            i_index -= 1
            j_index -= 1
        else:
            if alignment_matrix[i_index][j_index] == alignment_matrix[i_index-1][j_index] + scoring_matrix[seq_x[i_index-1]]["-"]:
                align_x = seq_x[i_index-1] + align_x
                align_y = "-" + align_y
                i_index -= 1
            else:
                align_x = "-" + align_x
                align_y = seq_y[j_index-1] + align_y
                j_index -= 1
    while i_index != 0:
        align_x = seq_x[i_index-1] + align_x
        align_y = "-" + align_y
        i_index -= 1
    while j_index != 0:
        align_x = "-" + align_x
        align_y = seq_y[j_index - 1] + align_y
        j_index -= 1

    return (score, align_x, align_y)

def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    This function computes a local alignment of seq_x and seq_y using the local alignment matrix, alignment_matrix.The
    function returns a tuple of the form (score,align_x,align_y) where score is the score of the local alignment align_x
    and align_y. Note that align_x and align_y should have the same length and may include the padding character ’-’.
    :param seq_x: string
    :param seq_y: string
    :param scoring_matrix: dict of dict from build_scoring_matrix()
    :param alignment_matrix: list of lists from compute_alignment_matrix() with global_flag = False
    :return: tuple (score,align_x,align_y)
    """
    score = 0
    i_index = 0
    j_index = 0
    align_x = ""
    align_y = ""
    for row in range(len(alignment_matrix)):
        for col in range(len(alignment_matrix[0])):
            if alignment_matrix[row][col] > score:
                score = alignment_matrix[row][col]
                i_index = row
                j_index = col

    while i_index != 0 and j_index != 0:
        # print "i=", i_index, " j=", j_index, " value=", alignment_matrix[i_index][j_index]
        if alignment_matrix[i_index][j_index]==0:
            return (score,align_x,align_y)
        if alignment_matrix[i_index][j_index] == alignment_matrix[i_index-1][j_index-1] + scoring_matrix[seq_x[i_index-1]][seq_y[j_index-1]]:
            align_x = seq_x[i_index-1] + align_x
            align_y = seq_y[j_index-1] + align_y
            i_index -= 1
            j_index -= 1
        else:
            if alignment_matrix[i_index][j_index] == alignment_matrix[i_index-1][j_index] + scoring_matrix[seq_x[i_index-1]]["-"]:
                align_x = seq_x[i_index-1] + align_x
                align_y = "-" + align_y
                i_index -= 1
            else:
                align_x = "-" + align_x
                align_y = seq_y[j_index-1] + align_y
                j_index -= 1

    while i_index != 0:
        # print "i=", i_index, " j=", j_index, " value=", alignment_matrix[i_index][j_index]
        if alignment_matrix[i_index][j_index]==0:
            return (score,align_x,align_y)
        align_x = seq_x[i_index-1] + align_x
        align_y = "-" + align_y
        i_index -= 1

    while j_index != 0:
        # print "i=", i_index, " j=", j_index, " value=", alignment_matrix[i_index][j_index]
        if alignment_matrix[i_index][j_index]==0:
            return (score,align_x,align_y)
        align_x = "-" + align_x
        align_y = seq_y[j_index - 1] + align_y
        j_index -= 1

    return (score, align_x, align_y)

def print_matrix(matrix):
    """
    prints a matrix that is a 2d-list
    :param matrix: list of list
    :return: None
    """
    height = len(matrix)
    width = len(matrix[0])
    for row in range(height):
        for col in range(width):
            print matrix[row][col],"|",
        print "\n----------------"

# M = build_scoring_matrix(("A","T","C","G"), 10, 2, -4)
# print M

# print_matrix(compute_alignment_matrix("AA", "TAAT",M,True))
# alm = compute_alignment_matrix("ACC", "TTTACACGG",M, False)
# print_matrix(alm)

# print compute_global_alignment("AA","TAAT",M,alm)
# print compute_global_alignment('ATG', 'ACG',
#                                {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2},
#                                 'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
#                                 '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
#                                 'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
#                                 'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}},
#                                [[0, -4, -8, -12], [-4, 6, 2, -2], [-8, 2, 8, 4], [-12, -2, 4, 14]])


# print compute_local_alignment('ACTACT', 'AGCTA',
#                         {'A': {'A': 2, 'C': 1, '-': 0, 'T': 1, 'G': 1},
#                          'C': {'A': 1, 'C': 2, '-': 0, 'T': 1, 'G': 1},
#                          '-': {'A': 0, 'C': 0, '-': 0, 'T': 0, 'G': 0},
#                          'T': {'A': 1, 'C': 1, '-': 0, 'T': 2, 'G': 1},
#                          'G': {'A': 1, 'C': 1, '-': 0, 'T': 1, 'G': 2}},
#                         [[0, 0, 0, 0, 0, 0], [0, 2, 2, 2, 2, 2], [0, 2, 3, 4, 4, 4],
#                          [0, 2, 3, 4, 6, 6], [0, 2, 3, 4, 6, 8], [0, 2, 3, 5, 6, 8], [0, 2, 3, 5, 7, 8]])

# print compute_local_alignment("AA", "TAAT", M, alm)
# print compute_local_alignment("ACC", "TTTACACGG",M, alm)