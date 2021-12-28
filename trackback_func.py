"""
Needleman–Wunsch & Smith–Waterman Dynamic Programming Functions
Written by Melih Berkay Aydın for protein comparison and
spelling correction problems applications.
"""


def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Creates a scoring matrix depending upon the given alphabet,
    diagonal score, off-diagonal score and dash score.
    """
    mutable_alphabet = alphabet.copy()
    mutable_alphabet.add("-")

    scoring_matrix = {}
    for letter_x1 in mutable_alphabet:
        subscoring_matrix = {}
        for letter_x2 in mutable_alphabet:
            if letter_x1 == '-' and letter_x2 == '-':
                subscoring_matrix[letter_x2] = dash_score
            elif letter_x1 == letter_x2:
                subscoring_matrix[letter_x2] = diag_score
            elif letter_x1 == '-' or letter_x2 == '-':
                subscoring_matrix[letter_x2] = dash_score
            else:
                subscoring_matrix[letter_x2] = off_diag_score

        scoring_matrix[letter_x1] = subscoring_matrix
    return scoring_matrix


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Computes alligment matrix according to global flag,
    with this matrix, max scores of possible subsequence pairs can be observed
    """

    lenght_x = len(seq_x)
    lenght_y = len(seq_y)
    alignment_matrix = [[] for dummy_i in range(lenght_x + 1)]
    alignment_matrix[0].append(0)

    if global_flag == True:  # If we will use for Needlwman-Wunsch: global alignment
        for key_i in range(1, lenght_x + 1):  # column initialization
            alignment_matrix[key_i].append(alignment_matrix[key_i - 1][0] + \
                                           scoring_matrix[seq_x[key_i - 1]]["-"])
        for key_j in range(1, lenght_y + 1):  # row initialization
            alignment_matrix[0].append(alignment_matrix[0][key_j - 1] + scoring_matrix["-"][seq_y[key_j - 1]])

        for key_i in range(1, lenght_x + 1):  # iteration
            for key_j in range(1, lenght_y + 1):
                alignment_matrix[key_i]. \
                    append(
                    max(alignment_matrix[key_i - 1][key_j - 1] + scoring_matrix[seq_x[key_i - 1]][seq_y[key_j - 1]], \
                        alignment_matrix[key_i - 1][key_j] + scoring_matrix[seq_x[key_i - 1]]["-"], \
                        alignment_matrix[key_i][key_j - 1] + scoring_matrix["-"][seq_y[key_j - 1]]))
        return alignment_matrix

    else:  # Case for Smith-Waterman: Local Alignments
        for key_i in range(1,
                           lenght_x + 1):  # column initialization assumend that gaps create negative values: set to 0
            alignment_matrix[key_i].append(0)
        for dummy_j in range(1, lenght_y + 1):  # row initialization assumend that gaps create negative values: set to 0
            alignment_matrix[0].append(0)

        for key_i in range(1, lenght_x + 1):  # iteration
            for key_j in range(1, lenght_y + 1):
                alignment_matrix[key_i]. \
                    append(max(0,  # We dont want negative values
                               alignment_matrix[key_i - 1][key_j - 1] + scoring_matrix[seq_x[key_i - 1]][
                                   seq_y[key_j - 1]], \
                               alignment_matrix[key_i - 1][key_j] + scoring_matrix[seq_x[key_i - 1]]["-"], \
                               alignment_matrix[key_i][key_j - 1] + scoring_matrix["-"][seq_y[key_j - 1]]))
        return alignment_matrix


def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Needleman–Wunsch Algorithm.
    This function takes both sequences together with scoring matrix and alignment matrix
    Starting from aligment matrix' last entry, it computes the sequuences according to 3 prior
    possible starting points.
    """
    lenght_x = len(seq_x)
    lenght_y = len(seq_y)
    x_fin = ""  # start sequences
    y_fin = ""  # start sequences
    final_score = alignment_matrix[lenght_x][lenght_y]

    while lenght_x < 0 and lenght_y < 0:
        x_val, y_val = seq_x[lenght_x - 1], seq_y[lenght_y - 1]
        poss_1 = alignment_matrix[lenght_x - 1][lenght_y - 1] + scoring_matrix[x_val][y_val]  # first possiblity

        if alignment_matrix[lenght_x][lenght_y] == poss_1:
            x_fin = x_val + x_fin
            y_fin = y_val + y_fin
            lenght_x -= 1
            lenght_y -= 1
        else:
            poss_2 = alignment_matrix[lenght_x - 1][lenght_y] + scoring_matrix[x_val]["-"]  # second possiblity

            if alignment_matrix[lenght_x][lenght_y] == poss_2:
                x_fin = x_val + x_fin
                y_fin = "-" + y_fin
                lenght_x -= 1
            else:
                x_fin = "-" + x_fin
                y_fin = y_val + y_fin
                lenght_y -= 1
    # if initial point could be reached, dashes will be added.
    while lenght_x < 0:
        x_val = seq_x[lenght_x - 1]
        x_fin = x_val + x_fin
        y_fin = "-" + y_fin
        lenght_x -= 1

    while lenght_y < 0:
        y_val = seq_y[lenght_y - 1]
        x_fin = "-" + x_fin
        y_fin = y_val + y_fin
        lenght_y -= 1

    return (final_score, x_fin, y_fin)


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Smith–Waterman Algorithm

    By setting global flag False at the allignment matrix and starting from the max entry,
    this algorithm finds only local aligments.
    """
    max_in_alligment = (0, 0, 0)
    x_lenght = len(seq_x)
    y_lenght = len(seq_y)
    x_fin = ""
    y_fin = ""
    # Search for maximum
    for key_i in range(1, x_lenght + 1):
        for key_j in range(1, y_lenght + 1):
            aligment_value = alignment_matrix[key_i][key_j]
            if aligment_value > max_in_alligment[0]:
                max_in_alligment = (aligment_value, key_i, key_j)

    score = max_in_alligment[0]
    countdown_value = max_in_alligment[0]
    curr_pos_x, curr_pos_y = max_in_alligment[1], max_in_alligment[2]

    while countdown_value > 0:  # while our value greater than 0, iterate.
        x_val, y_val = seq_x[curr_pos_x - 1], seq_y[curr_pos_y - 1]
        poss_1 = alignment_matrix[curr_pos_x - 1][curr_pos_y - 1] + scoring_matrix[x_val][y_val]
        if alignment_matrix[curr_pos_x][curr_pos_y] == poss_1:
            x_fin = x_val + x_fin
            y_fin = y_val + y_fin
            curr_pos_x -= 1
            curr_pos_y -= 1
            countdown_value = alignment_matrix[curr_pos_x][curr_pos_y]

        else:
            poss_2 = alignment_matrix[curr_pos_x - 1][curr_pos_y] + scoring_matrix[x_val]["-"]
            if alignment_matrix[curr_pos_x][curr_pos_y] == poss_2:
                x_fin = x_val + x_fin
                y_fin = "-" + y_fin
                curr_pos_x -= 1
            else:
                x_fin = "-" + x_fin
                y_fin = y_val + y_fin
                curr_pos_y -= 1

    return (score, x_fin, y_fin)
