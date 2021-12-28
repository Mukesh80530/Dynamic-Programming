import random
import numpy as np
import helper_func
import matplotlib.pyplot as plt

# URLs for data files
PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
HUMAN_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"


class SequenceAlignment:
    """
    for global alignments by Needleman–Wunsch algorithm & for local alignments by Smith–Waterman algorithm
    using Dynamic Programming Functions
    """

    def __init__(self):
        pass

    def build_scoring_matrix(self, alphabet, diag_score, off_diag_score, dash_score):
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

    def compute_alignment_matrix(self, seq_x, seq_y, scoring_matrix, global_flag):
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
                alignment_matrix[key_i].append(alignment_matrix[key_i - 1][0] + scoring_matrix[seq_x[key_i - 1]]["-"])

            for key_j in range(1, lenght_y + 1):  # row initialization
                alignment_matrix[0].append(alignment_matrix[0][key_j - 1] + scoring_matrix["-"][seq_y[key_j - 1]])

            for key_i in range(1, lenght_x + 1):  # iteration
                for key_j in range(1, lenght_y + 1):
                    alignment_matrix[key_i]. \
                        append(
                        max(alignment_matrix[key_i - 1][key_j - 1] + scoring_matrix[seq_x[key_i - 1]][seq_y[key_j - 1]],
                            alignment_matrix[key_i - 1][key_j] + scoring_matrix[seq_x[key_i - 1]]["-"],
                            alignment_matrix[key_i][key_j - 1] + scoring_matrix["-"][seq_y[key_j - 1]]))
            return alignment_matrix

        else:  # Case for Smith-Waterman: Local Alignments
            for key_i in range(1,
                               lenght_x + 1):  # column initialization assumend that gaps create negative values: set to 0
                alignment_matrix[key_i].append(0)
            for dummy_j in range(1,
                                 lenght_y + 1):  # row initialization assumend that gaps create negative values: set to 0
                alignment_matrix[0].append(0)

            for key_i in range(1, lenght_x + 1):  # iteration
                for key_j in range(1, lenght_y + 1):
                    alignment_matrix[key_i].append(max(0, alignment_matrix[key_i - 1][key_j - 1] +
                                                       scoring_matrix[seq_x[key_i - 1]][seq_y[key_j - 1]],
                                                       alignment_matrix[key_i - 1][key_j] +
                                                       scoring_matrix[seq_x[key_i - 1]]["-"],
                                                       alignment_matrix[key_i][key_j - 1] + scoring_matrix["-"][
                                                           seq_y[key_j - 1]]))
            return alignment_matrix

    def compute_global_alignment(self, seq_x, seq_y, scoring_matrix, alignment_matrix):
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

        while lenght_x != 0 and lenght_y != 0:
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
        while lenght_x != 0:
            x_val = seq_x[lenght_x - 1]
            x_fin = x_val + x_fin
            y_fin = "-" + y_fin
            lenght_x -= 1

        while lenght_y != 0:
            y_val = seq_y[lenght_y - 1]
            x_fin = "-" + x_fin
            y_fin = y_val + y_fin
            lenght_y -= 1

        return (final_score, x_fin, y_fin)

    def compute_local_alignment(self, seq_x, seq_y, scoring_matrix, alignment_matrix):
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

    def compare_seq_cons(self, seq1, consensus, matrix):
        """
        Function that takes sequences and consensus sequence
        Remove their dashs, align them globally and measure the similarity percentage.
        :param seg1:
        :param consenus:
        :return: their similarity percentage
        """
        seq1_dashless = helper_func.delete_dashes(seq1)
        align_m = self.compute_alignment_matrix(seq1_dashless, consensus, matrix, True)
        alignments = self.compute_global_alignment(seq1_dashless, consensus, matrix, align_m)
        percentage = helper_func.calculate_seqs_agreement(alignments[1], alignments[2])
        return percentage

    def generate_null_distribution(self, seq_x, seq_y, scoring_matix, num_trials):
        """

        :param seq_x:
        :param seq_y:
        :param scoring_matix:
        :param num_trials:
        :return:
        """
        result_dict = {}
        for trial_num in range(num_trials):
            print("Number count--->", trial_num)
            seq_y_copy = (seq_y + ".")[:-1]

            seq_y_copy = list(seq_y_copy)
            random.shuffle(seq_y_copy)
            ".".join(seq_y_copy)

            align_m = self.compute_alignment_matrix(seq_x, seq_y_copy, scoring_matix, False)
            trail_result = self.compute_local_alignment(seq_x, seq_y_copy, scoring_matix, align_m)
            if trail_result[0] in result_dict.keys():
                result_dict[trail_result[0]] += 1
            else:
                result_dict[trail_result[0]] = 1
        return result_dict

    def normalize_distribution(self, dist):
        """
        normalize the entered distribution given as dictionary
        :param dist:
        :return:
        """
        total = 0
        dist_copy = dist.copy()
        for key in dist:
            total += dist[key]

        for key in dist:
            dist_copy[key] = dist[key] / float(total)

        return dist_copy


if __name__ == '__main__':
    sequence_alignment_obj = SequenceAlignment()

    human_gene = helper_func.read_protein(HUMAN_EYELESS_URL)
    fruitfly_gene = helper_func.read_protein(FRUITFLY_EYELESS_URL)
    pam_50_matix = helper_func.read_scoring_matrix(PAM50_URL)
    consensus_protein = helper_func.read_protein(CONSENSUS_PAX_URL)

    print("Human Eyeless Gene ", human_gene)
    print("Fruitfly Eyeless Gene ", fruitfly_gene)
    print("consensus_protein ", consensus_protein)

    alignment_matrix = sequence_alignment_obj.compute_alignment_matrix(human_gene, fruitfly_gene, pam_50_matix, False)
    alignment_loc = sequence_alignment_obj.compute_local_alignment(human_gene, fruitfly_gene, pam_50_matix,
                                                                   alignment_matrix)

    for element in alignment_loc:
        print(element)

    human_cons_perc = sequence_alignment_obj.compare_seq_cons(alignment_loc[1], consensus_protein, pam_50_matix)
    fly_cons_perc = sequence_alignment_obj.compare_seq_cons(alignment_loc[2], consensus_protein, pam_50_matix)

    print("Human-Consensus Agreement Percentage: ", human_cons_perc)
    print("Fly-Consensus Agreement Percentage: ", fly_cons_perc)

    distribution_unnormalized = sequence_alignment_obj.generate_null_distribution(human_gene, fruitfly_gene, pam_50_matix, 100)
    distribution_normalized = sequence_alignment_obj.normalize_distribution(distribution_unnormalized)

    # Creating the plot
    keys_normalized = distribution_normalized.keys()
    values_normalized = [round(val, 3) for val in distribution_normalized.values()]

    fig = plt.figure()
    ax = fig.add_axes([0,0,2,1])
    ax.bar(keys_normalized, values_normalized, color="blue")
    ax.set_xticks(np.arange(keys_normalized[0], keys_normalized[-1], 2))
    ax.set_ylabel("Dist, fraction of total trials")
    ax.set_ylabel("Scores")
    ax.set_title("Score by normalized dist. fraction of the experiment")
    plt.show()
