#!/usr/bin/python

"""

Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

"""

import os
import sys

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        """
        Input: (String) A path to a scoring matrix file.
        Output: (Dictionary) A nested dictionary where keys are strings
                and elements are scores as integers.

        Example:

        >>> matrix, gap_penalty = NWAligner.load_score_matrix('/home/bioe131/BLOSUM62')
        >>> matrix['A']['A']
        4
        >>> matrix['W']['W']
        11
        >>> gap_penalty
        -4

        """

        score_matrix = {}
        gap_penalty = None

        with open(fname) as fp:
            temp=[]
            for line_num, line in enumerate(fp):
                # ignore comments in matrix file
                if line.startswith("#"):
                    continue
                else:
                    temp.append(line)

        gap_penalty=int(temp[-1].split('\n')[0])
        matrix=temp[:-1]

        order = matrix[0].split('\n')[0]
        name = []
        for i in order.split(' '):
            name.append(i)
            score_matrix[i] = {}

        rest = matrix[1:]

        for i in range(len(rest)):
            val = rest[i].split('\n')[0]
            for s in range(len(val.split(' '))):
                score_matrix[name[i]][name[s]]=int(val.split(' ')[s])

        return score_matrix, gap_penalty

    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """
        #seq_x='YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        #seq_y='WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        seqs = []

        with open(fname, 'rt') as fp:
            temp=[]
            for i in fp:
                temp.append(i)
            sequences = ''.join(temp).split('>')
            if len(sequences) != 3:
                raise error
            else:
                for i in sequences[1:]:
                    seqs.append(''.join(i.split('\n')[1:]))

        return seqs

    def align(self, seq_x, seq_y, print_matrix = False):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        for i in range(len(matrix)):
            if i == 0:
                for y in range(len(matrix)):
                    if y != 0:
                        matrix[i][y]=matrix[i][y-1] + self.gap_penalty
            else:
                matrix[i][0] = matrix[0][i]
        for i in range(len(pointers)):
            if i == 0:
                for y in range(len(pointers)):
                    if pointers[i][y] == 0:
                        pointers[i][y] = 'left'
            else:
                pointers[i][0] = 'above'

        for x in range(1, len(seq_x) + 1):
            for y in range(1, len(seq_y) + 1):
                match_score = self.score_matrix[seq_x[x - 1]][seq_y[y - 1]]
                left = matrix[x-1][y] + self.gap_penalty
                above = matrix[x][y-1] + self.gap_penalty
                diagonal = matrix[x-1][y-1] + match_score
                options = [left,above,diagonal]
                #print('the options are ' + str(left) + ' ' + str(above) + ' ' + str(diagonal) + ' and the max is ' + str(max(options)))
                if max(options) == left:
                    #print('left')
                    #matrix[x][y] = left
                    pointers[x][y] = 'left'

                elif max(options) == above:
                    #print('above')
                    #matrix[x][y] = above
                    pointers[x][y] = 'above'

                elif max(options) == diagonal:
                    #print('diagonal')
                    #matrix[x][y] = diagonal
                    pointers[x][y] = 'diagonal'
                    #print(match_score)
                matrix[x][y] = max(options)

        # print the dynamic programming matrix
        if print_matrix:
            for x in range(len(seq_x) + 1):
                print " ".join(map(lambda i: str(int(i)), matrix[x]))
        #print(matrix[-1][-1])
        #print(matrix)
        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x), len(seq_y)

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []

        while x > 0 or y > 0:
            #break
            move = pointers[x][y]
            #print(move)
            if move == 'left':
                #print('left')
                align_y.append('-')
                align_x.append(seq_x[x-1])
                x, y = x-1, y
            elif move == 'above':
                #print('above')
                align_y.append(seq_y[y-1])
                align_x.append('-')
                x, y = x, y-1
            elif move == 'diagonal':
                #print('diagonal')
                align_y.append(seq_y[y-1])
                align_x.append(seq_x[x-1])
                x, y = x-1, y-1

        # flip the alignments, as they're reversed
        return ("".join(align_x[::-1]), "".join(align_y[::-1]))

###                                      ###
### NO NEED TO EDIT CODE BELOW THIS LINE ###
###                                      ###

if __name__ == '__main__':
    def usage():
        print 'usage: %s matrixfilename stringfilename'
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print 'Can not open %s' % (fname,)
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))
