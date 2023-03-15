__author__ = 'Wombat'


class smithwaterman (object):
    """ smithwaterman.py takes a query string and a search string and uses
    a basic version (with no gaps) of the Smith-Waterman algorithm to find regions of alignment between the two.
    Version 2.0: 1 July 2015
    VERSION 2.0: 30 October 2019 -- conversion to Python 3.x
    """

    def __init__(self, query, target, match_reward=1, mismatch_penalty=-2):

        self.match_reward = match_reward    # This is the window size
        self.mismatch_penalty = mismatch_penalty
        self.query = query      # This is the string corresponding to the query sequence
        self.target = target    # This is the string corresponding to the target sequence

        self.querylen = len(query)
        self.targetlen = len(target)

        self.table = [[0 for i in range(self.targetlen)] for j in range(self.querylen)]
        # The above "list comprehension" in two dimensions will populate the Smith Waterman matrix with zeros

    def score(self):

        highscore = high_i = high_j = 0       # highest scores encountered so far in the matrix

        best_q_alignment = []  # best alignment for the query sequence
        best_t_alignment = []  # best alignment for the target sequence

        for i in range(self.querylen):

            for j in range(self.targetlen):

                queryword = self.query[i: i + 1]  # an array slice is perhaps more natural in python than a substring
                targetword = self.target[j: j + 1]

                if queryword == targetword:     # Did we have a match at this position?

                    increment = self.match_reward       # increment will contain either a positive reward or a
                    # negative penalty depending on whether we matched
                else:

                    increment = self.mismatch_penalty

                if i and j:  # We must test if we are at an edge to avoid looking back at non-existent table elements

                    self.table[i][j] = max(0, self.table[i - 1][j - 1] + increment)  # consider value on diagonal

                else:   # if we are at an edge, just record the match reward, as it must be the beginning

                    self.table[i][j] = max(0, increment)

                if self.table[i][j] > highscore:

                    highscore = self.table[i][j]    # record the new high score
                    high_i = i                      # also record the i and j positions associated with that score
                    high_j = j

        i = high_i
        j = high_j

        while self.table[i][j] and i > -1 and j > -1:

            best_q_alignment.append(self.query[i])
            best_t_alignment.append(self.target[j])

            i -= 1
            j -= 1

        best_q_alignment.reverse()
        best_t_alignment.reverse()

        return '\nBest alignment had a score of ' + str(highscore) + ' and is:\n\nTarget:\t' + \
               str(j + 2) + '\t' + ''.join(best_t_alignment) + '\n\t\t\t' + \
               '|' * len(best_t_alignment) + '\n' + \
               'Query:\t' + \
               str(i + 2) + '\t' + ''.join(best_q_alignment) + '\n'

        # The above long statement just concatenates together a multiline string that will correctly display
        # the best alignment when it is subsequently printed. Displaying the correct alignment is greatly eased
        # by the fact that there are no gaps -- we'll need a different approach later!

    def __str__(self):
        """ This is a "special method attribute" that controls what the string representation of objects
        of the class will look like it is invoked by print statements.
        """

        lineout = '\t' + '\t'.join(self.target) + '\n'
        # The above is just a fancy looking way to break the target string into tab delimited individual characters

        for i in range(self.querylen):
            lineout = lineout + self.query[i] + '\t'
            for j in range(self.targetlen):

                lineout = lineout + str(self.table[i][j]) + '\t'

            lineout += '\n'

        return lineout

###  MAIN PROGRAM

A = smithwaterman('AAGCGCGTCGCGATTGACATTATAGAAGGCGC', 'ATGTTCGCGATTGACATTCAAAAGCGAGTCAG', 1, -2)

print (A.score())

print (A)