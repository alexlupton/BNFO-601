_author__ = 'Wombat'


class smithwaterman:
    """  smithwaterman2.py takes a query string and a search string and uses
        a version the Smith-Waterman algorithm supporting gaps to find regions of alignment between the two.
        Version 1.0: 1 July 2015
        VERSION 2.0: 30 October 2019 -- conversion to Python 3.x
    """

    def __init__(self, query, target, match_reward=1, mismatch_penalty=-2, gap_penalty=-5):

        self.match_reward = match_reward    # This is the window size
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

        self.query = query      # This is the string corresponding to the query sequence
        self.target = target    # This is the string corresponding to the target sequence

        self.querylen = len(query)
        self.targetlen = len(target)

        self.table = [[0 for i in range(self.targetlen)] for j in range(self.querylen)]
        # The above "list comprehension" in two dimensions will populate the Smith Waterman matrix with zeros

        self.traceback = [[(0, 0) for i in range(self.targetlen)] for j in range(self.querylen)]
        # Again, a list comprehension that initializes the traceback table with (0,0) tuples

        # this is a two dimensional list that will contain tuples controlling how the gapped traceback should proceed
        # entries will have the form (-1,-1) for ungapped tracebacks on the diagonal, (0,-1) for gaps "up"
        # (i.e. gap is in the target sequence) and (-1, 0) for gaps "left" (i.e. a gap is present in the query sequence)
        # (0,0) is also possible when there is no traceback path

        return

    def score(self):    # This method performs SW scoring and returns a string describing the resulting alignment

        highscore = high_i = high_j = 0       # highest scores encountered so far in the matrix

        best_q_alignment = []  # best alignment for the query sequence
        best_t_alignment = []  # best alignment for the target sequence

        ## It simplifies the logic to score the top and left side edges first,
        ## since these are always 0 or $matchReward and do not have a traceback.

        for i in range(self.querylen):  # score the left edge

            if self.query[i] == self.target[0]:   # is it a match?

                self.table[i][0] = self.match_reward

        for j in range(self.targetlen):  # score the top edge

            if self.target[j] == self.query[0]:   # is it a match?

                self.table[0][j] = self.match_reward

        # Now we can nested loop through the remainder of the matrix, scoring as we go

        for i in range(1, self.querylen):      # start these iterations from 1, not 0, as we have already done edges

            for j in range(1, self.targetlen):

                queryword = self.query[i: i + 1]  # An array slice is perhaps more natural in python than a substring
                targetword = self.target[j: j + 1]

                if queryword == targetword:     # Did we have a match at this position?

                    increment = self.match_reward

                else:

                    increment = self.mismatch_penalty

                matchscore = self.table[i - 1][j - 1] + increment   # increment will contain either a positive reward
                # or a negative penalty depending on whether we matched

                target_gap_score = self.table[i][j - 1] + self.gap_penalty    # scores associated with gapping
                query_gap_score = self.table[i - 1][j] + self.gap_penalty     # in either the target or query

                best_score = max(
                    (0, (0, 0)),                        # a 0 score will never have a traceback
                    (matchscore, (1, 1)),               # A match corresponds to a -1,-1 traceback
                    (target_gap_score, (0, 1)),         # A target gap corresponds to a 0, -1 traceback
                    (query_gap_score, (1, 0))           # A query gap corresponds to a -1, 0 traceback

                )
                # I am bundling my possible scores into a tuple with *another* tuple representing traceback offsets!

                self.table[i][j] = best_score[0]    # The first element in the tuple is the actual score to be recorded
                self.traceback[i][j] = best_score[1]    # The traceback offsets associated with the score are in a tuple

                # The above relies on the fact that the python max function will make its choice just based on the
                # first value in a tuple. The other information (offsets) just conveniently come along for the ride
                # this "trick" is in fact idiomatically pythonic and is commonly encountered.
                # Note how much cleaner this is than doing a series of if-elses to determine what to record for the
                # traceback value...

                if self.table[i][j] > highscore:    # This represents the "high road" approach.
                # "low road" would be >=

                    highscore = self.table[i][j]    # record the new high score
                    high_i = i                      # also record the i and j positions associated with that score
                    high_j = j

        ## Now we can go ahead and produce an output string corresponding to the best alignment

        i = high_i          # our approach is start at the high scoring box, and to trace our way back
        j = high_j

        while self.table[i][j] and i > -1 and j > -1:

            i_offset, j_offset = self.traceback[i][j]       # unpack the offset tuples stored in the traceback table

            if i_offset:
                best_q_alignment.append(self.query[i])
            else:
                best_q_alignment.append('-')                # if the value is a zero, we are gapping!

            if j_offset:

                best_t_alignment.append(self.target[j])

            else:
                best_t_alignment.append('-')                # if the value is a zero, we are gapping, now the other way

            i -= i_offset
            j -= j_offset

        best_q_alignment.reverse()  # flip 'em both once we are done, since we built them "end-to-beginning"
        best_t_alignment.reverse()

        # Alternatively, we could have built our alignments by adding things at the beginning using statements like
        # best_t_alignment.insert(0,'-') etc. But in Python inserting items at the beginning of a list is much slower
        # than appending at the end. We are better off appending at the end, then reversing the whole mess when done.

        return_string = '\nBest alignment had a score of ' + str(highscore) + ' and is:\n\nTarget:\t' + \
                        str(j + 2) + '\t' + ''.join(best_t_alignment) + '\n\t\t\t'

        for k in range(len(best_t_alignment)):     # t and q alignments should be the same length!

            if best_t_alignment[k] == best_q_alignment[k]:

                return_string += '|'    # Only put a bar if the two characters are identical at this position

            else:

                return_string += ' '    # otherwise just insert a space

        return_string += '\nQuery:\t' + str(i + 2) + '\t' + ''.join(best_q_alignment) + '\n'

        # The above statements just concatenate together a multi-line string that will correctly display
        # the best alignment when it is subsequently printed.

        return return_string

    def __str__(self):                          # This is a "special method attribute" that controls what the
    # string representation of objects of the SW class will look like.
    # It is invoked by print statements, which will print the return value

        lineout = 'Scoring table:\n\t' + '\t'.join(self.target) + '\n'
        # The above is just a fancy looking way to break the target string into tab-delimited individual characters

        for i in range(self.querylen):
            lineout += self.query[i] + "\t"
            for j in range(self.targetlen):

                lineout += str(self.table[i][j]) + "\t"

            lineout += '\n'

        lineout += '\n\nTraceback table:\n\t' + '\t'.join(self.target) + '\n'

        for i in range(self.querylen):

            lineout += self.query[i] + "\t"

            for j in range(self.targetlen):

                lineout += ''.join([str(k) for k in self.traceback[i][j]]) + "\t"  # prettying up the traceback tuples

            lineout += '\n'

        return lineout


# MAIN PROGRAM

A = smithwaterman('AAGCGCGTCGCGTTTGACATTCTAGAAGGCGC', 'ATGTTCGCGATTGACAATTCTTGAACGAGTCAG', 1, -2, -4)

print (A.score())

print (A)
