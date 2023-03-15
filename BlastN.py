_author__ = 'Wombat'


class BLAST:

    """ BlastN.py takes a query string and a target string and uses the BLAST (Basic Local Alignment Search Tool)
    algorithm to find regions of local alignment between the two. This toy version of BLAST employs the word heuristic,
    but has greatly simplified logic for identifying and extending seed matches.
    Version 5, by Paul Fawcett, 31 November 2019
    Changes from Version 1:
        - Fixes a bug identified by John Drake in the _extend_alignment method relating to pre-filling the first column
        - Fixes a bug identified by Kristen Kelly relating to the choice of query words;
        - Fixes a bug identified by William Potts where the word_score was being inappropriately multiplied
            by match_reward (word_score is already the product of word_size and match_reward);
        - Word size has been adjusted to 11, the default for NCBI BlastN
        - Conversion from Python 2.x to 3.x
    """

    FORWARD = 1         # These are class variables shared by all instances of the BLAST class
    BACKWARD = -1
    ROW = (0, 1)
    COLUMN = (1, 0)

    def __init__(self, query=None, target=None, word_size=13, match_reward=1, mismatch_penalty=-2,
                 gap_open=-5, gap_extend=-2):

        self.query = query          # This is the string corresponding to the query sequence
        self.target = target        # This is the string corresponding to the target sequence
        self.word_size = word_size  # Size of the seed word for initiating extensions
        self.match_reward = match_reward
        self.word_score = self.word_size * self.match_reward    # something different required for PBLAST!

        self.mismatch_penalty = mismatch_penalty    # These are our key scoring parameters
        self.gap_open = gap_open
        self.gap_extend = gap_extend

        self.querylen = len(query)
        self.targetlen = len(target)

        self.blast_table = {}       # Our main dynamic programming table containing scores
        self.traceback_table = {}   # A corresponding table for recording the tracebacks
        self.target_index = {}

        return

    def score(self):    # This method performs BLAST scoring and returns a string describing the resulting alignment

        result_summary = []         # A list, for now, that will store results of the alignments

        if not self.target_index:   # if this is the first time scoring we should index the target

            for i in range(len(self.target) - self.word_size + 1):

                word = self.target[i: i + self.word_size]

                if word in self.target_index:

                    self.target_index[word].append(i)       # A dict of lists is an efficient structure for this index.
                    # The list items are word coordinates in the target.
                else:

                    self.target_index[word] = [i]

        #print self.target_index
        ## First we must iterate through words in the query:

        query_position = 0
        while query_position < self.querylen - self.word_size + 1:

            #print "Query position is", query_position
            query_word = self.query[query_position:query_position + self.word_size]

            # for each word we then ask if it is present in the target index
            if query_word in self.target_index:

                for target_position in self.target_index[query_word]:

                    print ('Searching for seed', query_word, 'at target position', target_position)

                    #print "Extending forward"
                    forward_score, forward_extension_q, forward_extension_t = \
                        self._extend_alignment(query_position, target_position, self.FORWARD)

                    #print "Extending backwards"
                    backward_score, backward_extension_q, backward_extension_t = \
                        self._extend_alignment(query_position, target_position, self.BACKWARD)

                    q_result = backward_extension_q[:-1] + query_word + forward_extension_q[1:]
                    t_result = backward_extension_t[:-1] + query_word + forward_extension_t[1:]

                    # Note that the last character of a backward extension, and the zeroth character of a forward
                    # extension overlap with the query word and should therefore be discarded - thus the slice notation.

                    score = forward_score + backward_score - self.word_score
                    # We need to make sure that we don't double count the seed score!

                    query_begin = query_position - len(backward_extension_q) + 2

                    target_begin = target_position - len(backward_extension_t) + 2

                    result_summary.append((score, q_result, t_result, query_begin, target_begin))

                    alignment_string = '\nAlignment had a score of ' + str(score) + ' and is:\n\nTarget:\t' + \
                                       str(target_begin) + '\t' + str(t_result) + '\n\t\t\t'

                    for k in range(len(t_result)):     # t and q alignments should be the same length!

                        if t_result[k] == q_result[k]:

                            alignment_string += '|'
                            # Only put a bar if the two characters are identical at this position
                        else:
                            alignment_string += ' '    # otherwise just insert a space

                    alignment_string += '\nQuery:\t' + str(query_begin) + '\t' + str(q_result) + '\n'

                    print (alignment_string)
                    # The above statements just concatenate together a multi-line string that will correctly display
                    # the best alignment when it is subsequently printed.

            query_position += 1

        return result_summary

    def _extend_alignment(self, query_start, target_start, direction):

        """ This private method attempts to extend an alignment in the forward and backward direction
        depending on the value of the direction flag, which here takes the value 1 (for forward extension) or
        -1 for backward.For clarity these constants are defined by the class variables self.FORWARD and self.BACKWARD
        """

        self.high_score = self.word_score

        # highest scores encountered so far will always initially be the word_score * match_reward

        self.high_q_pos = self.high_t_pos = 0

        if direction == self.FORWARD:            # We start with the 0,0 position representing the last character
            query_start += self.word_size - 1    # of the seed word for forward extensions.
            target_start += self.word_size - 1   # For backward extensions, leave it as it is (i.e. zeroth character)

        self.blast_table = dict()
        # The BLAST table is a dict of tuples. Each tuple represents a (query, target) position
        # this sparse representation will be much more efficient than using a 2D list

        self.blast_table[0, 0] = self.high_score   # initialize the top left corner with the word score
        self.high_q_pos = 0
        self.high_t_pos = 0

        self.traceback_table[0, 0] = (1, 1)
        # There is no traceback path for the origin, but the program logic elsewhere dictates that we provide one

        cur_t_pos = 1   # we are going to score the edges first (top and left), which can *only* ever be gaps back
        # to the origin.  i.e. the question of matching or not matching is completely irrelevant here.
        # We start by scoring the top edge, beginning with position 1..

        cur_score = max(0, self.blast_table[(0, 0)] + self.gap_open)    # first one always a gap open

        while cur_score:        # only keep going as long as we have non-zero values

            self.blast_table[(0, cur_t_pos)] = cur_score    # only record non-zero values
            self.traceback_table[(0, cur_t_pos)] = (0, 1)   # record a target gap in the traceback table

            cur_score = max(0, self.blast_table[(0, cur_t_pos)] + self.gap_extend)  # any subsequent are extends
            cur_t_pos += 1

        cur_t_pos = 0   # Now we do the same thing for the left edge as we just did for the top edge
        cur_q_pos = 1

        cur_score = max(0, self.blast_table[(0, 0)] + self.gap_open)    # first one always a gap open

        while cur_score:        # only keep going as long as we have non-zero values

            self.blast_table[(cur_q_pos, 0)] = cur_score    # only record non-zero values
            self.traceback_table[(cur_q_pos, 0)] = (1, 0)   # record a query gap in the traceback table
            cur_score = max(0, self.blast_table[(cur_q_pos, 0)] + self.gap_extend)
            cur_q_pos += 1

        #print "blast table 0,0 is", self.blast_table[0, 0], "and high score is", self.high_score

        # alright, finished with edges. Note that high scores can NEVER occur in an edge so these were not considered.
        # Henceforth, however, we will need to think about this.

        cur_t_pos = 0   # Start at the first position
        cur_q_pos = 0

        # Now we will score the table, proceeding according to the algorithm description: first incrementing along
        # the diagonal, then scoring the adjacent row, then the column below
        # Unlike Smith Waterman, the matrix is no longer of defined size, so we need to use while loops instead of for

        while True:  # I think it's cleaner to affirmatively break out of this main loop. Too bad Python has no do-while

            cur_t_pos += 1  # Advance along the diagonal by incrementing
            cur_q_pos += 1  # Remember, these refer to coordinates in our table, not in the actual target or query

            # Probably we need to do some bounds checking here too with respect to absolute position in the query and
            # target similar to what is done in the _fill_in_row_or_column method

            #print "Beginning row starting at", cur_q_pos, cur_t_pos, "of the blast table"
            max_in_row = self._fill_in_row_or_column(cur_q_pos, cur_t_pos, query_start, target_start,
                                                     direction, self.ROW)

            #print "Max in row was ", max_in_row
            #print "Beginning column starting at", cur_q_pos, cur_t_pos, "of the blast table"
            max_in_column = self._fill_in_row_or_column(cur_q_pos, cur_t_pos, query_start,
                                                        target_start, direction, self.COLUMN)

            #print "Max in column was ", max_in_column

            if not max(max_in_row, max_in_column):
                break   # If the maximum value we encounter in both the rows and columns is zero, we are done building

        #print "Finished building a matrix"

        best_q_alignment = []  # best partial alignment for the query sequence
        best_t_alignment = []  # best partial alignment for the target sequence

        ## Now we can go ahead and produce an output string corresponding to the best alignment

        cur_q_pos = self.high_q_pos          # our approach is start at the high scoring box, and to trace our way back
        cur_t_pos = self.high_t_pos

        while cur_q_pos >= 0 and cur_t_pos >= 0 and self.blast_table.setdefault((cur_q_pos, cur_t_pos), 0):

            q_offset, t_offset = self.traceback_table[cur_q_pos, cur_t_pos]
            # unpack the offset tuples stored in the traceback table

            if q_offset:

                try:
                    best_q_alignment.append(self.query[query_start + cur_q_pos * direction])

                except IndexError:

                    print ('YO!', query_start, cur_q_pos, direction, query_start + cur_q_pos * direction)
                    print ('Best_q_alignment', best_q_alignment)
                    quit()   # There has been some kind of problem

            else:
                best_q_alignment.append('-')                # if the value is a zero, we are gapping!

            if t_offset:

                best_t_alignment.append(self.target[target_start + cur_t_pos * direction])

            else:
                best_t_alignment.append('-')                # if the value is a zero, we are gapping, now the other way

            cur_q_pos -= q_offset       # Note that we are subtracting positively valued offsets.
            cur_t_pos -= t_offset       # This design choice makes later printing a traceback table a lot prettier.

        # Alternatively, we could have built our alignments by adding things at the beginning using statements like
        # best_t_alignment.insert(0,'-') etc. But in Python inserting items at the beginning of a list is much slower
        # than appending at the end. We are better off appending at the end, then reversing the whole mess when done.

        #print "Returning information about a partial alignment", self.high_score, best_q_alignment, best_t_alignment

        # flip 'em both once we are done, since we built them "end-to-beginning". Note that we don't need to flip
        # sequences corresponding to backwards extensions!

        if direction == self.FORWARD:
            best_q_alignment.reverse()
            best_t_alignment.reverse()

        return self.high_score, ''.join(best_q_alignment), ''.join(best_t_alignment)

    def _fill_in_row_or_column(self, cur_q_pos, cur_t_pos, query_start, target_start, direction, row_or_column):

        """This private method will fill in a row or column, depending on the tuple passed in the row_or_column argument
        Each row or column is filled in until a zero-valued result is obtained.
        """

        #print "filling in a row or column"

        max_in_current_row_or_column = 0
        q_add, t_add = row_or_column

        # These variables will control whether we fill in a row or a column.  If the argument row_or_column = (0,1)
        # we will end filling in a row.  If the argument is assigned (1,0) we will fill a column

        while True:

            query_position = query_start + cur_q_pos * direction    # remember, direction here is either -1 or 1
            target_position = target_start + cur_t_pos * direction  # so is a positive or negative offset multiplier

            # query and target position variables here refer to the actual (absolute) position within the query
            # and target sequences respectively

            if (query_position < 0) or (target_position < 0):

                #print "Ran out of query or target sequence while attempting backwards extension"
                break       # we can go no further

            if (query_position >= self.querylen) or (target_position >= self.targetlen):
                #print "Ran out of q or t while attempting forwards extension", query_position, target_position

                break       # again, we can go no further

            q_char = self.query[query_position]
            t_char = self.target[target_position]

            #print "comparing", q_char, query_position, "to", t_char, target_position

            if q_char == t_char:

                increment = self.match_reward

            else:

                increment = self.mismatch_penalty

            match_score = self.blast_table[(cur_q_pos - 1, cur_t_pos - 1)] + increment

            # improvement for later - decide whether to apply gap opening or gap extension penalties
            # for the moment just set gap increment to the gap_open value

            increment = self.gap_open

            # scores associated with gapping in either the target or query
            target_gap_score = self.blast_table.setdefault((cur_q_pos, cur_t_pos - 1), 0) + increment
            query_gap_score = self.blast_table.setdefault((cur_q_pos - 1, cur_t_pos), 0) + increment

            best_score = max(
                (0, (0, 0)),                            # a 0 score will never have a traceback
                (match_score, (1, 1)),                  # A match corresponds to a -1,-1 traceback
                (target_gap_score, (0, 1)),             # A target gap corresponds to a 0, -1 traceback
                (query_gap_score, (1, 0))               # A query gap corresponds to a -1, 0 traceback
            )

            if not best_score[0]:

                break

            self.blast_table[cur_q_pos, cur_t_pos] = best_score[0]
            # The first element in the tuple is the actual score to be recorded
            #print "Recording", best_score[0], "at position", cur_q_pos, cur_t_pos
            self.traceback_table[cur_q_pos, cur_t_pos] = best_score[1]
            # The traceback offsets associated with the score are in a tuple as described earlier

            if best_score[0] >= self.high_score:
            # This represents the "high road" approach. "low road" would simply be >

                self.high_score = best_score[0]  # record the new high score
                self.high_q_pos = cur_q_pos      # also record the i and j positions associated with that score
                self.high_t_pos = cur_t_pos

            if best_score[0] > max_in_current_row_or_column:
                max_in_current_row_or_column = best_score[0]

            # The maximum in a particular row or column is different from the overall high score! We actually
            # only care if this value is non-zero, as this will tell us that another iteration along the diagonal is
            # required.

            cur_t_pos += t_add      # We end up adding either a zero or a one to these depending on
            cur_q_pos += q_add      # whether we are filling in a row or a column, setting us up for the next iteration

        return max_in_current_row_or_column

    def __str__(self):

        """ This is a "special method attribute" overwriting the __str__ method defined in object.
        __str__ controls what the string representation of objects of the BLAST class will look like.
        It is invoked by print statements, which will print the return value. The bad news is that the routine here
        was more-or-less just lifted from the old Smith Waterman program. However, BLAST uses a fundamentally
        different sort of data structure for representing the blast and traceback tables.
        Can you fix this method so that it does something useful?
        """

        lineout = 'Scoring table:\n\t' + '\t'.join(self.target) + '\n'
        # The above is just a fancy looking way to break the target string into tab-delimited individual characters

        for i in range(self.querylen):
            lineout += self.query[i] + "\t"
            for j in range(self.targetlen):

                lineout += str(self.blast_table[i, j]) + "\t"

            lineout += '\n'

        lineout += '\n\nTraceback table:\n\t' + '\t'.join(self.target) + '\n'

        for i in range(self.querylen):

            lineout += self.query[i] + "\t"

            for j in range(self.targetlen):

                lineout += ''.join([str(k) for k in self.traceback_table[i, j]]) + "\t"
                # just prettying up the traceback tuples

            lineout += '\n'

        return lineout

# MAIN PROGRAM

A = BLAST("AAGCGCGTCGCGTTTGACATTCTAGAAGGCGC", "ATGTTCGCGATTGACAATTCTTGAACGAGTCAG")

#A = BLAST('AACCGGCCCTTCCGGCCG', 'AACCGGCCCAACCGGCCG')

#A = BLAST('AAAA', 'AAAA')

print (A.score())
