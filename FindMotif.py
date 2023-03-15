import math
import re
from math import log, log2, sqrt

class FindMotif:

    """
    FindMotif
    AUTHOR: Paul Fawcett
    VERSION: Python Version 2.0 (October 2019) closely inspired by an earlier PERL program by Jeff Elhai
                Version 2.0 is a conversion to Python 3.x
    PURPOSE: Calculates position-specific scoring matrix from aligned sequence
             Finds sequences in genome sequence that the PSSM fits best
    INPUT FILES:
             motif_input: input file of aligned motifs
                            alignments may have gaps, represented by "."
             sequence_file: input file containing sequences in FA format
                            May contain multiple contigs
    OUTPUT FILE:
                motif_hits: output file of best motifs within sequence file
                    Each line is name of sequence, score, start_coord, sequence
    """

    def __init__(self, motif_input=None, sequence_file=None, bg_freq=None):

        if not motif_input:     # If the user does not specify a motif file path, try and use this one as default
            self.motif_input = '71NpNtSm.txt'
        else:
            self.motif_input = motif_input

        if not sequence_file:     # If the user does not specify a sequence file path, try and use this one as default
            self.sequence_file = \
                'AnabaenaChromosome.nt'
        else:
            self.sequence_file = sequence_file

        if not bg_freq:
            self.bg_freq =  {               # Background frequency of nucleotides in genome
                                            'A' : 0.33264,  # If the user does not specify a background frequency,
                                            'T' : 0.33264,  # then default to these ones, taken from Anabaena intergenic sequences
                                            'G' : 0.16736,
                                            'C' : 0.16736
            }
        #Comprehension- python shortcut for building dictionaries
        self.log_p = {symbol: log(self.bg_freq[symbol]) for symbol in self.bg_freq.keys()}
        #print(self.log_p)

        self.log_q = {}
        self.q={}

        # The above is a dictionary comprehension.  These are conceptually very similar to the idea of a list
        # comprehension.  Note that this statement will only work in Python 2.7+, it is NOT supported in Python
        # versions 2.6x and below! If it doesn't work for you, you will need to unroll this into a standard 'for' loop

        # print (self.log_p) # uncomment this line to explore what the dict comprehension does!

        # Note that the log_p variable is a dict  frequency, at a particular position, of symbols in
        # that are read in from our aligned motif file as specified in motif_input
        # like log_p above, these values are log transformed for reasons that will subsequently become clear

        self.N = self.motif_len = 0
        self.motif_sequences = self.informational = []

        return

    def GetMotif(self):

        """
        Return array of aligned sequences
        Set the instance variable motif_length
        """

        motif_pattern = re.compile('\s([ACGTacgt]+)$')     # Note that our symbols here are only nucleotides
        motif_len = None
        self.motif_sequences = []

        print ('Opening motif file', self.motif_input)

        with open(self.motif_input, 'r') as f:

            for line in f:

                line = line.strip()        # Get rid of any extraneous whitespace

                match_object = motif_pattern.search(line)

                if match_object:
                    motif = match_object.group(1)

                    if not motif_len:
                        motif_len = len(motif)
                        print ('Motif of length', motif_len, 'has been detected.')

                    if motif_len == len(motif):
                        self.motif_sequences.append(motif.upper())  # add an uppercase version of the motif to the list

                    else:
                        print ('There seems to be a problem. Not all motifs appear to be the same length')
                        exit()

        self.N = len(self.motif_sequences)
        print ('Found', self.N, 'motifs of length', motif_len)

        self.motif_len = motif_len

        return self.motif_sequences   # returns a list of motif sequences in case this is used "standalone"

    def ReadContigs(self):

        """
        ReadContigs is a rather crude method for reading in FastA formatted files. Each annotation line is stored in a
        list called 'headers', while each associated sequence is stored in list called 'sequences'. Both are returned.
        There are several disadvantages to this method. First, it is not smart about multi-line annotations, and
        using a file with these would break it.  Also, it is crude in the sense that it reads an entire file into memory
        all at once, a practice that I have discouraged. For the Pythonistas in the audience, consider that a smarter
        version of this method might be set up as a class that has the behaviour of an iterator (i.e. implements
        next() and __iter__() methods).  In that scenario, a ReadContigs object would simply return a tuple of the next
        available header/sequence pair, until the file runs out, at which time it would throw a StopIteration exception.
        Non-Pythonistas: ignore the preceding two sentences!
        """

        headers = []
        sequences = []
        sequence = ''
        total_len = 0

        print ('Opening sequence file', self.sequence_file)

        with open(self.sequence_file, 'r') as f:

            for line in f:

                line = line.strip()

                if line.startswith('>'):

                    headers.append(line)

                    if sequence:
                        sequences.append(sequence)
                        total_len += len(sequence)
                        sequence = None

                else:

                    sequence += line.upper()

            sequences.append(sequence)          # Process the last sequence
            total_len += len(sequence)

        print ('Read in a total of', total_len, 'nucleotides of sequence data')

        return headers, sequences
        # keep in mind that each of these returned items is a list!

    def MakePSSM(self, B=0.1):

    # B is the total Pseudocount constant. Note the default value of 0.1

        counts = {}                     # A dict with keys that are tuples of nucleotide symbol and position in a motif
        # The values will correspond to the number of occurrences of that symbol at that
        # particular position in the motif alignment

        if not self.motif_sequences:    # If the we haven't got any motif sequences yet, go ahead and grab 'em

            self.GetMotif()

        if B == 'SQRTN':

            B = sqrt(self.N)



        # informational is a list containing positions in the
        # alignment that should be considered for scoring on the basis
        # that they are sufficiently information-rich.

        #   First, we calculate raw nucleotide counts at each position in the motif alignment

        for motif in self.motif_sequences:      # We do this by iterating over motifs, and then positions
            for position, symbol in (enumerate(motif)):

                #print (position, symbol)        # Uncomment this line to get a feel for what enumerate must be doing!

                if (symbol, position) in counts:        # It's a minor annoyance of Python that applying the
                    counts[(symbol, position)] += 1     # += operator to 'None' throws an error.
                    # Instead, we test for existence of the key in the dict
                    # before trying to increment the value.
                else:                                   # We could also use exception handling to deal with this
                    counts[(symbol, position)] = 1      # by detecting and handling a KeyError

        #print(counts)

        #    Now that we have raw counts for symbols at each position, we convert these to frequencies
        #    However, this is not as straightforward as you might guess!
        #
        #    The actual frequency is
        #
        #               counts[(symbol, position)] / self.N
        #
        #    where:
        #       counts(symbol, position) is the number of times the nucleotide appears at the position
        #       while self.N is the number of motifs found in the training session,
        #
        #    However, this approach will cause trouble for us later in cases where a symbol was not observed at
        #    a particular position.  So we modify our formula to include "pseudocounts" as described in class.
        #
        #     The adjusted formula is a little more complicated, see class PDF
        #
        #    but we will need the following -- counts and N are as before, and:
        #       q is the adjusted frequency of the nucleotide at the position (what data type?!)
        #       B is the weighting constant for pseudocounts
        #       bg_freq is the fraction that that nucleotide symbol appears in total DNA (again what data type?!)
        #
        #    Note that we will actually store log_q[(symbol, position)], the log of the adjusted
        #    frequency, rather than the frequency itself.  Why?

        for symbol in self.bg_freq.keys(): # This should be all the symbols that we know about -- hopefully just
        # 'A', 'C', 'G', and 'T' in the current example
            for position in range(self.motif_len):
                #update the cur_count so its never zero
                cur_count = counts.setdefault((symbol, position),0) # if no tupel give value of zero

                # The above makes sure that cur_count always has a value, even if the key specified by a the combination
                # of symbol and position didn't exist in the dict (i.e. we did not observe it in the training data).
                # The dict attribute .setdefault(key, default) returns either the value associated with the specified
                # key, or, if that does not exist, the specified default value.  This is handy, because we do not need
                # to first test for existence of the key to avoid the possibility of a KeyError.
                #print (self.N)
                #print (cur_count, position)

                #base frequency calculation

                # self.pssm = (cur_count + B*self.bg_freq[symbol])
                self.log_q[(symbol, position)] = log((cur_count + B*self.bg_freq[symbol])/ (self.N + B)) #add in psudo codes
                #print(self.log_q)

                #calc the frequency in real space instead of log space
                self.q[(symbol, position)] = (cur_count + B * self.bg_freq[symbol]) / (self.N + B)

                # Oh, oh!  This version of the PSSM calculation doesn't handle pseudocounts! Can you fix it?

                # Note that by default that the log function in Python uses e, the base of the natural logarithms

        self._GetInformational()

        return()


    def AnalyzeSequences(self, threshold, number_of_hits):

        """ AnalyzeSequences will score each PSSM-lengthed region in each of the sequences contained in the input file
        It will invoke the ReadContigs() method to populate the sequences list.
        Scores (along with other useful information such as contig name and position) that exceed some threshold value
        will be stored in the scores list. Note that the threshold value argument defaults to zero, so by default all
        scores will be recorded regardless of "goodness". We can pass a value of our choosing if we wish to limit the
        output to only "good" scores.  The 'scores' list is returned.
        """

        scores = []

        if not self.log_q:          # If we haven't yet constructed a PSSM manually, force it to happen

            print ('You have not yet created a PSSM using the MakePSSM method, so attempting to do so automatically..')
            self.MakePSSM()

        headers, sequences = self.ReadContigs()

        for seq_num, sequence in enumerate(sequences):
        # iterate over all sequences. Note that the same indices are used for sequences and headers, so I am here
        # also enumerating the seq_num for use later.

            for start_position in range(len(sequence) - self.motif_len):
            # Iterate over each possible starting position within a particular sequence

                if start_position % 20000 == 0: # Print updates of progress every so often. % is the modulus operator
                    print ('Sequence number:', seq_num, 'position', start_position)

                #slice notation
                region = sequence[start_position: start_position + self.motif_len]
                # Extract from the current sequence a region of size motif beginning at the current starting position

                region_score = 0

                for position in self.informational:      # iterate over each position and symbol in region
                    symbol = region[position]
                    # print(self.log_q[("A", position)])
                    # quit()
                    #get position and score based on the background score
                    #symbol = region[position]
                    #then add to region score
                    region_score += self.log_q[(symbol, position)] - self.log_p[symbol]

                    #these regions scores are not probs they are odds log ratios

                    # um, this doesn't look right! Can you fix this so that the score is calculated correctly?
                    # Don't forget that we also need to adjust for background frequencies!

                    # Note that 'position' here will also correspond to positions in our PSSM!

                if region_score > threshold:

                    scores.append((region_score, headers[seq_num], start_position, region))
                    # As you can see the 'scores' list has elements that are tuples bundling up all the things we might
                    # wish to know about a "hit".

                if (start_position > 50000) or (seq_num > 100):       # This is strictly for testing!
                    print ('TEST MODE ONLY!  Exceeded 50,000 nucleotides or 100 total sequences, so quiting analysis')
                    break

        #  Now that we have collected all the scores, we will sort them then chop out the extraneous hits beyond
        #  the number specified in the method argument 'number_of_hits'

        print ('Sorting hits, please wait.')

        scores.sort(reverse=True)
        # Python has a very powerful list method .sort() that does an "in-place" sort (the return value is None).
        # An alternative is the newlist = sorted(list) function, which will return a NEW list that is sorted while
        # preserving the old one. It is slightly less efficient.
        # Note also the use of the reverse parameter, which accepts a boolean value (default value is 'False')
        # Finally, keep in mind that we are sorting tuples.. Python does this by considering the first value in the
        # tuple first, which in our case is the score. It will use the next item, etc. if necessary to break ties.
        # if we wanted to sort on the second (or arbitrary elements) in the tuple we can say: .sort(key=itemgetter(i))
        # where i is the desired element. This requires 'from operator import itemgetter' at the beginning of the code.

        del scores[number_of_hits:]
        # Get rid of extraneous hits in the reverse sorted scores array.
        # I think this is a reasonably efficient way to to do this.  There are number of alternative ways this
        # might be done -- for instance 'scores[number_of_hits:] = []' should work, and there are other good ways too.
        # The obvious method using slice notation scores = scores[0:number_of_hits] is probably way less efficient!
        # Ambitious students could use the timeit module to figure out which alternative is fastest.

        return scores

    def _GetInformational(self):
        #input is counts or frequency of individual cols, they are
        #output is cols that we want to include

        #start our list that we will return with all valuable cols

        #we will use goodlist in analyze method
        self.goodlist = []


        #start with nested for loop that looks at all cols
        for position in range(self.motif_len): #iterate through cols
            #when we get to each col, calc prob for each symbol
            information = 2
            for symbol in self.bg_freq.keys(): #this is just using the keys as symbols ATCG

                #print(self.q[symbol, position])
                information += self.q[symbol, position] * log2(self.q[symbol, position])

                #make an if else statment that says that if number is greater than or less than threshold it is appended to valuable base pairs

                print(information)
            if information > 1.25: #replace the literal 1 with something defined earlier
                self.informational.append(position)

        print(self.informational)

        return   # We only want to return the appended list of things that we deemd most valuable goodlist
    #calculating based on input data and only adding to a list if it surpasses a specific number

def main():

    myMotif = FindMotif()
    # Instantiate a FindMotif object.  For now, we'll just use the default values for the various arguments.

    myMotif.MakePSSM(0.2)
    # Tell our FindMotif object to create a PSSM. We can pass as an argument a pseudocount constant
    # we can specify this either as a number, or as the string 'SQRTN', in which case B will be assigned the
    # square root of the number of rows in the training set alignment (i.e. the motif alignment). The default value
    # if nothing is specified is 0.1

    results = myMotif.AnalyzeSequences(-2**30, 10)

    # I am here initially specifying a very small threshold score of -2**30 so that all hits should be reported.
    # We can subsequently change this value to something larger to make our hit list more stringent.
    # I am also requesting a list of just the top 10 scoring hits. This also can & should be changed!
    myMotif._GetInformational()
    for i, hit in enumerate(results):

        print ('Hit', i + 1, '\t{}\t{}\t{}\t{}'.format(*hit)) # Again the "splat" operator breaks apart the hit tuples

if __name__ == '__main__':
    main()


