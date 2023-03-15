import re

class Sequence(object):
    """Version 2.0 February 2020 -- conversion to Python 3
    A single sequence as read from a FASTA-formatted input file.

    Stores header data and sequence data.
    """

    def __init__(self, input_file):
        self.header = ''
        self.sequence = ''

        self.read_sequence(input_file)

    def __len__(self):
        return len(self.sequence)

    def read_sequence(self, input_file):
        """Read header and sequence data from the specified input file.

        Parameters:
            input_file (str): the path to the FASTA format sequence file.
        """
        with open(input_file, 'r') as ifh:
            for line in ifh:
                line = line.strip()
                if line.startswith('>'):
                    self.header = line.lstrip('>')
                else:
                    self.sequence += line

    def find_matches(self, pattern, mismatches_allowed):
        """Find the portions of the sequence that match the pattern.

        Searches both the forward and reverse strand for matches.  Allows for a user-specified number of mismatches.
        Currently, mismatches > 1 are not supported.

        Parameters:
            pattern (str): the regex pattern to search sequence with
            mismatches_allowed (int): the max number of permitted mismatches

        Returns:
            a list of matches
        """
        results = []
        matches = None
        if mismatches_allowed == 0:
            matches = re.finditer(pattern, self.sequence, re.IGNORECASE)
            if matches:
                for m in matches:
                    start = m.start()
                    end = m.end()
                    results.append('Match: {}...{} D'.format(start, end))

            matches = re.finditer(self.revcomp(pattern), self.sequence, re.IGNORECASE)
            if matches:
                for m in matches:
                    start = m.start()
                    end = m.end()
                    results.append('Match: {}...{} R'.format(start, end))

        elif mismatches_allowed == 1:
            print ('mismatches allowed!')
            mismatch_pattern = []

            for i in range(len(pattern)):
                new_pattern = pattern[:i]
                new_pattern = '{}.{}'.format(new_pattern, pattern[i + 1:])
                mismatch_pattern.append(new_pattern)
            mismatch_pattern = '|'.join(mismatch_pattern)
            mismatch_pattern = '({})'.format(mismatch_pattern)
            print ('Mismatch pattern is', mismatch_pattern)
            matches = re.finditer(mismatch_pattern, self.sequence, re.IGNORECASE)
            if matches:
                for m in matches:
                    start = m.start()
                    end = m.end()
                    results.append('Match: {}...{} D'.format(start, end))
            print ('Reverse pattern is', self.revcomp(mismatch_pattern))
            matches = re.finditer(self.revcomp(mismatch_pattern), self.sequence, re.IGNORECASE)
            if matches:
                for m in matches:
                    start = m.start()
                    end = m.end()
                    results.append('Match: {}...{} R'.format(start, end))
        else:
            print('Error: not implemented for mismatches > 1')

        return results

    def revcomp(self, pattern):
        """Reverse complement a given sequence.

        If the sequence comes in the form of a regex with parentheses, they will be stripped and replaced provided
        they are located only at the beginning and end of the pattern.

        Parameters:
            pattern (str): the sequence or regex pattern to reverse complement

        Returns:
            a string containing the reverse complement of pattern
        """
        if '(' in pattern or ')' in pattern:
            pattern = pattern.replace('(', '')
            pattern = pattern.replace(')', '')
        trans_table = str.maketrans('ACTGactg', 'TGACtgac')
        #print ('({})'.format(pattern.translate(trans_table)[::-1]))
        return '({})'.format(pattern.translate(trans_table)[::-1])


def main():
    s = Sequence(r'anabseq.nt')
    print(s.header)
    print(len(s))
    results = s.find_matches('AATA', 1)
    print('\n'.join(results))


if __name__ == '__main__':
    main()
