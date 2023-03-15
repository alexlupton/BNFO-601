# import find_matches
#
# mySeq = find_matches.Sequence(r'tester_file')
#
#
# print(mySeq.find_matches("acttgt", 0))
#
# print(mySeq.sequence)
# print(mySeq.revcomp(mySeq.sequence))


observed = 3
B=0.1
baseFreq= 0.2
totalCount= 5

#print((observed+B* baseFreq)/(5+B))
from math import log, sqrt
#log calculations
#score = -log(score/background frequency)
score= .004
backFreq= .2

print(-log(score/backFreq))

#for pssm start with alligned data, counting integers,