#random sequence genration
#alex Lupton
#021723
#BNFO 601

import random
#define weighted list for characters to be randomly grabbed
myChars = ['A'] * 100 + ['C'] * 400 + ['G'] * 400 + ['T'] * 100
#define length of sequence
seqLen= 20
'''
for i in range(seqLen):
    #make an empty list
    strand=""
    #generates a random ineger every loop
    myint = random.randint(0, 999)
    #prints a random character atcg based on the random integer and adds it to list
    strand += (myChars[myint])
    print(myChars[myint])
print(strand)
'''
#or use numPy to generate list


    def randSequence(self, seqLength) :
        from numpy.random import choice

        sampleList = ['a', 't', 'c', 'g']
        randomNumberList = choice(
            sampleList, seqLen, p=[0.1, 0.1, 0.4, 0.4])
        return(randomNumberList)

print(randSequence(100))

from find_matches import Sequence
object = Sequence.revcomp()

#generate the reverse complimentary
reverse= randomNumberList.revcomp
print(reverse)


#list that we are searching for
searchList = ['c','c','g']
numMisMatch = 1
match= randomNumberList.find_matches(searchList, numMisMatch )
#print(find_matches(randomNumberList,searchList , 1))
print(match)
