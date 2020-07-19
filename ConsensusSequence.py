"""
This python script makes a position weight matrix of the consensus sequence around TIS.
It reads a file of genomic sequences and gives another file of consensus sequence.
"""

#readFasta reads the fasta file given and
# formats it into the consensus sequence (length of 20 bp)
# The file has a ATG in the middle and has same number of nucleotides up/downstream of the ATG
def readFasta(file):

    #list to save the consensus sequences
    cons = []
    #read the fasta file line by line
    for seq in open(file, 'r').readlines():
        #find the length of the up/downstream sequences of TIS
        #   assumes that the upstream and downstream have the same length
        l = int((len(seq)-3)/2)
        #find the upstream consensus sequence (length 10)
        up = seq[l-10:l]
        #print(len(up))
        #find the downstream consensus sequence (length 10)
        down = seq[l+3:l+13]
        #print(len(down))
        #add the consensus sequence to the list
        cons += [up+down]

    return cons


#The consensus function takes in a file of TIS sequences and produces the position weight matrix of the consensus sequence
#   The sequence has equal number of up/downstream divided by ATG
def Consensus(inputfile, outputfile):
    import json

    #get the consensus sequences of each sequence in list form
    cons = readFasta(inputfile)
    l = len(cons)
    #Each list in consensus represent each row of pwm; each column represents the position on the sequence
    #consensus = [['A'],['C'],['G'],['T']]
    consensus = dict()
    position = [p for p in range(-10,13) if p not in [0,1,2]]

    #Create a loop ranging along the length of the consensus
    for i in range(0, 20):
        #Reset the values for each nucleotide to 0
        #   These will count the number of occurrences of the nucleotide on a certain position
        A,C,G,T = 0,0,0,0
        #Create another loop that reads each consensus sequence
        #   A nucleotide of certain position of all the sequences in the file will be put forward
        for j in cons:
            #Add to the nucleotide value if that nucleotide exists
            if j[i] == 'A':
                A += 1
            elif j[i] == 'C':
                C += 1
            elif j[i] == 'G':
                G += 1
            elif j[i] == 'T':
                T += 1

        #After adding all the occurences of each nucleotide, divide it by the total number of seqneuces
        A = round(A/l*100)
        C = round(C/l*100)
        G = round(G/l*100)
        T = round(T / l * 100)
        #add all to dictionary
        consensus[position[i]] = [A,C,G,T]

    with open(outputfile, 'w') as out:
        json.dump(consensus, out, indent=10)

Consensus('at_pos_dic2013.txt','testfile.txt')
