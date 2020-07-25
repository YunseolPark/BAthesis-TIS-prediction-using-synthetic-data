"""
July 12, 2020
"""

"""
This python script generates a synthetic TIS dataset for training a prediction model.
It makes use of consensus sequence, upstream start codon, downstream stop codon, splice site, and nucleotide frequency to generate the dataset.
It takes a file containing the position weight matrix (pwm).
It creates two files as outputs, each containing the positive and negative samples.

It is a more simplified, randomized version of GenerateTIS.py and also has codon frequency added as a feature.
"""

import random
import json

'''
This function forms the basic structure of the synthetic sequence into a list of nucleotides.
It takes in the start codon and a certain length where len(seq) = length*2+3.
'''
def BasicStructure(start, length):
    seq = []
    seq += 'u' * length     #upstream sequence
    seq += start            #start codon
    seq += 'd' * length     #downstream sequence

    #print(len(seq),"".join(i for i in seq))
    return seq


'''
This function takes in either four numbers (probabilities) or a list of the four numbers
and gives a certain nucleotide that could be selected based on the probability.
'''
def PWMtoBase(A,C=None,G=None,T=None):
    #If only a list is given, then separate it into the corresponding nucleotide.
    if C == None and type(A) is list:
        C,G,T = A[1],A[2],A[3]
        A = A[0]

    total = A+C+G+T
    prob = random.uniform(1,total)
    if prob <= A:
        return 'A'
    elif prob > A and prob <= A+C:
        return 'C'
    elif prob > A+C and prob <= A+C+G:
        return 'G'
    else:
        return 'T'


'''
This function adds the consensus sequence around the start codon for the positive sample.
It takes a certain length, a sequence, a file with consensus, and l which is the length of consensus.
'''
def ConsensusSequence(length, seq, conFile):
    with open(conFile, 'r') as fp:
        consensus = json.load(fp)       #consensus=dictionary with positions as keys and list of probabilities as value

    #keys span from negative to positive number, with ATG being 0,1,2
    for i in consensus.keys():
        seq[length+int(i)] = PWMtoBase(consensus[i][0],consensus[i][1],consensus[i][2],consensus[i][3])

    #print(len(seq),"".join(i for i in seq))
    return seq


'''
This function adds the upstream start codon to the negative samples.
It takes in the start codon, length, sequence, and the length of consensus sequence.
'''
def UpstreamStart(start, length, seq, l):
    # there can be up to two upstream start codons
    n = random.randint(0, 2)

    # set in_len so that consensus sequence and any codon overlaping it will be disregarded for insert site.
    in_len = int((length - (l + 5)) / 3)  # divide by 3 to select each codon not nucleotide

    # select random codon position to insert start codon
    if n == 1:
        pos1 = random.randint(1, in_len)  # the first codon place is not considered to avoid cutting it out later on
        seq[pos1 * 3:pos1 * 3 + 3] = start
    elif n == 2:
        pos1 = random.randint(1, in_len - 1)
        pos2 = random.randint(pos1 + 1, in_len)
        seq[pos1 * 3:pos1 * 3 + 3] = start
        seq[pos2 * 3:pos2 * 3 + 3] = start

    #print(len(seq),"".join(i for i in seq))
    return seq


'''
This function adds a downstream or upstream stop codon for negative or positive samples.
It takes a list of stop codons, length, sequence, and length of consensus, and either pos or neg.
'''
def Stop(stop_list, length, seq, l):
    stop = random.choice(stop_list)     #select a random stop codon

    in_len = int(length / 3)
    stop_site = random.randint(in_len + int((l+5)/3), in_len*2-2)  # last codon place not considered to avoid cutting it out
    seq[stop_site * 3:stop_site * 3 + 3] = stop

    #print(len(seq),"".join(i for i in seq))
    return seq


'''
This function adds splice site consensus sequence to upstream for negative and downstream for positive sample.
It takes in length, sequence, length of consensus, 'pos' or 'neg' to indicate the sample.
'''
def SpliceSite(length, seq, l, train, cut1, cut2):
    #Donor splice site
    splice = [[352,361,154,132],[621,131,87,159],[85,45,777,91]]
    #add splice site downstream for positive and upstream for negative sample
    if train == 'pos':
        n = random.randint(length+l+3, length*2+2-cut2-len(splice))
    else:
        n = random.randint(cut1, length-11-len(splice))

    for i,s in enumerate(splice):
        seq[n+i] = PWMtoBase(s)

    #print(len(seq),"".join(i for i in seq))
    return seq


"""
This function adds codons to the rest of the sequence.
It takes in length of sequence, the sequence, and two files that give the codon frequency of upstream and downstream.
"""
def CodonUsage(length, seq, CodonFile):
    with open(CodonFile[0],'r') as fp:
        up = json.load(fp)      #dictionary of upstream codon frequency
    with open(CodonFile[1],'r') as fp:
        down = json.load(fp)    #dictionary of downstream codon frequency
    dna = ['A', 'C', 'G', 'T']
    up_key,up_value,down_key,down_value = [],[],[],[]      #list for storing codon and its frequency for up- and downstream
    for (k1,v1),(k2,v2) in zip(up.items(),down.items()):
        up_key.append(k1)
        up_value.append(v1)
        down_key.append(k2)
        down_value.append(v2)

    #add codon for upstream sequence
    for u in range(0, length, 3):
        #only add codons in position where no nucleotide is added to
        check = any(i in seq[u:u+3] for i in dna)
        if check == True:
            continue
        codon = random.choices(up_key, up_value)[0]
        seq[u:u+3] = codon

    # add codon for downstream sequence
    for d in range(length+3, len(seq), 3):
        # only add codons in position where no nucleotide is added to
        check = any(i in seq[d:d + 3] for i in dna)
        if check == True:
            continue
        codon = random.choices(down_key, down_value)[0]
        seq[d:d+3] = codon

    #print(len(seq),"".join(i for i in seq))
    return seq


'''
This function is used to fill all the bases in the sequence other than those already filled.
It takes in length, sequence, and either pos or neg to indicate the sample.
'''
def NucleotideFrequency(length, seq, train):
    #negative training set
    if train == 'neg':
        for i in range(0,len(seq)):
            if seq[i] not in ['u','d']:     #do not add bases in positions already filled
                continue
            #fill with probability for both upstream and downstream of start codon
            if i in range(0, length):
                seq[i] = PWMtoBase(31.4,18.4,18.5,31.7)
            elif i in range(length+3, len(seq)):
                seq[i] = PWMtoBase(31.3,18.1,18.5,31.4)
        #print("".join(i for i in seq))
        return seq

    #for positive training set
    for j in range(0, len(seq)):
        if seq[j] not in ['u','d']:     #do not add bases in positions already filled
            continue
        # fill with probability for both upstream and downstream of start codon
        if j in range(0, length):
            seq[j] = PWMtoBase(31.1,19.6,15.5,33.8)
        elif j in range(length + 3, len(seq)):
            seq[j] = PWMtoBase(26.1,23.0,20.4,29.7)
    #print(len(seq),"".join(i for i in seq))
    return seq


'''
This function generates TIS sequence by combining the above functions.
It takes in either 'pos' or 'neg, length, file containing consensus, length of consensus sequence,
codon usage files, donor, start and stop arguemtns.
'''
def Generate(train, fL, conFile, l, CodonFile, donor, startcheck=None, stopcheck=None):
    #convert final length (fL) to two lengths that are equal and are divisible by 3
    fL -= 3     #remove the length for start codon
    cut = 0     #value needed to cut the sequence to fit final length
    #add nucleotide to make length in upstream and downstream equal
    if fL%2 == 0:
        length = int(fL/2)
    elif fL%2 == 1:
        length = (fL//2)+1
        cut += 1
    #add nucleotide(s) to make length multiple of 3 (codon)
    if length%3 == 1:
        length += 2
        cut += 2*2
    elif length%3 == 2:
        length += 1
        cut += 1*2
    #cut1 is for upstream and cut2 is for downstream
    cut1 = cut//2
    cut2 = cut - cut1

    start = ['A','T','G']
    stop_list = [['T', 'A', 'A'], ['T', 'A', 'G'], ['T', 'G', 'A']]

    #generate the sequence
    seq = BasicStructure(start, length)
    if conFile!=None and train=='pos':
        seq = ConsensusSequence(length, seq, conFile)       #add consensus sequence
    if startcheck==True and train=='pos':
        seq = UpstreamStart(start, length, seq, l)       #add upstream start codon
    if stopcheck==True and train=='neg':
        seq = Stop(stop_list, length, seq, l)     #add downstream stop
    if donor==True:
        seq = SpliceSite(length, seq, l, train, cut1, cut2)     #add splice site
    if CodonFile!=None:
        seq = CodonUsage(length, seq, CodonFile)    #fill with codons
    seq = NucleotideFrequency(length, seq, train)       #fill the rest with nucleotides

    #Remove any downstream stop codon in the positive set
    if train == 'pos':
        for p in range(0, len(seq), 3):
            if p >= length+3 and seq[p:p + 3] in stop_list:
                seq = seq[:p] + [PWMtoBase(26.1, 23.0, 20.4, 29.7)] \
                            + [PWMtoBase(26.1, 23.0, 20.4, 29.7)] \
                            + [PWMtoBase(26.1, 23.0, 20.4, 29.7)] + seq[p+3:]

    #cut the sequence to match the final length (fL)
    seq = seq[cut1:len(seq)-cut2]
    #print(len(seq),"".join(i for i in seq))
    return seq


'''
This function writes the generated TIS dataset into a file.
It takes in the number of sequences, output files for positive and negative set, final length, length of consensus,
and six arguments for features: a file containing consensus, boolean for start and stop codon, donor splice sites,
and a list of files containing codon usage.
It outputs the two file containing synthetic datasets.
'''
def WriteTIS(rows, posFile, negFile,fL=300,l=10,conFile=None,startcheck=True,stopcheck=True,donor=True,CodonFile=None):
    tis_pos = open(posFile, "w+")
    tis_neg = open(negFile, "w+")

    print("consensus: {}, start: {}, stop: {}, donor: {}, codon: {}".format(
        isinstance(conFile,str),startcheck,stopcheck,donor,isinstance(CodonFile,list)))

    if CodonFile!=None:
        posCodon = CodonFile[0:2]      #CodonFile for positive set
        negCodon = CodonFile[2:4]      #CodonFile for negative set
    else:
        posCodon, negCodon = None, None

    #Generate TIS dataset
    for i in range(0, rows):
        positive = "".join(i for i in Generate('pos',fL,conFile,l,posCodon,donor,startcheck=startcheck,stopcheck=stopcheck))
        negative = "".join(j for j in Generate('neg',fL,conFile,l,negCodon,donor,startcheck=startcheck,stopcheck=stopcheck))

        tis_pos.write(positive+'\n')
        tis_neg.write(negative+'\n')
    tis_pos.close()
    tis_neg.close()



print("########## {:^50s} ##########".format("Start generating synthetic datasets"))
#Dataset with all features
WriteTIS(rows=27102,posFile='arabTIS.pos',negFile='arabTIS.neg',conFile='consensus_sequence.txt')
print("########## {:^50} ##########".format("Generated TIS"))
print()
#Dataset without consensus sequence
print("########## {:^50s} ##########".format("Start generating synthetic dataset without consensus sequence"))
WriteTIS(rows=27102,posFile='arabTIScon.pos',negFile='arabTIScon.neg')
print("########## {:^50} ##########".format("Generated TIS without consensus sequence"))
print()
#Dataset without upstream start codon
print("########## {:^50s} ##########".format("Start generating synthetic dataset without upstream start codon"))
WriteTIS(rows=27102,posFile='arabTISstart.pos',negFile='arabTISstart.neg',conFile='consensus_sequence.txt',startcheck=False)
print("########## {:^50} ##########".format("Generated TIS without upstream start codon"))
print()
#Dataset wtihout stop codon
print("########## {:^50s} ##########".format("Start generating synthetic dataset without stop codon"))
WriteTIS(rows=27102,posFile='arabTISstop.pos',negFile='arabTISstop.neg',conFile='consensus_sequence.txt',stopcheck=False)
print("########## {:^50} ##########".format("Generated TIS without stop codon"))
print()
#Dataset without donor splice sites
print("########## {:^50s} ##########".format("Start generating synthetic dataset without donor splice site"))
WriteTIS(rows=27102,posFile='arabTISdonor.pos',negFile='arabTISdonor.neg',conFile='consensus_sequence.txt',donor=False)
print("########## {:^50} ##########".format("Generated TIS without donor splice site"))
print()
#Dataset with codon usage
print("########## {:^50s} ##########".format("Start generating synthetic dataset with codon usage"))
c = ['codonUp_pos.txt','codonDown_pos.txt','codonUp_neg.txt','codonDown_neg.txt']
WriteTIS(rows=27102,posFile='arabTIScodon.pos',negFile='arabTIScodon.neg',conFile='consensus_sequence.txt',CodonFile=c)
print("########## {:^50} ##########".format("Generated TIS with codon usage"))

print()
print()

print("########## {:^50s} ##########".format("Start generating synthetic dataset with only consensus"))
#Dataset with only consensus
WriteTIS(rows=27102,posFile='arabTIS_wcon.pos',negFile='arabTIS_wcon.neg',conFile='consensus_sequence.txt',startcheck=False,stopcheck=False,donor=False)
print("########## {:^50} ##########".format("Generated TIS with consensus"))
print()
print("########## {:^50s} ##########".format("Start generating synthetic dataset with only upstream start"))
#Dataset with only upstream start codon
WriteTIS(rows=27102,posFile='arabTIS_wstart.pos',negFile='arabTIS_wstart.neg',stopcheck=False,donor=False)
print("########## {:^50} ##########".format("Generated TIS with upstream start"))
print()
print("########## {:^50s} ##########".format("Start generating synthetic dataset with only dowstream stop"))
#Dataset with only stop codon
WriteTIS(rows=27102,posFile='arabTIS_wstop.pos',negFile='arabTIS_wstop.neg',startcheck=False,donor=False)
print("########## {:^50} ##########".format("Generated TIS with stop"))
print()
print("########## {:^50s} ##########".format("Start generating synthetic dataset with only donor"))
#Dataset with only donor splice site
WriteTIS(rows=27102,posFile='arabTIS_wdonor.pos',negFile='arabTIS_wdonor.neg',startcheck=False,stopcheck=False)
print("########## {:^50} ##########".format("Generated TIS with only donor"))
print()
print("########## {:^50s} ##########".format("Start generating synthetic dataset with only codon"))
#Dataset with only codon
c = ['codonUp_pos.txt','codonDown_pos.txt','codonUp_neg.txt','codonDown_neg.txt']
WriteTIS(rows=27102,posFile='arabTIS_wcodon.pos',negFile='arabTIS_wcodon.neg',startcheck=False,stopcheck=False,donor=False,CodonFile=c)
print("########## {:^50} ##########".format("Generated TIS with only codon"))
print()
print("########## {:^50s} ##########".format("Start generating synthetic dataset with only nucleotide"))
#Dataset with only nucleotide frequnecy
WriteTIS(rows=27102,posFile='arabTIS_wnuc.pos',negFile='arabTIS_wnuc.neg',startcheck=False,stopcheck=False,donor=False)
print("########## {:^50} ##########".format("Generated TIS with only nucleotide"))
print()
