"""
This function gives out the codon frequency of a file of sequences.
It takes a file of sequences, two files for output and a boolean to indicate if file is positive set or not.
It gives two files, each for upstream and downstream codon frequency.
"""
def Codon(inputfile, upfile, downfile, positive):
    import json

    lines = [l for l in open(inputfile).readlines()]
    u,d = 0,0       #counter for number of codons
    dna = ['A','C','T','G']
    upstream = dict()       #dictionary for upstream codon frequency
    downstream = dict()     #dictionary for downstream codon frequency
    for line in lines:
        line = line.strip()
        # only get sequence with codons in frame with central ATG
        if positive==False:
            line1 = line[2:149]
            line2 = line[152:299]
        #for positive set, same as negative, but does not take into account the codons for consensus
        else:
            line1 = line[2:137]
            line2 = line[164:299]

        #get each codon from the sequence
        for n in range(0,len(line1),3):
            codon = line1[n:n+3]
            #do not take it into account if the codon contains non DNA letters
            if line1[n] not in dna or line1[n+1] not in dna or line1[n+2] not in dna:
                continue
            u += 1
            if codon in upstream:
                upstream[codon] += 1
            else:
                upstream[codon] = 1

        for n in range(0, len(line2), 3):
            codon = line2[n:n + 3]
            # do not take it into account if the codon contains non DNA letters
            if line2[n] not in dna or line2[n + 1] not in dna or line2[n + 2] not in dna:
                continue
            d += 1
            if codon in downstream:
                downstream[codon] += 1
            else:
                downstream[codon] = 1

    #remove upstream start codons from positive data
    if positive == True:
        u -= upstream['ATG']
        upstream.pop('ATG', None)
    #remove all downstream stop codons
    stop = ['TAA','TAG','TGA']
    for j in stop:
        d -= downstream[j]
        downstream.pop(j, None)

    print("number of upstream codons: {:>10}".format(u))
    print("number of downstream codons: {:>10}".format(d))
    #divide by number of codons and multiply 1000 for all values
    upstream = {k: round(v/u*1000, 2) for k,v in upstream.items()}
    downstream = {k: round(v/d*1000, 2) for k, v in downstream.items()}

    with open(upfile, 'w') as up:
        json.dump(upstream, up, indent=4)
    with open(downfile, 'w') as down:
        json.dump(downstream, down, indent=4)
