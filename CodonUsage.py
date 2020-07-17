'''
July 13, 2020
'''

def CodonUsage(inputfile, upfile, downfile, positive):
    import json

    lines = [l for l in open(inputfile).readlines()]

    dna = ['A','C','T','G']
    upstream = dict()
    downstream = dict()
    for line in lines:
        #149->147:2, 148->147:298, 149 ATG
        line = line.strip()
        line = line[2:149]+line[152:299]
        for n in range(0,len(line),3):
            codon = line[n:n+3]
            if line[n] not in dna or line[n+1] not in dna or line[n+2] not in dna:
                continue
            if n <= 146:
                if codon in upstream:
                    upstream[codon] += 1
                else:
                    upstream[codon] = 1
            else:
                if codon in downstream:
                    downstream[codon] += 1
                else:
                    downstream[codon] = 1

    if positive == True:
        upstream.pop('ATG', None)
    stop = ['TAA','TAG','TGA']
    for j in stop:
        downstream.pop(j, None)

    with open(upfile, 'w') as up:
        json.dump(upstream, up, indent=4)
    with open(downfile, 'w') as down:
        json.dump(downstream, down, indent=4)

CodonUsage('at_neg_dic2013.txt','codonUp_neg.txt','codonDown_neg.txt',False)
