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

"""
def writeCodon():

    table = {('TTT','TTC'):'Phe', ('TTA','TTG','CTT','CTC','CTA','CTG'):'Leu', ('ATT','ATC','ATA'):'Ile',('ATG'):'Met',
             ('GTT','GTC','GTA','GTG'):'Val',('TCT','TCC','TCA','TCG','AGT','AGC'):'Ser',('CCT','CCC','CCA','CCG'):'Pro',
             ('ACT','ACC','ACA','ACG','TAT','TAC'):'Thr',('GCT','GCC','GCA','GCG'):'Ala',('TAT','TAC'):'Try',
             ('TAA','TAG','TGA'):'Stop',('CAT','CAC'):'His',('CAA','CAG'):'Gln',('AAT','AAC'):'Asn',('AAA','AAG'):'Lys',
             ('GAT','GAC'):'Asp',('GAA','GAG'):'Glu',('TGT','TGC'):'Cys',('CGT','CGC','CGA','CGG','AGA','AGG'):'Arg',
             ('GGT','GGC','GGA','GGG'):'Gly'}

    out = open(outputfile, 'w+')
    UpCodon = "Upstream \n"
    DownCodon = "Downstream \n"
    for k1, k2 in zip(upstream.keys(), downstream.keys()):
        #upstream[k1] = round(upstream[k1]/27102, 2)
        #downstream[k2] = round(downstream[k2]/27102, 2)
        for i in table.keys():
            if k1 in i:
                n1 = table[i]
            elif k2 in i:
                n2 = table[i]
        UpCodon += "{:^9}".format(n1)
        DownCodon += "{:^9}".format(n2)
    UpCodon += "\n"
    DownCodon += "\n"
    for k1, k2 in zip(upstream.keys(), downstream.keys()):
        UpCodon += "{:^9}".format(k1)
        DownCodon += "{:^9}".format(k2)
    UpCodon += "\n"
    DownCodon += "\n"
    for v1,v2 in zip(upstream.values(),downstream.values()):
        UpCodon += "{:^9}".format(str(v1))
        DownCodon += "{:^9}".format(str(v2))
    out.write(UpCodon + "\n" + DownCodon)
    out.close()
"""

CodonUsage('arabTIS_codon.pos','testUp.txt','testDown.txt',True)
#CodonUsage('at_neg_dic2013.txt','codonUp_neg.txt','codonDown_neg.txt',False)