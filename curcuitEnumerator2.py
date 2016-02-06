#-------------------------------------------------------------------------------
# Name:        circuitEnumerator
# Purpose:     enumerate all legal circuits in the 5-state machine architecture
# by performing depth-first search over a grammar which describes the circuits
#
# Author:      mquintin
#
# Created:     19/05/2015
# Copyright:   (c) mquintin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import operator, collections, copy

#maximum number of parts
nmax = 12 #default=12

#for when we need to iterate
parttypes = 20
parts = range(0,parttypes)
recs = [3,4,5] #recombinase sites
genes = range(6,parttypes)
runcounter = [0]

#because these aren't built in to python
def prod(factors):
    return reduce(operator.mul, factors, 1)
def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]


seqs = []



##ID	PartName
##0	BiTerminator
##1	UniTerminator
##2	Promoter
##3	RecA1
##4	RecA2
##5	RecB
##6	Gene1
##7	Gene2
##8	Gene3
##9	Gene4
##10	Gene5
##11	Gene6
##12	Gene7
##13	Gene8
##14	Gene9
##15	Gene10
##16	Gene11
##17	Gene12
##18	Gene13
##19	Gene14


#flip/delete recombinaseA sites
#IDs of interest are 3 & 4 (and -3,-4) for recA, 5 for recB
#There is guaranteed to be no overlap between 3&4
#There is guaranteed to be either 0 or 2 of each
#
#seq=[int]
#x= absolute value of part number to flip
#return:[int]
def permute(seq,x):
    s = copy.deepcopy(seq)
    #location of sites
    idx = [i for i,j in enumerate(s) if x==abs(j)]
    if(len(idx) == 2):
        if(s[idx[0]] + s[idx[1]] == 0): #opposite: invert
            subseq = s[idx[0]:idx[1]+1]
            subseq.reverse()
            s[idx[0]:idx[1]+1] = subseq
        else: #excise
            s[idx[0]:idx[1]+1] = []
    return s

def permuteA(seq):
    s = permute(seq,3)
    s = permute(s,4)
    return s

def permuteB(seq):
    s = permute(seq,5)
    return s

#increment the last element in the given circuit
#if it overflows, add an element to the end and start the count over
def getNextSeq(seq):
    runcounter[0] += 1
    if runcounter[0] == 500000:
        runcounter[0] = 0
        print(seq)
    if (seq[-1] <= (-1 * genes[1])): #these are all going to fail anyway
        seq[-1] = (-1 * genes[1]) + 1
        return seq
    pointer = len(seq)-1
    return nextSeq(seq,pointer)

def nextSeq(seq,pointer):
    if seq[pointer] >= parts[-1]:
        if pointer <= 0:
            #print(seq, pointer)
            #start at the lowest value seq of length +1
            newseq = [-parts[-1]]*(len(seq)+1)
            newseq[0] = -5
            return newseq
        else:
            return nextSeq(seq,pointer-1)
    else:
        seq[pointer]+=1 #increment at pointer
        #reset every value after the pointer
        seq[pointer+1:len(seq)]=[-parts[-1]]*(len(seq)-pointer-1)
        return seq

#enforces rules about part numbers
def validateFormat(seq):
    #each recombinase site is used 0 or 2 times
    pos = map(abs,seq)
    for r in recs: #recs = [3,4,5]: recombinase sites
        if pos.count(r) not in [0,2]:
            return 0
    #a gene and a promoter are present
    if pos.count(6) < 1:
        return 0
    if pos.count(2) <1:
        return 0
    #the circuit must use only the lowest numbered genes.
    #get the highest gene number, make sure all genes between
    #it and geneA (6) are present
    genesUsed = range(6,max(pos)) #doesn't include the max, we know it's used
    for g in genesUsed:
        if pos.count(g) <1:
            return 0
    #adjacency rules
    for i in range(0,len(seq)-1):
        #no adjacent identical parts
        if seq[i]==seq[i+1]:
            return 0
        #no adjacent opposing unidirectional terminators
        if sorted(seq[i:i+2]) == [-1,1]:
            return 0
    #recombinaseA sites are not improperly staggered, where the result will
    #depend on if RecA1 or RecA2 switches first
    #ex: 1,2,1,2  1,2,B,2,1,B
    if pos.count(3) == 2 & pos.count(4) == 2:
        x = [pos.index(3),len(pos - pos[::-1].index(3))] #first and last occurrence
        y = [pos.index(4),len(pos - pos[::-1].index(4))]
        xy = (x+y).sort
        if (x[0] == xy[0] & x[1] == xy[2]) | (x[0] == xy[1] & x[1] == xy[2]):
            return 0
        #in terms of RecB, the only invalid arrangements are x,y,B,y,x,B and its
        #reverse, or x,B,y,x,B,y. If the B sites are 2 away, make sure the sites
        #between them are matching targets (x,B,y,y,B,x  B,x,x,B,y,y)
        if pos.count(5) == 2:
            z = [pos.index(5),len(pos - pos[::-1].index(5))]
            if z[1] - z[0] == 2:
                if sorted([pos[z[0]+1],pos[z[0]+2]]) == [3,4]:
                    return 0
    return 1

#enforces rules about interaction of parts
#only checks left to right; must also validate the inverse of a sequence
#returns a list of integers >=1 iff the corresponding part didn't misbehave
#   >1 if the part us properly used
def validateExpression(seq):
    res = [1] * len(seq)
    for i in range(0,len(seq)):
        if abs(seq[i]) in[3,4,5]:
            res[i] *=2 #already know the count of recombinases is good
        if seq[i] == 2: #promoter
            express = [] #idx of genes that this promoter expresses
            for j in range(i+1,len(seq)):
                if seq[j] in [0,1]: #terminators
                    if len(express) > 0:
                        for k in ([i] + [j] + express):
                            res[k]*=2 #good score
                    break; #this promoter's done
                if seq[j] in genes:
                    express.append(j)
                if seq[j] == 2: #redundant promoter
                    if len(express) == 0:
                        break;
    return res;

#returns two lists of lists. The first is
def validatePermutation(seq):
    s = copy.deepcopy(seq)
    scores = [validateExpression(s)]
    imap = range(0,len(s)) #map indexes to positions in the score list
    #shuffle A
    recA1 = [i for i,j in enumerate(s) if 3==abs(j)]
    if len(recA1) == 2:
        imap = permuteIdx(imap,recA1[0],recA1[1],seq[recA1[0]]==seq[recA1[1]])
    recA2 = [i for i,j in enumerate(s) if 4==abs(j)]
    if len(recA2) == 2:
        imap = permuteIdx(imap,recA2[0],recA2[1],seq[recA2[0]]==seq[recA2[1]])
    if len(recA1 + recA2) > 0:
        val = validateExpression(permuteA(copy.deepcopy(seq)))
        for i in range(0,len(imap)):
            score = [1] * len(seq)
            score[imap[i]] *= val[i]
            scores.append(score)
    #shuffle B
    imap = range(0,len(s))
    recB = [i for i,j in enumerate(s) if 5==abs(j)]
    if len(recB) == 2:
        imap = permuteIdx(imap,recB[0],recB[1],seq[recB[0]]==seq[recB[1]])
        val = validateExpression(permuteB(copy.deepcopy(seq)))
        for i in range(0,len(imap)):
            score = [1] * len(seq)
            score[imap[i]] *= val[i]
            scores.append(score)
    if (len(recA1 + recA2) > 0) & (len(recB) > 0):
        #shuffle A,B
        imap = range(0,len(s))
        recA1 = [i for i,j in enumerate(s) if 3==abs(j)]
        if len(recA1) == 2:
            imap = permuteIdx(imap,recA1[0],recA1[1],seq[recA1[0]]==seq[recA1[1]])
        recA2 = [i for i,j in enumerate(s) if 4==abs(j)]
        if len(recA2) == 2:
            imap = permuteIdx(imap,recA2[0],recA2[1],seq[recA2[0]]==seq[recA2[1]])
        recB = [i for i,j in enumerate(s) if 5==abs(j)]
        imap = permuteIdx(imap,recB[0],recB[1],seq[recB[0]]==seq[recB[1]])
        val = validateExpression(permuteA(permuteB(copy.deepcopy(seq))))
        for i in range(0,len(imap)):
            score = [1] * len(seq)
            score[imap[i]] *= val[i]
            scores.append(score)
        #shuffle B,A
        imap = range(0,len(s))
        recB = [i for i,j in enumerate(s) if 5==abs(j)]
        imap = permuteIdx(imap,recB[0],recB[1],seq[recB[0]]==seq[recB[1]])
        recA1 = [i for i,j in enumerate(s) if 3==abs(j)]
        if len(recA1) == 2:
            imap = permuteIdx(imap,recA1[0],recA1[1],seq[recA1[0]]==seq[recA1[1]])
        recA2 = [i for i,j in enumerate(s) if 4==abs(j)]
        if len(recA2) == 2:
            imap = permuteIdx(imap,recA2[0],recA2[1],seq[recA2[0]]==seq[recA2[1]])
        val = validateExpression(permuteB(permuteA(copy.deepcopy(seq))))
        for i in range(0,len(imap)):
            score = [1] * len(seq)
            score[imap[i]] *= val[i]
            scores.append(score)
    #print(scores)
    return scores


#rearrange the given index list with recombination sites at the given positions
def permuteIdx(imap, left, right, match):
    i = copy.deepcopy(imap)
    if(not match): #opposite: invert
        subseq = i[left:right+1]
        subseq.reverse()
        i[left:right+1] = subseq
    else: #excise
        i[left:right+1] = []
    return i

def getUseScores(seq):
    s = copy.deepcopy(seq)


#check that the machine is legal in all states
def validate(seq):
    formatOK = validateFormat(seq)
    if formatOK == 0:
        return 0
    scores = []
    v0 = validatePermutation(seq)
    if not all([x>=1 for x in v0]):
        return 0
    pa = permuteA(seq)
    pb = permuteB(seq)
    pab = permuteB(pa)
    pba = permuteA(pb)
    va = validateExpression(pa)
    if not all([x>=1 for x in va]): #no part did something illegal here
        return 0
    vb = validateExpression(pb)
    if not all([x>=1 for x in vb]):
        return 0
    vab = validateExpression(permuteB(pa))
    if not all([x>=1 for x in vab]):
        return 0
    vba = validateExpression(permuteA(vb))
    if not all([x>=1 for x in vba]):
        return 0
    rev = map(lambda x : x*-1,seq[::-1]) #reverse seq and flip the signs
    v0r = validateExpression(rev)
    par = permuteA(rev)
    pbr = permuteB(rev)
    pabr = permuteB(par)
    pbar = permuteA(pbr)
    var = validateExpression(par)
    if not all([x>=1 for x in v0r]):
        return 0
    if not all([x>=1 for x in var]): #no part did something illegal here
        return 0
    vbr = validateExpression(pbr)
    if not all([x>=1 for x in vbr]):
        return 0
    vabr = validateExpression(permuteB(par))
    if not all([x>=1 for x in vabr]):
        return 0
    vbar = validateExpression(permuteA(vbr))
    if not all([x>=1 for x in vbar]):
        return 0
    #check that all parts are used
    allScores= validatePermutation(seq) #[v0,va,vb,vab,vba,v0r,var,vbr,vabr,vbar]
    total = [1] * len(seq)
    for score in allScores:
        total = [a * b for a,b in zip(total,score)] #a part will get >1 if it's
            #used properly once
    #print("total: ",total)
    #print("flat: ", flatten(total))
    if not all([x>1 for x in flatten(total)]):
        return 0
    return 1

#output a list of ints of length = # of genes, in order
#0 if not expressed, 1 if expressed
def truth(seq):
    res = [0] * len(genes)
    offset = genes[0] #part number of the first gene
    rev = map(lambda x : x*-1,seq[::-1]) #reverse seq and flip the signs
    for s in [seq,rev]:
        transcribing = False
        idx = 0
        while idx < len(s):
            if s[idx]==2:
                transcribing = True
            if s[idx] in [0,1]:
                transcribing = False
            if transcribing & (s[idx] in genes):
                res[s[idx]-offset] = 1
            idx += 1
    return res

#create/append line in sequence file to represent the machine
def writeSeq(seq,sfile):
    with open(sfile, 'a') as f:
        f.write("\t".join(map(str,seq)))
        f.write("\n")

#create append line in truth file
#each line contains nGenes*5 tab-delimited digits
#these correspond to each gene in numerical order for each orientation
#in the order [no change, recA, recB, recA/B, recB/A]
def writeTruth(seq,tfile):
    pa = permuteA(seq)
    pb = permuteB(seq)
    pab = permuteB(pa)
    pba = permuteA(pb)
    t = truth(seq) + truth(pa) + truth(pb) + truth(pab) + truth (pba)
    with open(tfile, 'a') as f:
        f.write("\t".join(map(str,t)))
        f.write("\n")



#build all the circuits up to the max length, print them if they validate
#shortest legal circuits start at length 2 (promoter+gene)
def generate(seq, seqfile, truthfile):
        pa = permuteA(seq)
        pb = permuteB(seq)
        pab = permuteB(pa)
        pba = permuteA(pb)
        if validate(seq):
            writeSeq(seq,seqfile)
            writeTruth(seq,truthfile)

def test():
    #seq = getNextSeq([1,-5,-6])
    #writeSeq([1,2,3,4],'testseq.txt')
    #writeTruth([0,1,0,1],'testtruth.txt')
    #print(seq)
    seq = [4,2,6,-4]
    print(validate(seq))
    print('seq: ', seq)
    pa = permuteA(seq)
    print('pa: ', pa)
    print('seq: ', seq)
    pb = permuteB(seq)
    print('pb: ', pb)
    print('seq: ', seq)
    pab = permuteB(pa)
    print('pab: ', pab)
    print('seq: ', seq)
    pba = permuteA(pb)
    print('pba: ', pba)
    print('seq: ', seq)
    print([pa,pb,pab,pba])
    print([validateFormat(seq),
        validatePermutation(seq),
        validatePermutation(pa),
        validatePermutation(pb),
        validatePermutation(pab),
        validatePermutation(pba),
        validate(seq)])
    print('validate expression:')
    print(validateExpression(seq))
    print(pab)
    print(validateExpression(pab))
    print(truth(seq))
    print(not all([x>1 for x in flatten(truth(seq))]))
    generate(seq,'testseq.txt','testtruth.txt')



def main():
    testing = 0
    if testing:
        test()
    else:
        #seq = [-1,-1*genes[0],-2] #the shortest circuit is 3 elements: pro/gene/term
        seq = [-5,-6,-19,-19,-19,-19,-19]
        while len(seq) <= nmax:
            generate(seq, "circuits3.txt", "truth3.txt")
            seq=getNextSeq(seq)


if __name__ == '__main__':
    main()