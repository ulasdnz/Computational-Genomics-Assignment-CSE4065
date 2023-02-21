from copy import copy
import random
import math
import time


f = open("dnabase.txt","r")

KMERS= []
Dna = []

for x in f:
  x = x.replace("\n","")
  x = list(x)
  Dna.append(x)

f.close()


def getLetter(numb):
   if numb % 4 == 0: return "A"
   if numb % 4 == 1: return "C"
   if numb % 4 == 2: return "G"
   if numb % 4 == 3: return "T"


def getAllKmerComb(kMerNumb):
    loop = pow(4,kMerNumb)
    arr = list(range(kMerNumb))
    for i in range(loop):
        for x in range(kMerNumb):
            arr[x] = getLetter(math.floor(i / pow(4, abs(x+1-kMerNumb))))
        KMERS.append(arr.copy())


def getProfile(Motifs, kMerNumb, isGibbs):
    profile={"A":[],"C":[],"G":[],"T":[]}
    rowLength = len(Motifs)
    divide = len(Motifs)
    if isGibbs:
        divide = divide+4
    for i in range(kMerNumb):
        a=0;c=0;g=0;t=0
        if isGibbs:
            a=1;c=1;g=1;t=1
        for x in range(rowLength):
            if Motifs[x][i]=="A":
                a +=1
            elif Motifs[x][i]=="C":
                c+=1
            elif Motifs[x][i]=="G":
                g+=1
            elif Motifs[x][i]=="T":
                t += 1
        profile["A"].append([a/divide])
        profile["C"].append([c/divide])
        profile["G"].append([g/divide])
        profile["T"].append([t/divide])

    return profile


def getConsensus(Motif, kMerNumb):
    profile = getProfile(Motif, kMerNumb, False)
    ltrs = ["A","C","G","T"]
    consensus = []
    for x in range(kMerNumb):
        highestInColumn = 0
        letter = ""
        for y in range(4):
            if highestInColumn < profile[ltrs[y]][x][0]:
                highestInColumn = profile[ltrs[y]][x][0]
                if y==0:
                    letter = "A"
                elif y==1:
                    letter = "C"
                elif y==2:
                    letter = "G"
                elif y==3:
                    letter = "T"
        consensus.append(letter)
    return consensus


def getDistance(kMer, dnaStrings):
    result = 0
    motif = []
    columnLength = len(dnaStrings[0])
    kmerLength = len(kMer)
    for dna in dnaStrings:
        localDistance = 999
        localMotif = []
        for i in range(columnLength - kmerLength + 1):
            dst = 0;
            for indx, letter in enumerate(kMer):
                if letter != dna[i+indx]:
                    dst +=1
            if dst < localDistance:
                localDistance = dst
                localMotif = dna[i:i+kmerLength]
        result += localDistance
        motif.append(localMotif)

    return {"score":result,"motif":motif}


def getScore(Motif, kMerNumb):
    consensus = getConsensus(Motif, kMerNumb)
    return getDistance(consensus, Motif)["score"]


def MedianString(kMerNumb):
    getAllKmerComb(kMerNumb)
    bestKmer = KMERS[0]
    res = getDistance(bestKmer, Dna)
    bestDist = res["score"]
    bestMotif =  copy(res["motif"])
    for kMer in KMERS:
        newRes = getDistance(kMer, Dna)
        if newRes["score"] < bestDist:
            bestKmer = kMer
            bestDist = newRes["score"]
            bestMotif = copy(newRes["motif"])
    print(bestKmer)
    print(bestDist)
    print(bestMotif)



def selectRandomMotifs(dnaStrings, kMerNumb):
    motifs = []
    columnLength = len(dnaStrings[0])
    rowLength = len(dnaStrings)
    for i in range(rowLength):
        rand = random.randrange(0,columnLength-kMerNumb+1)
        sliced = dnaStrings[i][rand:rand+kMerNumb]
        motifs.append(copy(sliced))
    return motifs



def getMotif(profile, kMerNumb, dnaStrings):
    motif = []
    columnLength = len(dnaStrings[0])
    for dna in dnaStrings:
        optimum=0
        motifRow=0
        for x in range(columnLength-kMerNumb+1):
            probability = 1
            row = []
            for i in range(kMerNumb):
                probability = profile[dna[x+i]][i][0]*probability
                row.append(dna[x+i]) 
            if optimum < probability:
                optimum = probability
                motifRow = row

        motif.append(copy(motifRow))
    return motif



def RandomizedMotifSearch(kMerNumb):
    Motifs = selectRandomMotifs(Dna,kMerNumb)
    bestMotifs = copy(Motifs)
    while True:
        profile = getProfile(Motifs, kMerNumb, False)
        Motifs = getMotif(profile, kMerNumb, Dna)
        if getScore(Motifs,kMerNumb) < getScore(bestMotifs,kMerNumb):
            bestMotifs = copy(Motifs)
        else:
            return {"score":getScore(bestMotifs,kMerNumb), "consensus":getConsensus(bestMotifs,kMerNumb), "motif":bestMotifs} 




def getMotifGibbs(profile, kMerNumb, string):
    motif = {
        "motifs":[],
        "probabilities":[]
    }
    columnLength = len(string)
    for x in range(columnLength-kMerNumb+1):
        probability = 1
        row = []
        for i in range(kMerNumb):
            probability = profile[string[x+i]][i][0]*probability
            row.append(string[x+i])

        motif["motifs"].append(copy(row))
        motif["probabilities"].append(probability)
    return motif



def GibbsSampler(kMerNumb):
    Motifs = selectRandomMotifs(Dna,kMerNumb)
    bestMotifs = copy(Motifs)
    bestScore = getScore(bestMotifs,kMerNumb)
    iterationNumber = 0
    while True:
        rnd = random.randrange(0,10)
        Motifs.pop(rnd)
        profile = getProfile(Motifs, kMerNumb,True)
        MotifRemovedRow = getMotifGibbs(profile, kMerNumb, Dna[rnd])
        motifRollDie = random.choices(MotifRemovedRow["motifs"], weights= MotifRemovedRow["probabilities"],k=1)[0]
        Motifs.insert(rnd,motifRollDie)
        currentScore = getScore(Motifs,kMerNumb) 
        iterationNumber +=1
        if currentScore <= bestScore:
            bestMotifs = copy(Motifs)
            bestScore = currentScore
            iterationNumber = 0        
        if iterationNumber == 50 and  bestScore<=currentScore:
            return {"score":getScore(bestMotifs,kMerNumb), "consensus":getConsensus(bestMotifs,kMerNumb), "motif":bestMotifs}





inp = input("Enter the k value:\n")
kValue = int(inp)

print("\nGibbs Sampler Algorithm is running!")
start = time.time()
resultGibbsSampler = GibbsSampler(kValue)
print(resultGibbsSampler["score"])
print(resultGibbsSampler["consensus"])
print(resultGibbsSampler["motif"])
end = time.time()
print("Total time elapsed: ",(end - start))



print("\n\n\nRandomized Motif Search Algorithm is running!")
start = time.time()
randomizedResult = RandomizedMotifSearch(kValue)
end = time.time()
print("Total time elapsed: ",(end - start))
print("Randomized Motif Search Score: ",randomizedResult["score"])
print("Randomized Motif Search Consensus: ",randomizedResult["consensus"])
print("Randomized Motif Search Motif: ",randomizedResult["motif"])


print("\n\n\nMedian String Algorithm is running!")
start = time.time()
MedianString(kValue)
end = time.time()
print("Total time elapsed: ",(end - start))