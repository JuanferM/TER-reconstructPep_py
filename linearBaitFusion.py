import csv
import sys
import pandas as pd
import matplotlib.pyplot as plt
from math import exp, trunc
from ctypes import ArgumentError

if len(sys.argv) != 3:
    raise FileNotFoundError('File not found. Please provide the stats filename\
 and the path to the mass table.')
infilename, massfilename = sys.argv[1], sys.argv[2]

# ---------------- PARAMETERS -----------------
verbose = False
onlythisbait = ""
minNumBaits = 1
maxNumBaits = float('inf') # included
reconstruct = 0.8
uncertainty, trace = 0.01, 1.0
valid, probation, invalid = 4, 1, 0
canreverse, reconstructFromBoth = True, True
# ---------------------------------------------

histo, baits, massTable = {}, {}, {}
nBait, nBaitOne, numBait = 0, 0, 0
solvedbaits, unsolvedbaits, nobaits = 0, 0, 0
massDispCount, minCount, maxCount, moyCount = 0, 0, 0, 0
stopcount, charcount, totalchar, totalInBM = 0, 0, 0, 0
stopAA = "KR"
NUMS = "-.0123456789"
mono = {"A" : 71.03,
        "R" : 156.10,
        "N" : 114.04,
        "D" : 115.02,
        "C" : 103.00,
        "E" : 129.04,
        "Q" : 128.05,
        "G" : 57.02,
        "H" : 137.05,
        "I" : 113.08,
        "L" : 113.08,
        "K" : 128.09,
        "M" : 131.04,
        "F" : 147.06,
        "P" : 97.05,
        "S" : 87.03,
        "T" : 101.04,
        "W" : 186.07,
        "Y" : 163.06,
        "V" : 99.06,
        "U" : 150.95
        }

# truncate to `decimals` decimals
def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return trunc(number)

    factor = 10.0 ** decimals
    return trunc(number * factor) / factor

# Compare bait with fusedBait
def compare(bait, fusedBait):
    same, i = 0, 0
    lenBait, lenFusedBait = len(bait), len(fusedBait)
    while i < lenBait and i < lenFusedBait:
        same += 1 if bait[i] == fusedBait[i] else 0
        i += 1
    # Return True if the sequences are the same (false otherwise)
    # Also return the number of equal characters
    return (same == lenBait and same == lenFusedBait), same

# Compute score for a baitModel
def scoreBM(stats):
    _, _, _, LS, _, GSC, GMC, GUM = stats
    return exp(LS/(10+10*GSC+100*GMC+1000*GUM))

# ------------- READING STATS FILE ------------
print("Reading stats file...")
with open(infilename, 'r') as file:
    i, rows = 1, file.read().split('\n')
    nBait, nBaitOne, massDispCount, minCount, maxCount, moyCount = rows[0].split()
    moyCount = float(moyCount)
    t = (nBait, nBaitOne, massDispCount, minCount, maxCount)
    nBait, nBaitOne, massDispCount, minCount, maxCount = map(int, t)
    while i < len(rows)-1:
        bait, nBaitModel, meanMass, sdMass = rows[i].split()
        nBaitModel = int(nBaitModel)
        baitStats = (float(meanMass), float(sdMass))
        baitModels, baitModelsStats = [], []
        for i in range(i+1, i+1+nBaitModel):
            baitModel, SPC, SGscore, baitMass, LS, nMass, GSC, GMC, GUM = rows[i].split()
            baitModels.append(baitModel)
            st = [int(SPC), float(SGscore), float(baitMass)]
            st = st + list(map(int, [LS, nMass, GSC, GMC, GUM]))
            baitModelsStats.append(st)
        i += 1
        if minNumBaits < nBaitModel <= maxNumBaits:
            numBait += 1
            baits[bait] = (baitModels, baitStats, baitModelsStats)
print("Done\n")
# ---------------------------------------------

# ------------- READ MASS TABLE ---------------
if massfilename != "":
    print("Reading mass table...")
    with open(massfilename, 'r') as file:
        reader = csv.reader(file, delimiter=',')
        numrow = 0
        for row in reader:
            if numrow != 0:
                massTable[row[0]] = []
                for combi in row[1:len(row)]:
                    if combi != '':
                        massTable[row[0]].append(combi)
            numrow += 1
    print("Done\n")
# ---------------------------------------------

# ------------ FUSION BAIT MODELS -------------
print("Verbose                                  : ", verbose)
print("Trace                                    : ", trace, " Da")
print("Uncertainty                              : ", uncertainty, " Da")
print("# of baits                               : ", numBait)
print("# of baits in stats file                 : ", nBait)
print("# of baits with at least one bait model  : ", nBaitOne)
print("# of baits with mass dispersion > 1.0 Da : ", massDispCount)
print("Min # of baitModels                      : ", minCount)
print("Max # of baitModels                      : ", maxCount)
print("Avg # of baitModels                      : ", moyCount)

for bait, data in baits.items():
    startchar, stopchar = '[', ']'
    keepgoing, stopped, reverse = True, False, False
    fusedBait, befRevBait, candidate, = "", "", "_"
    baitModels, baitStats, baitModelsStats = data
    lenBaitModels = len(baitModels)
    masses = [0.0 for _ in range(lenBaitModels)]
    indices = [0 for _ in range(lenBaitModels)]
    validation = [valid for _ in range(lenBaitModels)]
    if onlythisbait != "" and bait != onlythisbait:
        continue
    if verbose:
        print("\nBait   : ", bait)

    # Find sequence
    while (candidate not in stopAA or reverse) and keepgoing:
        # STEP1 : Elect candidate
        candidate = "_"
        candidates, scoresA, scoresB = {}, {}, {}
        for i in range(lenBaitModels):
            frommass, c = False, "_"

            if 0 <= indices[i] < len(baitModels[i]) and validation[i] != invalid:
                c = baitModels[i][indices[i]]

                # check if not mass
                if c == startchar:
                    indices[i] += -1 if reverse else 1
                    currentMass, c = "", baitModels[i][indices[i]]
                    while c != stopchar and c in NUMS:
                        currentMass = c+currentMass if reverse else currentMass+c
                        indices[i] += -1 if reverse else 1
                        c = baitModels[i][indices[i]]
                    c = '_'
                    masses[i] += truncate(float(currentMass), 2)

                # If there is enough excess mass from previous iteration
                if abs(masses[i]) >= trace:
                    j, mass, ncombi = 0, truncate(masses[i], 2), -1
                    while 0 <= j < 3:
                        mass = truncate(mass + (-1)*(j%2)*uncertainty*j, 2)
                        j += 1
                        if str(abs(mass)) in massTable:
                            j, ncombi = -abs(mass), len(massTable[str(abs(mass))])

                    if ncombi == 1:
                        combi = massTable[str(abs(j))][0]
                        if len(combi) == 1: # excess mass is an amino acid
                            c, frommass = combi[0], True
                else:
                    masses[i] = 0.0

            if c not in "_[]":
                if c not in candidates:
                    candidates[c], scoresA[c], scoresB[c] = 0, 0, 0
                candidates[c] += 0.25 if frommass else 1
                scoresA[c] += validation[i]
                scoresB[c] += validation[i]*scoreBM(baitModelsStats[i])

        mostpresent, doubt = 0, True
        for k, v in candidates.items():
            doubt = False
            if v > mostpresent:
                candidate, mostpresent = k, v
            elif v == mostpresent:
                if scoresA[k] > scoresA[candidate]:
                    candidate, mostpresent = k, v
                elif scoresB[k] > scoresB[candidate]:
                    candidate, mostpresent = k, v
                else:
                    sameA = scoresA[k] == scoresA[candidate]
                    sameB = scoresB[k] == scoresB[candidate]
                    doubt = sameA and sameB
        # No clear candidate
        if candidate not in mono or mostpresent == 0 or doubt:
            j = 0
            # Check if we reached the last character in reverse
            # If so then we're done and we may quit
            while canreverse and j >= 0 and j < lenBaitModels:
                # print(indices[j], validation[j])
                if indices[j] <= 0 or validation[j] == invalid:
                    j += 1
                else:
                    j = -1
            if canreverse and j == lenBaitModels:
                keepgoing = False
                break
            elif canreverse and not reverse:
                befRevBait = fusedBait
                reverse, fusedBait = True, ""
                startchar, stopchar = stopchar, startchar
                masses = [0.0 for _ in range(lenBaitModels)]
                indices = [len(bm)-1 for bm in baitModels]
                validation = [valid for _ in range(lenBaitModels)]
                continue
            else:
                stopcount += 1
                stopped, keepgoing = reconstructFromBoth, False
                break
        fusedBait += candidate # clear candidate

        # STEP2 : validate baitModels
        for i in range(lenBaitModels):
            if 0 <= indices[i] < len(baitModels[i]):
                c = baitModels[i][indices[i]]
                # Come back (previously invalidated but might match now)
                if c == candidate and abs(masses[i]) < trace:
                    if validation[i] == invalid:
                        validation[i] = probation
                else:
                    # excess mass...
                    if masses[i] >= trace:
                        if masses[i] - mono[candidate] >= 0:
                            # we can substract mass of candidate from excess mass
                            masses[i] -= mono[candidate]
                            masses[i] = truncate(masses[i], 2)
                        else:
                            # if excess mass corresponds to an amino acid then we
                            # put the baitModel in probation. The baitModel will be
                            # invalidated otherwise.
                            j, mass, ncombi = 0, truncate(masses[i], 2), -1
                            while j < 3 and j != -1:
                                mass = truncate(mass + (-1)*(j%2)*uncertainty*j, 2)
                                j += 1
                                if str(abs(mass)) in massTable:
                                    j, ncombi = -1, len(massTable[str(abs(mass))])
                            if ncombi == -1:
                                validation[i] = invalid
                    elif c in mono and validation[i] == valid: # no match but it's an amino acid
                        validation[i] = probation
                    else:
                        validation[i] = invalid

                if validation[i] != invalid and abs(masses[i]) < trace:
                    indices[i] += -1 if reverse else 1

    # STEP3: determine best sequence computed
    if reverse:
        if stopped:
            fusedBait = befRevBait + fusedBait[::-1]
        else:
            cond = len(fusedBait) > len(befRevBait)
            fusedBait = fusedBait[::-1] if cond else befRevBait

    # Some stats
    inBM = bait in baitModels
    totalInBM += 1 if inBM else 0
    fun = lambda l1, l2: min(l1, l2)/max(l1, l2) >= reconstruct
    rec = fun(len(bait), len(fusedBait))
    isequal, numMatch = compare(bait, fusedBait)

    if lenBaitModels not in histo:
        histo[lenBaitModels] = [0, 0, 0, 0, 0, 0, 0, 0]
    # Found bait data
    solvedbaits += 1 if isequal else 0
    histo[lenBaitModels][0] += 1 if isequal else 0
    histo[lenBaitModels][1] += 1 if isequal and inBM else 0
    # Not found but fusedBait is not empty
    cond = not isequal and rec
    unsolvedbaits += 1 if cond else 0
    histo[lenBaitModels][2] += 1 if cond else 0
    histo[lenBaitModels][3] += 1 if cond and inBM else 0
    # fusedBait is empty
    cond = not isequal and not rec
    nobaits += 1 if cond else 0
    histo[lenBaitModels][4] += 1 if cond else 0
    histo[lenBaitModels][5] += 1 if cond and inBM else 0
    # Totals
    histo[lenBaitModels][6] += 1 if inBM else 0
    histo[lenBaitModels][7] += 1
    charcount += numMatch
    totalchar += len(bait)
    if verbose:
        print("Fusion : ", fusedBait)
# --------------------------------------------

# ------------------ PLOTS -------------------
options = ""
if canreverse and not reconstructFromBoth:
    options = " (plus inverse)"
elif canreverse and reconstructFromBoth:
    options = " (deux sens)"
lableg = ["reconstitués", "reconstitués mais erronés", "non reconstitués"]
lengthBaitModels, tuples = tuple(zip(*[t[0] for t in sorted(zip(histo.items()))]))

propFound = list(map(lambda t: t[0]/t[7], tuples))
propNotEqual = list(map(lambda t: t[2]/t[7], tuples))
propNoSeq = list(map(lambda t: t[4]/t[7], tuples))

plt.bar(lengthBaitModels, propFound, width=1, color='g')
plt.bar(lengthBaitModels, propNotEqual, width=1, color='y', bottom=propFound)
plt.bar(lengthBaitModels, propNoSeq, width=1, color='r',
        bottom=[x + y for x, y in zip(propFound, propNotEqual)])
plt.legend(lableg, loc=4)
plt.title("Proportion de baits retrouvés selon le nombre de baitModels")
plt.xlabel("Nombre de baitModels")
plt.savefig("proportion baits reconstitues{}.png".format(options))
# plt.show()

propFoundInBM = list(map(lambda t: 0 if not t[6] else t[1]/t[6], tuples))
propNotEqualInBM = list(map(lambda t: 0 if not t[6] else t[3]/t[6], tuples))
propNoSeqInBM = list(map(lambda t: 0 if not t[6] else t[5]/t[6], tuples))

plt.clf()
plt.bar(lengthBaitModels, propFoundInBM, width=1, color='g')
plt.bar(lengthBaitModels, propNotEqualInBM, width=1, color='y', bottom=propFoundInBM)
plt.bar(lengthBaitModels, propNoSeqInBM, width=1, color='r',
        bottom=[x + y for x, y in zip(propFoundInBM, propNotEqualInBM)])
plt.legend(lableg, loc=4)
plt.title("Proportion de baits retrouvés selon le nombre de baitModels\n\
et où la séquence du bait fait partie des baitModels")
plt.xlabel("Nombre de baitModels")
plt.savefig("proportion baits reconstitues inclus dans baitmodels{}.png".format(options))
# plt.show()

propFoundInBM = list(map(lambda t: t[1]/t[7], tuples))
propNotEqualInBM = list(map(lambda t: t[3]/t[7], tuples))
propNoSeqInBM = list(map(lambda t: t[5]/t[7], tuples))

plt.clf()
plt_1 = plt.figure(figsize=(6.4, 5.2))
plt.bar(lengthBaitModels, propFoundInBM, width=1, color='g')
plt.bar(lengthBaitModels, propNotEqualInBM, width=1, color='y', bottom=propFoundInBM)
plt.bar(lengthBaitModels, propNoSeqInBM, width=1, color='r',
        bottom=[x + y for x, y in zip(propFoundInBM, propNotEqualInBM)])
plt.legend(lableg, loc=4)
plt.title("Proportion de baits retrouvés selon le nombre de baitModels\n\
et où la séquence du bait fait partie des baitModels\npar rapport à l'ensemble \
des baits retrouvés")
plt.xlabel("Nombre de baitModels")
plt.savefig("proportion baits reconstitues inclus dans baitmodels ensemble {}.png".format(options))
# plt.show()

# propBaitInBMFound = list(map(lambda t: t[6]/t[7], tuples))
# plt.bar(lengthBaitModels, propBaitInBMFound, width=1)
# plt.title("Proportion de baits selon le nombre de baitModels et\n\
# où la séquence du bait fait partie des baitModels")
# plt.xlabel("Nombre de baitModels")
# plt.show()
# --------------------------------------------

# ---------------- RESULTS -------------------
print("# of bait sequence incl. in bait models  : ", totalInBM)
print("\nSolved baits\t\t : {} / {} ({:.2f} %)".format(solvedbaits,
                                                     numBait,
                                                     (solvedbaits/numBait)*100))
print("Close enough baits\t : {} / {} ({:.2f} %)".format(unsolvedbaits,
                                                     numBait,
                                                     (unsolvedbaits/numBait)*100))
print("Wrong baits\t\t : {} / {} ({:.2f} %)".format(nobaits,
                                                     numBait,
                                                     (nobaits/numBait)*100))
print("# of matching characters : {} / {} ({:.2f} %)".format(charcount,
                                                             totalchar,
                                                             (charcount/totalchar)*100))
print("# of stopped resolution  : {} / {} ({:.2f} %)".format(stopcount,
                                                            numBait,
                                                            (stopcount/numBait)*100))
# --------------------------------------------


