import sys
from functions import *
from math import exp, trunc

# ------------ ARGUMENTS CHECK ----------------
if len(sys.argv) != 3:
    raise FileNotFoundError('File not found. Please provide the stats filename\
 and the path to the mass table.')
infilename, massfilename = sys.argv[1], sys.argv[2]
# ---------------------------------------------

# ---------------- PARAMETERS -----------------
verbose = False
fulltable = True
onlythisbait = ""
minNumBaits = 1
maxNumBaits = float('inf') # included
uncertainty, trace = 0.01, 1
valid, probation, invalid = 4, 1, 0
canreverse, reconstructFromBoth = True, True
cansimplify, simplifyBothWays = True, True
# ---------------------------------------------

# ---------- VARIABLES DEFINITION -------------
results = [0] * 21
solvedbaits, totalBaitInBM = 0, 0
totalBait, totalBaitWithOneBM, numBait = 0, 0, 0
numBaitWithMassDispersion, minBMcount, maxBMcount, meanBMcount = 0, 0, 0, 0
histo, baits, massTable, resultsPerBM = {}, {}, {}, {}
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
# ---------------------------------------------

# ------------- READ MASS TABLE ---------------
print("Reading mass table... ", end='')
massTable = readMassTable(massfilename)
print("Done")
# ---------------------------------------------

# ------------- READING STATS FILE ------------
print("Reading stats file... ", end='')
baits, totalBait, totalBaitWithOneBM, numBait, totalBaitInBM, numBaitWithMassDispersion,\
minBMcount, maxBMcount, meanBMcount = readStatsFile(infilename, minNumBaits, maxNumBaits)
print("Done\n")
# ---------------------------------------------

# --------------- PRINT STATS -----------------
printStats(verbose, trace, uncertainty, numBait, totalBait, totalBaitWithOneBM,
           numBaitWithMassDispersion, minBMcount, maxBMcount, meanBMcount, totalBaitInBM)
# ---------------------------------------------

# ------------ FUSION BAIT MODELS -------------
iteration = 1 # Iteration counter to print progress bar
for bait, data in baits.items():
    wholebaitmodel = False
    startchar, stopchar = '[', ']'
    keepgoing, stopped, reverse = True, False, False
    fusedBait, befRevBait, candidate, = "", "", "_"
    baitModels, baitStats, baitModelsStats = data
    lenBaitModels = len(baitModels)
    originalBaitModels = []
    masses = [0.0 for _ in range(lenBaitModels)]
    indices = [0 for _ in range(lenBaitModels)]
    validation = [valid for _ in range(lenBaitModels)]

    # If onlythisbait is defined then we only run the algorithm on onlythisbait
    # if it is found
    if onlythisbait != "" and bait != onlythisbait:
        continue
    if verbose:
        print("\nBait   : ", bait)

    # STEP0 : simplify baitmodels
    if cansimplify:
        originalBaitModels = baitModels.copy()
        baitModels = simplifyBM(mono, baitModels, baitModelsStats, massTable, uncertainty)

    # Find sequence
    while (candidate not in stopAA or reverse) and keepgoing:
        # STEP1 : Elect candidate
        candidate = "_"
        candidates, scoresA, scoresB = {}, {}, {}
        for i in range(lenBaitModels):
            if baitModelsStats[i][4] == 0:
                if baitModels[i] != bait:
                    wholebaitmodel = True

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
                if masses[i] < -trace:
                    validation[i] = invalid
                    continue
                elif abs(masses[i]) >= trace:
                    j, mass, ncombi = 0, truncate(masses[i], 2), -1
                    while 0 <= j < 3:
                        mass = truncate(mass + (-1)*(j%2)*uncertainty*j, 2)
                        j += 1
                        if str(abs(mass)) in massTable:
                            j, ncombi = -abs(mass), len(massTable[str(abs(mass))])

                    if ncombi == 1:
                        combi = massTable[str(-j)][0]
                        if len(combi) == 1:
                            # excess mass is an amino acid so c equals that
                            # amino acid
                            c, frommass = combi[0], True
                else:
                    masses[i] = 0.0

            # If c is an amino acid then we can use it as candidate
            if c not in "_[]":
                if c not in candidates:
                    candidates[c], scoresA[c], scoresB[c] = 0, 0, 0
                candidates[c] += 0.25 if frommass else 1
                scoresA[c] += validation[i]
                scoresB[c] += validation[i]*scoreBM(baitModelsStats[i])

        # Determine who is the elected candidate
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
        # If there is no clear elected candidate
        if candidate not in mono or mostpresent == 0 or doubt:
            j = 0
            # Check if we reached the last character in reverse
            # If so then we're done and we may quit
            while canreverse and j >= 0 and j < lenBaitModels:
                if indices[j] <= 0 or validation[j] == invalid:
                    j += 1
                else:
                    j = -1
            if canreverse and j == lenBaitModels:
                keepgoing = False # we're done
                break
            elif canreverse and not reverse:
                # the first pass didn't work out so we start the reverse pass
                befRevBait = fusedBait
                reverse, fusedBait = True, ""
                startchar, stopchar = stopchar, startchar
                masses = [0.0 for _ in range(lenBaitModels)]
                indices = [len(bm)-1 for bm in baitModels]
                validation = [valid for _ in range(lenBaitModels)]
                if cansimplify and not simplifyBothWays:
                    baitModels = originalBaitModels
                continue
            else:
                # the reverse pass (2nd pass) didn't work out so the algorithm
                # is stopped
                stopped, keepgoing = True, False
                break
        # there is a clear elected candidate so we add it to the sequence
        fusedBait += candidate

        # STEP2 : validate baitModels
        for i in range(lenBaitModels):
            if 0 <= indices[i] < len(baitModels[i]):
                c = baitModels[i][indices[i]]
                # Come back (previously invalidated baitModel but it might match now)
                if c == candidate and abs(masses[i]) <= trace:
                    if validation[i] == invalid:
                        validation[i] = probation
                else:
                    # excess mass...
                    if masses[i] >= trace:
                        # we can substract mass of candidate from excess mass
                        if masses[i] - mono[candidate] >= -trace:
                            masses[i] -= mono[candidate]
                            masses[i] = truncate(masses[i], 2)
                        else: # we can't substract mass of candidate so baitModel is invalidated
                            indices[i] = -1
                            validation[i] = invalid
                    elif c in mono and validation[i] == valid: # no match but it's an amino acid
                        validation[i] = probation
                    else:
                        validation[i] = invalid

                if validation[i] != invalid and abs(masses[i]) < trace:
                    indices[i] += -1 if reverse else 1

    # STEP3: determine best sequence computed
    bothWays = reverse and stopped and reconstructFromBoth
    if reverse:
        if bothWays:
            # Compute the mass M1 of the concatenation of the sequences
            # If that mass is greater than the mean of the baitModels masses by
            # a margin (M2) then we remove an amino acid from the right
            # sequence until M1 == M2
            seqmass = getMass(mono, befRevBait) + getMass(mono, fusedBait)
            while seqmass > baitStats[0] + trace and fusedBait != "":
                seqmass -= mono[fusedBait[-1]]
                fusedBait = fusedBait.rstrip(fusedBait[-1])
            # If M1 if less than M2 by at least m(G) (the mass of the smallest
            # amino acid) then we insert an '-' meaning that there is a gap of
            # unknown mass in the sequence
            if seqmass < baitStats[0] - mono['G'] + trace:
                fusedBait = befRevBait + '-' + fusedBait[::-1]
            else:
                fusedBait = befRevBait + fusedBait[::-1]
        else:
            # If concatenation is not ON then we return the longest sequence
            # between the two passes
            cond = len(fusedBait) > len(befRevBait)
            fusedBait = fusedBait[::-1] if cond else befRevBait

    # Print result if verbose else print progress bar
    if verbose:
        print("Fusion : ", fusedBait)
    else:
        if iteration == 1:
            print()
        printProgressBar(solvedbaits, iteration, numBait, prefix = 'Progress:', suffix
             = 'Solved')

    # Some stats
    lenbait = len(bait)
    inBM = bait in baitModels
    isequal, numMatch = compare(bait, fusedBait)
    if lenBaitModels not in resultsPerBM:
        resultsPerBM[lenBaitModels] = [0] * 11
    if not stopped:
        if isequal:
            resultsPerBM[lenBaitModels][0] += 1
            resultsPerBM[lenBaitModels][4] += 1
            results[0:3] = fillResults(results[0:3], inBM, wholebaitmodel)
        else:
            ratio = numMatch/lenbait
            resultsPerBM[lenBaitModels][1] += 1
            resultsPerBM[lenBaitModels][5] += 1 if ratio >= 0.8 else 0
            resultsPerBM[lenBaitModels][6] += 1 if ratio < 0.8 else 0
            results[3:6] = fillResults(results[3:6], inBM, wholebaitmodel)
    else:
        if isequal:
            resultsPerBM[lenBaitModels][2] += 1
            resultsPerBM[lenBaitModels][4] += 1
            results[6:9] = fillResults(results[6:9], inBM, wholebaitmodel)
        else:
            ratio = numMatch/lenbait
            resultsPerBM[lenBaitModels][3] += 1
            resultsPerBM[lenBaitModels][5] += 1 if ratio >= 0.8 else 0
            resultsPerBM[lenBaitModels][6] += 1 if ratio < 0.8 else 0
            if ratio >= 0.8:
                results[9:12] = fillResults(results[9:12], inBM, wholebaitmodel)
            elif ratio >= 0.5:
                results[12:15] = fillResults(results[12:15], inBM, wholebaitmodel)
            elif ratio >= 0.3:
                results[15:18] = fillResults(results[15:18], inBM, wholebaitmodel)
            else:
                results[18:21] = fillResults(results[18:21], inBM, wholebaitmodel)
    if not inBM:
        if isequal:
            resultsPerBM[lenBaitModels][7] += 1
        else:
            resultsPerBM[lenBaitModels][8] += 1
    else:
        if isequal:
            resultsPerBM[lenBaitModels][9] += 1
        else:
            resultsPerBM[lenBaitModels][10] += 1

    # Totals
    solvedbaits += 1 if isequal else 0
    iteration += 1
# ---------------------------------------------

# ----------------- RESULTS -------------------
if solvedbaits != 0:
    printResults(solvedbaits, numBait, results, fulltable)
else:
    print("\nNO SOLUTION!")
# ---------------------------------------------

# ------------------- PLOTS -------------------
options = " (" if canreverse or cansimplify else ""
if canreverse:
    if reconstructFromBoth:
        options += "deux sens"
    else:
        options += "plus inverse"
options += " + " if canreverse and cansimplify else ""
if cansimplify:
    if simplifyBothWays:
        options += "simplif deux sens"
    else:
        options += "simplif gauche-droite"
options += ")" if canreverse or cansimplify else ""

if solvedbaits != 0:
    lengthBaitModels, respBM = tuple(zip(*[t[0] for t in sorted(zip(resultsPerBM.items()))]))
    plotResults(options, lengthBaitModels, respBM)
# ---------------------------------------------
