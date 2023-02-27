import os
import sys
import tempfile
from functions import *
import subprocess as sbp
from math import trunc, ceil

if len(sys.argv) != 3:
    raise FileNotFoundError('File not found. Please provide the stats filename\
 and the path to the mass table.')
infilename, massfilename = sys.argv[1], sys.argv[2]

# ---------------- PARAMETERS -----------------
verbose = False
fulltable = False
onlythisbait = ""
minNumBaits = 1
maxNumBaits = float('inf') # included
uncertainty, trace = 0.01, 1
useMUSCLE, cansimplify = True, True
clustalopt = "-QUICKTREE -MATRIX=GONNET -GAPOPEN=5 -GAPEXT=1 -NOHGAP \
-NOWEIGHTS -CLUSTERING=UPGMA"
muscleopt = ""
# ---------------------------------------------

results = [0] * 21
charcount, totalchar, totalInBM = 0, 0, 0
nBait, nBaitOne, numBait, solvedbaits = 0, 0, 0, 0
histo, baits, massTable, resultsPerBM = {}, {}, {}, {}
massDispCount, minCount, maxCount, moyCount = 0, 0, 0, 0
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
clustalcmd = "./clustalw2 -ALIGN -QUIET -OUTPUT=FASTA"
musclecmd = "./muscle3.8"

# ------------- READING STATS FILE ------------
nBait, nBaitOne, numBait, totalInBM, baits = readStatsFile(infilename,
                                                           minNumBaits,
                                                           maxNumBaits)
# ---------------------------------------------

# ------------- READ MASS TABLE ---------------
massTable = readMassTable(massfilename)
# ---------------------------------------------

# --------------- PRINT STATS -----------------
printStats(verbose, trace, uncertainty, numBait, nBait, nBaitOne,
           massDispCount, minCount, maxCount, moyCount, totalInBM)
# ---------------------------------------------

# ------------ FUSION BAIT MODELS -------------
for bait, data in baits.items():
    wholebaitmodel = False
    keepgoing, stopped = True, False
    fusedBait, candidate = "", "_"
    baitModels, baitStats, baitModelsStats = data
    lenBaitModels = len(baitModels)
    indices = [0 for _ in range(lenBaitModels)]
    convertedBaitModels = []
    if onlythisbait != "" and bait != onlythisbait:
        continue
    if verbose:
        print("\nBait   : ", bait)

    # STEP0 : simplify baitmodels
    if cansimplify:
        baitModels = simplifyBM(baitModels, baitModelsStats, massTable,
                                uncertainty, NUMS)

    # STEP1 : convert baitmodels
    for i in range(lenBaitModels):
        converted, currentMass = "", ""
        for c in baitModels[i]:
            if c == ']' and currentMass != "":
                mass = float(currentMass)

                # If there is enough excess mass
                if abs(mass) >= trace:
                    j, mass, ncombi = 0, truncate(mass, 2), -1
                    while 0 <= j < 3:
                        mass = truncate(mass + (-1)*(j%2)*uncertainty*j, 2)
                        j += 1
                        if str(abs(mass)) in massTable:
                            j, ncombi = -abs(mass), len(massTable[str(abs(mass))])

                    if ncombi >= 1:
                        moylen = 0
                        for combi in massTable[str(-j)]:
                            lencombi = len(combi)
                            moylen += lencombi
                        converted += '-'*trunc(moylen/ncombi)
                    else:
                        converted += '-'*(ceil(-j/mono['G']))

                currentMass = ""
            elif c not in NUMS and c not in '[]':
                converted += c
            elif c in NUMS:
                currentMass += c
            else:
                pass
        convertedBaitModels.append(converted)

    # STEP2 : create FASTA file
    fd, path = tempfile.mkstemp()
    with open(path, 'wb') as tmp:
        for i in range(lenBaitModels):
            tmp.write(">s{}\n{}\n".format(i, convertedBaitModels[i]).encode())

    # STEP3 : run Clustal on the file and create an output file
    IN, OUT = path, path+".out"
    if useMUSCLE:
        CMD = musclecmd+" {} -in {} -out {}".format(muscleopt, IN, OUT)
    else:
        CMD = clustalcmd+" {} -INFILE={} -OUTFILE={}".format(clustalopt, IN, OUT)
    with open(os.devnull, 'wb') as nil:
        sbp.check_call(CMD.split(), stdout=nil, stderr=sbp.STDOUT)

    # STEP4 : read outputfile and determine sequence
    with open(OUT, 'r') as f:
        lines = f.read().split('\n')
        sequences = [""]*lenBaitModels
        idx = -1
        for line in lines:
            if '>' in line:
                idx = int(line[2:])
            if '>' not in line and line != "":
                sequences[idx] = line

    candidate = "_"
    while candidate not in stopAA and keepgoing:
        # Elect candidate
        candidate = "_"
        candidates, scoresB = {}, {}
        for i in range(lenBaitModels):
            if baitModelsStats[i][4] == 0:
                if baitModels[i] != bait:
                    wholebaitmodel = True

            c = "_"
            if 0 <= indices[i] < len(sequences[i]):
                c = sequences[i][indices[i]]

                if c not in "_[]":
                    if c not in candidates:
                        candidates[c], scoresB[c] = 0, 0
                    candidates[c] += 0.01 if c == '-' else 1
                    scoresB[c] += scoreBM(baitModelsStats[i])

        mostpresent, doubt = 0, True
        for k, v in candidates.items():
            doubt = False
            if v > mostpresent:
                candidate, mostpresent = k, v
            elif v == mostpresent:
                if scoresB[k] > scoresB[candidate]:
                    candidate, mostpresent = k, v
                else:
                    doubt = scoresB[k] == scoresB[candidate]

        # - is candidate...
        if candidate == '-':
            # pass
            stopped = True # Mark as stopped but ignore - and keep going
        # No clear candidate
        elif candidate not in mono or mostpresent == 0 or doubt:
            stopped, keepgoing = True, False
            break
        if candidate != '-':
            fusedBait += candidate # clear candidate
        for i in range(lenBaitModels):
            indices[i] += 1

    # STEP5 : clean (remove unused files)
    os.close(fd) # Close file descriptor
    if os.path.exists(IN):
        os.remove(IN)
    if os.path.exists(OUT):
        os.remove(OUT)
    if os.path.exists(IN+".dnd"):
        os.remove(IN+".dnd")

    # Print result if verbose
    if verbose:
        print("Fusion : ", fusedBait)

    # Some stats
    lenbait = len(bait)
    inBM = bait in baitModels
    isequal, numMatch = compare(bait, fusedBait, True)
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
    totalchar += lenbait
    charcount += numMatch
    solvedbaits += 1 if isequal else 0
# ---------------------------------------------

# ----------------- RESULTS -------------------
printResults(solvedbaits, numBait, charcount, totalchar, results, fulltable)
# ---------------------------------------------

# ------------------- PLOTS -------------------
options = " (sans simplification)" if not cansimplify else ""
lengthBaitModels, respBM = tuple(zip(*[t[0] for t in sorted(zip(resultsPerBM.items()))]))
plotResults(options, lengthBaitModels, respBM)
# ---------------------------------------------
