import os
import sys
import tempfile
from functions import *
import subprocess as sbp
from math import trunc, ceil

# ------------ ARGUMENTS CHECK ----------------
if len(sys.argv) != 3:
    raise FileNotFoundError('File not found. Please provide the stats filename\
 and the path to the mass table.')
infilename, massfilename = sys.argv[1], sys.argv[2]
# ---------------------------------------------

# ---------------- PARAMETERS -----------------
# Display settings
verbose = False
fulltable = True

# Error margins and baitModels weights
trace = 1
tolerance = 0.02
sensitivity = 0.01

# Method options
simplification = True
ignoreDuplicateBM = False

# Settings about solvers
useMUSCLE = False
clustalopt = "-QUICKTREE -MATRIX=ID -GAPOPEN=5 -GAPEXT=1 -NOHGAP \
-NOWEIGHTS -CLUSTERING=UPGMA"
muscleopt = ""

# Settings about baitModels
onlythisbait = ""
minNumBaits = 2
maxNumBaits = float('inf')
# ---------------------------------------------

# ---------- VARIABLES DEFINITION -------------
musclecmd = "./muscle3.8"
clustalcmd = "./clustalw2 -ALIGN -QUIET -OUTPUT=FASTA"
results = [0] * 27
resultsPercent = [0] * 27
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
csvdata = []
csvheader = ["Bait", "Output", "Bait reconstructed", "Length of bait",
             "Length of longest common continuous sequence (LCCS)",
             "Number of amino acids not in LCCS"]
# ---------------------------------------------

# ------------- READ MASS TABLE ---------------
print("Reading mass table... ", end='')
massTable = readMassTable(massfilename)
print("Done")
# ---------------------------------------------

# ------------- READING STATS FILE ------------
print("Reading stats file... ", end='')
baits, totalBait, totalBaitWithOneBM, numBait, totalBaitInBM, numBaitWithMassDispersion,\
minBMcount, maxBMcount, meanBMcount = readStatsFile(infilename, minNumBaits,
                                                    maxNumBaits, ignoreDuplicateBM)
print("Done\n")
# ---------------------------------------------

# --------------- PRINT STATS -----------------
printStats(numBait, totalBait, totalBaitWithOneBM, numBaitWithMassDispersion, minBMcount,
           maxBMcount, meanBMcount, totalBaitInBM)
# ---------------------------------------------

# ---------- PRINT PARAMETERS ------------
printParameters(verbose, trace, tolerance, sensitivity, simplification, None,
                ignoreDuplicateBM, minNumBaits, maxNumBaits, useMUSCLE)
# ---------------------------------------------

# ------------ FUSION BAIT MODELS -------------
iteration = 1 # Iteration counter to print progress bar
for bait, data in baits.items():
    csvrow = [bait]
    wholebaitmodel = False
    keepgoing, stopped = True, False
    fusedBait, candidate = "", "_"
    baitModels, baitStats, baitModelsStats = data
    lenBaitModels = len(baitModels)
    originalBaitModels = []
    convertedBaitModels = []
    indices = [0 for _ in range(lenBaitModels)]

    # If onlythisbait is defined then we only run the algorithm on onlythisbait
    # if it is found
    if onlythisbait != "" and bait != onlythisbait:
        continue
    if verbose:
        print("\nBait   : ", bait)

    # STEP0 : simplify baitmodels
    if simplification:
        originalBaitModels = baitModels.copy()
        baitModels = simplifyBM(mono, baitModels, baitModelsStats, massTable,
                                tolerance, sensitivity)
        # For some cases simplifying will convert a baitModel into one big mass
        # (e.g. [14.15]G[140.3] ==> [211.47]). In such cases, if we try to
        # replace the mass by '-' we'll have a sequence with no amino acids and
        # straight '-'s (e.g. ------). Feeding such a sequence with no amino acid
        # to Clustal or MUSCLE will ALWAYS THROW AN ERROR. So to avoid that
        # problem we revert the baitModel to its unsimplified form if such
        # a case is encountered. NB : this is done only for method 2 (alignBaitFusion.py)
        for i in range(lenBaitModels):
            noaminoacid = True
            for c in baitModels[i]:
                if c not in "[]"+NUMS:
                    noaminoacid = False
                    break
            if noaminoacid:
                baitModels[i] = originalBaitModels[i]

    # STEP1 : convert baitmodels
    # Each mass is replaced by a certain number of '-' (see details below)
    for i in range(lenBaitModels):
        converted, currentMass = "", ""
        for c in baitModels[i]:
            if c == ']' and currentMass != "":
                mass = float(currentMass)

                # If there is enough excess mass
                if abs(mass) >= trace:
                    j, step, n = 0, 0, tolerance/sensitivity
                    mass, ncombi = truncate(mass, 2), -1
                    while 0 <= j < 3*n:
                        sign = -1 if j % 2 == 1 else 1
                        step += (1 if j % 2 == 1 else 0)
                        queryMass = truncate(mass + sign*sensitivity*step, 2)
                        j += 1
                        if str(abs(queryMass)) in massTable:
                            j = -abs(queryMass)
                            ncombi = len(massTable[str(abs(queryMass))])

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
    scoreBM = []
    for i in range(lenBaitModels):
        scoreBM.append(computeScoreBM(baitModelsStats[i]))
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

                # If c is an amino acid then we can use it as candidate
                if c not in "_[]":
                    if c not in candidates:
                        candidates[c], scoresB[c] = 0, 0
                    candidates[c] += 0.15 if c == '-' else 1
                    scoresB[c] += scoreBM[i]

        # Determine who is the elected candidate
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

        # - is candidate mark as stopped but ignore - and keep going
        if candidate == '-':
            stopped = True
        # If there is no clear elected candidate
        elif candidate not in mono or mostpresent == 0 or doubt:
            stopped, keepgoing = True, False
            break

        # clear elected candidate
        if candidate != '-':
            fusedBait += candidate
        # advance cursors
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

    # Print result if verbose else print progress bar
    if verbose:
        print("Fusion : ", fusedBait)
    else:
        # Print only if not redirecting to file
        if sys.stdout.isatty():
            if iteration == 1:
                print()
            if onlythisbait != "":
                numBait = 1
            printProgressBar(solvedbaits, iteration, numBait, prefix = 'Progress:', suffix
                 = 'Solved')

    # Some stats
    lenbait = len(bait)
    inBM = bait in baitModels
    isequal, numMatch = compare(bait, fusedBait)
    if lenBaitModels not in resultsPerBM:
        resultsPerBM[lenBaitModels] = [0] * 8
    if '-' not in fusedBait:
        if isequal:
            resultsPerBM[lenBaitModels][0] += 1
            results[0:3] = fillResults(results[0:3], inBM, wholebaitmodel)
            resultsPercent[0:3] = fillResults(resultsPercent[0:3], inBM, wholebaitmodel)
        else:
            resultsPerBM[lenBaitModels][1] += 1
            if numMatch >= 21:
                results[3:6] = fillResults(results[3:6], inBM, wholebaitmodel)
            elif numMatch >= 14:
                results[6:9] = fillResults(results[6:9], inBM, wholebaitmodel)
            elif numMatch >= 7:
                results[9:12] = fillResults(results[9:12], inBM, wholebaitmodel)
            else:
                results[12:15] = fillResults(results[12:15], inBM, wholebaitmodel)

            ratio = numMatch/lenbait
            if ratio >= 0.8:
                resultsPercent[3:6] = fillResults(resultsPercent[3:6], inBM, wholebaitmodel)
            elif ratio >= 0.5:
                resultsPercent[6:9] = fillResults(resultsPercent[6:9], inBM, wholebaitmodel)
            elif ratio >= 0.3:
                resultsPercent[9:12] = fillResults(resultsPercent[9:12], inBM, wholebaitmodel)
            else:
                resultsPercent[12:15] = fillResults(resultsPercent[12:15], inBM, wholebaitmodel)
    else:
        resultsPerBM[lenBaitModels][3] += 1
        if numMatch >= 21:
            results[15:18] = fillResults(results[15:18], inBM, wholebaitmodel)
        elif numMatch >= 14:
            results[18:21] = fillResults(results[18:21], inBM, wholebaitmodel)
        elif numMatch >= 7:
            results[21:24] = fillResults(results[21:24], inBM, wholebaitmodel)
        else:
            results[24:27] = fillResults(results[24:27], inBM, wholebaitmodel)

        ratio = numMatch/lenbait
        if ratio >= 0.8:
            resultsPercent[15:18] = fillResults(resultsPercent[15:18], inBM, wholebaitmodel)
        elif ratio >= 0.5:
            resultsPercent[18:21] = fillResults(resultsPercent[18:21], inBM, wholebaitmodel)
        elif ratio >= 0.3:
            resultsPercent[21:24] = fillResults(resultsPercent[21:24], inBM, wholebaitmodel)
        else:
            resultsPercent[24:27] = fillResults(resultsPercent[24:27], inBM, wholebaitmodel)
    if not inBM:
        if isequal:
            resultsPerBM[lenBaitModels][4] += 1
        else:
            resultsPerBM[lenBaitModels][5] += 1
    else:
        if isequal:
            resultsPerBM[lenBaitModels][6] += 1
        else:
            resultsPerBM[lenBaitModels][7] += 1

    # Data for output.csv file
    csvrow.append(fusedBait)
    csvrow.append(str(isequal))
    csvrow.append(str(len(bait)))
    csvrow.append(str(numMatch))
    offset = 1 if '-' in fusedBait else 0
    csvrow.append(str(len(fusedBait)-numMatch-offset))
    csvdata.append(csvrow)

    # Totals
    solvedbaits += 1 if isequal else 0
    iteration += 1
# ---------------------------------------------

# ----------------- RESULTS -------------------
if solvedbaits != 0:
    printResults(solvedbaits, numBait, results, resultsPercent, fulltable)
    version = "MUSCLE" if useMUSCLE else "Clustal"
    writeResults("output_align{}.csv".format(version), csvheader, csvdata)
else:
    print("\nNO SOLUTION!")
# ---------------------------------------------

# ------------------- PLOTS -------------------
options = " (sans simplification)" if not simplification else ""

if solvedbaits != 0:
    lengthBaitModels, respBM = tuple(zip(*[t[0] for t in sorted(zip(resultsPerBM.items()))]))
    plotResults(options, lengthBaitModels, respBM)
# ---------------------------------------------
