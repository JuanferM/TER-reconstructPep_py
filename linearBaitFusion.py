import csv
import sys
from math import exp, trunc
from ctypes import ArgumentError
import matplotlib.text as mtext
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    raise FileNotFoundError('File not found. Please provide the stats filename\
 and the path to the mass table.')
infilename, massfilename = sys.argv[1], sys.argv[2]

# ---------------- PARAMETERS -----------------
verbose = False
onlythisbait = ""
minNumBaits = 1
maxNumBaits = float('inf') # included
uncertainty, trace = 0.01, 1.0
valid, probation, invalid = 4, 1, 0
canreverse, reconstructFromBoth = True, True
cansimplify, simplifyBothWays = False, False
canreturnBMasis = False
# ---------------------------------------------

results = [0] * 21
nBait, nBaitOne, numBait = 0, 0, 0
charcount, totalchar, totalInBM = 0, 0, 0
solvedbaits, unsolvedbaits, nobaits = 0, 0, 0
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
def compare(bait, fusedBait, bothWays=False):
    same, i, lasteq = 0, 0, -1
    lenBait, lenFusedBait = len(bait), len(fusedBait)
    while i < lenBait and i < lenFusedBait:
        if bait[i] == fusedBait[i]:
            same += 1
            lasteq = i
        i = i+1

    j, k = lenBait-1, lenFusedBait-1
    while bothWays and j > lasteq and k >= 0:
        same += 1 if bait[j] == fusedBait[k] else 0
        j, k = j-1, k-1

    # Return True if the sequences are the same (false otherwise)
    # Also return the number of equal characters
    return (same == lenBait and same == lenFusedBait), same

# Compute score for a baitModel
def scoreBM(stats):
    _, _, _, LS, _, GSC, GMC, GUM = stats
    return exp(LS/(10+10*GSC+100*GMC+1000*GUM))

# Simplify baitModels
def simplifyBM(originalBaitModels, baitModelsStats):
    global massTable, uncertainty, NUMS
    baitModels = originalBaitModels.copy()
    for i in range(len(baitModels)):
        # If there is more than one unknown mass in the BM
        if baitModelsStats[i][7] > 1:
            firstunk, j, k = -1, 0, 0
            firstmass, masses, sequences = False, [], []
            seq, mas = "", ""

            for c in baitModels[i]:
                if j == 0 and c == '[':
                    firstmass = True
                j += 1
                if c == '[':
                    if seq != "":
                        sequences.append(seq)
                    seq = ""
                elif c == ']' and mas != "":
                    l, mass, ncombi = 0, truncate(float(mas), 2), -1
                    while 0 <= l < 3:
                        mass = truncate(mass + (-1)*(l%2)*uncertainty*l, 2)
                        l += 1
                        if str(abs(mass)) in massTable:
                            l, ncombi = -abs(mass), len(massTable[str(abs(mass))])

                    if firstunk == -1 and ncombi == -1:
                        firstunk = k
                    masses.append(float(mas))
                    mas, k = "", k+1
                elif c in NUMS:
                    mas += c
                else:
                    seq += c
            if seq != "":
                sequences.append(seq)
            if mas != "":
                masses.append(float(mas))
            lenM, lenS = len(masses), len(sequences)

            # Sum masses
            if firstunk == lenM-1:
                masses[firstunk-1] += masses[firstunk]
                masses[firstunk-1] += truncate(masses[firstunk-1], 2)
                masses[firstunk] = ''
            else:
                masses[firstunk] += masses[firstunk+1]
                masses[firstunk] = truncate(masses[firstunk], 2)
                masses[firstunk+1] = ''

            # Reverse corresponding sequence
            if firstunk+1 >= lenS:
                sequences[firstunk-1] = sequences[firstunk-1][::-1]
            elif firstmass:
                sequences[firstunk] = sequences[firstunk][::-1]
            else:
                sequences[firstunk+1] = sequences[firstunk+1][::-1]

            # Reconstruct new baitModel
            newBaitModel, m = "", min(lenM, lenS)
            for j in range(max(lenM, lenS)):
                if j < m:
                    if firstmass:
                        newBaitModel += '['+str(masses[j])+']' if masses[j] != '' else ''
                        newBaitModel += sequences[j]
                    else:
                        newBaitModel += sequences[j]
                        newBaitModel += '['+str(masses[j])+']' if masses[j] != '' else ''
                else:
                    if lenM > lenS:
                        newBaitModel += '['+str(masses[j])+']' if masses[j] != '' else ''
                    else:
                        newBaitModel += sequences[j]

            # Replace existing
            # print(baitModels[i], " => ", newBaitModel)
            baitModels[i] = newBaitModel
    return baitModels

# Fill results array according to booleans
def fillResults(R, inBM, wholebaitmodel):
    # bait sequence is amongst the baitmodels
    if inBM:
        R[0] += 1
    # bait sequence is not amongst the baitmodels
    # However, there is a sequence without mass != bait amongst the baitmodels
    elif wholebaitmodel:
        R[1] += 1
    # none of the above
    else:
        R[2] += 1

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
print("Avg # of baitModels                      : ", truncate(moyCount, 2))

for bait, data in baits.items():
    wholebaitmodel = False
    startchar, stopchar = '[', ']'
    keepgoing, stopped, reverse = True, False, False
    fusedBait, befRevBait, candidate, = "", "", "_"
    baitModels, baitStats, baitModelsStats = data
    lenBaitModels = len(baitModels)
    masses = [0.0 for _ in range(lenBaitModels)]
    indices = [0 for _ in range(lenBaitModels)]
    validation = [valid for _ in range(lenBaitModels)]
    originals = []
    if onlythisbait != "" and bait != onlythisbait:
        continue
    if verbose:
        print("\nBait   : ", bait)

    # STEP0 : simplify baitmodels
    if cansimplify:
        originals = baitModels
        baitModels = simplifyBM(baitModels, baitModelsStats)

    # Find sequence
    while (candidate not in stopAA or reverse) and keepgoing:
        # STEP1 : Elect candidate
        candidate = "_"
        candidates, scoresA, scoresB = {}, {}, {}
        for i in range(lenBaitModels):
            if baitModelsStats[i][5] == 0:
                if baitModels[i] != bait:
                    wholebaitmodel = True
                if canreturnBMasis:
                    fusedBait, keepgoing = baitModels[i], False
                    break

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
                if cansimplify and not simplifyBothWays:
                    baitModels = originals
                continue
            else:
                stopped, keepgoing = True, False
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
    bothWays = reverse and stopped and reconstructFromBoth
    if reverse:
        if bothWays:
            fusedBait = befRevBait + fusedBait[::-1]
        else:
            cond = len(fusedBait) > len(befRevBait)
            fusedBait = fusedBait[::-1] if cond else befRevBait

    # Print result if verbose
    if verbose:
        print("Fusion : ", fusedBait)

    # Some stats
    lenbait = len(bait)
    inBM = bait in baitModels
    isequal, numMatch = compare(bait, fusedBait, bothWays)
    if lenBaitModels not in resultsPerBM:
        resultsPerBM[lenBaitModels] = [0] * 11
    if not stopped:
        if isequal:
            resultsPerBM[lenBaitModels][0] += 1
            resultsPerBM[lenBaitModels][4] += 1
            fillResults(results[0:3], inBM, wholebaitmodel)
        else:
            ratio = numMatch/lenbait
            resultsPerBM[lenBaitModels][1] += 1
            resultsPerBM[lenBaitModels][5] += 1 if ratio >= 0.8 else 0
            resultsPerBM[lenBaitModels][6] += 1 if ratio < 0.8 else 0
            fillResults(results[3:6], inBM, wholebaitmodel)
    else:
        if isequal:
            resultsPerBM[lenBaitModels][2] += 1
            resultsPerBM[lenBaitModels][4] += 1
            fillResults(results[6:9], inBM, wholebaitmodel)
        else:
            ratio = numMatch/lenbait
            resultsPerBM[lenBaitModels][3] += 1
            resultsPerBM[lenBaitModels][5] += 1 if ratio >= 0.8 else 0
            resultsPerBM[lenBaitModels][6] += 1 if ratio < 0.8 else 0
            if ratio >= 0.8:
                fillResults(results[9:12], inBM, wholebaitmodel)
            elif ratio >= 0.5:
                fillResults(results[12:15], inBM, wholebaitmodel)
            elif ratio >= 0.3:
                fillResults(results[15:18], inBM, wholebaitmodel)
            else:
                fillResults(results[18:21], inBM, wholebaitmodel)
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
    totalInBM += 1 if inBM else 0
    solvedbaits += 1 if isequal else 0
# ---------------------------------------------

# ----------------- RESULTS -------------------
if verbose:
    print("\n")
print("# of bait sequence incl. in bait models  : ", totalInBM)
print("\nSolved baits\t\t : {} / {} ({:.2f} %)".format(solvedbaits,
                                                     numBait,
                                                     (solvedbaits/numBait)*100))
print("# of matching characters : {} / {} ({:.2f} %)".format(charcount,
                                                             totalchar,
                                                             (charcount/totalchar)*100))
# ---------------------------------------------

# ------------------- PLOTS -------------------
class LegendTitle(object):
    def __init__(self, text_props=None):
        self.text_props = text_props or {}
        super(LegendTitle, self).__init__()

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        title = mtext.Text(x0, y0, r'\underline{' + orig_handle + '}', usetex=True, **self.text_props)
        handlebox.add_artist(title)
        return title

options = " (" if canreverse or cansimplify else ""
if canreturnBMasis:
    options += "baitmodel as is"
options += " + " if canreturnBMasis and (canreverse or cansimplify) else ""
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

lengthBaitModels, respBM = tuple(zip(*[t[0] for t in sorted(zip(resultsPerBM.items()))]))

# First plot
lab = ["baits retrouvés", "baits non retrouvés"]
leg2 = ["nombre de baits", ""] + lab + [""] + lab
pTotal = list(map(lambda L: sum(L[0:4]), respBM))
pVP = list(map(lambda L: L[0], respBM))
pFP = list(map(lambda L: L[1], respBM))
pFN = list(map(lambda L: L[2], respBM))
pVN = list(map(lambda L: L[3], respBM))
sum2 = [x + y for x, y in zip(pVP, pFN)]

plt.clf()
kpat = plt.bar(lengthBaitModels, pTotal, width=1, color='k')
gpat = plt.bar(lengthBaitModels, pVP, width=1, color='tab:green')
bpat = plt.bar(lengthBaitModels, pFN, width=1, color='tab:blue', bottom=pVP)
opat = plt.bar(lengthBaitModels, pFP, width=1, color='tab:orange', bottom=sum2)
rpat = plt.bar(lengthBaitModels, pVN, width=1, color='tab:red',
        bottom=[x + y for x, y in zip(sum2, pFP)])
leg1 = [kpat, "Reconstitution complète", gpat, opat, "Reconstitution \
incomplète", bpat, rpat]
plt.legend(leg1, leg2, loc=1, handler_map={str: LegendTitle({'fontsize' : 12})})
plt.title("Nombre de baits retrouvés en fonction du nombre de baitModels")
plt.ylabel("Nombre de baits")
plt.xlabel("Nombre de baitModels")
plt.savefig("number_plot1{}.png".format(options))

# Proportion version
pVP = list(map(lambda t: t[0]/t[1], zip(pVP, pTotal)))
pFP = list(map(lambda t: t[0]/t[1], zip(pFP, pTotal)))
pFN = list(map(lambda t: t[0]/t[1], zip(pFN, pTotal)))
pVN = list(map(lambda t: t[0]/t[1], zip(pVN, pTotal)))
sum2 = [x + y for x, y in zip(pVP, pFN)]

plt.clf()
gpat = plt.bar(lengthBaitModels, pVP, width=1, color='tab:green')
bpat = plt.bar(lengthBaitModels, pFN, width=1, color='tab:blue', bottom=pVP)
opat = plt.bar(lengthBaitModels, pFP, width=1, color='tab:orange', bottom=sum2)
rpat = plt.bar(lengthBaitModels, pVN, width=1, color='tab:red',
        bottom=[x + y for x, y in zip(sum2, pFP)])
leg1 = ["Reconstitution complète", gpat, opat, "Reconstitution \
incomplète", bpat, rpat]
plt.legend(leg1, leg2[1:], loc=4, handler_map={str: LegendTitle({'fontsize' : 12})})
plt.title("Proportion de baits retrouvés en fonction du nombre de baitModels")
plt.ylabel("Proportion")
plt.xlabel("Nombre de baitModels")
plt.savefig("proportion_plot1{}.png".format(options))


# Second plot
labels = ["nombre de baits", "baits retrouvés",
          "baits non retrouvés mais\nsimilaire à 80% ou plus", "bait non retrouvés"]
# pTotal = list(map(lambda L: sum(L[4:7]), respBM))
pEqual = list(map(lambda L: L[4], respBM))
pClose = list(map(lambda L: L[5], respBM))
pDiff  = list(map(lambda L: L[6], respBM))

plt.clf()
plt.bar(lengthBaitModels, pTotal, width=1, color='k')
plt.bar(lengthBaitModels, pEqual, width=1, color='tab:green')
plt.bar(lengthBaitModels, pClose, width=1, color='tab:orange', bottom=pEqual)
plt.bar(lengthBaitModels, pDiff, width=1, color='tab:red',
        bottom=[x + y for x, y in zip(pEqual, pClose)])
plt.legend(labels, loc=1)
plt.ylabel("Nombre de baits")
plt.xlabel("Nombre de baitModels")
plt.title("Nombre de baits retrouvés en fonction du nombre de baitModels")
plt.savefig("number_plot2{}.png".format(options))

# Proportion version
pEqual = list(map(lambda t: t[0]/t[1], zip(pEqual, pTotal)))
pClose = list(map(lambda t: t[0]/t[1], zip(pClose, pTotal)))
pDiff = list(map(lambda t: t[0]/t[1], zip(pDiff, pTotal)))

plt.clf()
plt.bar(lengthBaitModels, pEqual, width=1, color='tab:green')
plt.bar(lengthBaitModels, pClose, width=1, color='tab:orange', bottom=pEqual)
plt.bar(lengthBaitModels, pDiff, width=1, color='tab:red',
        bottom=[x + y for x, y in zip(pEqual, pClose)])
plt.legend(labels[1:], loc=4)
plt.ylabel("Proportion")
plt.xlabel("Nombre de baitModels")
plt.title("Proportion de baits retrouvés en fonction du nombre de baitModels")
plt.savefig("proportion_plot2{}.png".format(options))


# Third plot
lab = ["baits retrouvés", "baits non retrouvés       "]
leg2 = ["nombre de baits", ""] + lab + [""] + lab
# pTotal = list(map(lambda L: sum(L[0:4]), respBM))
pEqNi = list(map(lambda L: L[7], respBM))
pDiNi = list(map(lambda L: L[8], respBM))
pEqIn = list(map(lambda L: L[9], respBM))
pDiIn = list(map(lambda L: L[10], respBM))
sum2 = [x + y for x, y in zip(pEqNi, pDiNi)]

plt.clf()
kpat = plt.bar(lengthBaitModels, pTotal, width=1, color='k')
gpat = plt.bar(lengthBaitModels, pEqNi, width=1, color='tab:green')
opat = plt.bar(lengthBaitModels, pDiNi, width=1, color='tab:orange', bottom=pEqNi)
bpat = plt.bar(lengthBaitModels, pEqIn, width=1, color='tab:blue', bottom=sum2)
rpat = plt.bar(lengthBaitModels, pDiIn, width=1, color='tab:red',
        bottom=[x + y for x, y in zip(sum2, pEqIn)])
leg1 = [kpat, "bait incl. dans baitmodels", gpat, opat, "bait non incl. dans \
baitmodels", bpat, rpat]
plt.legend(leg1, leg2, loc=1, handler_map={str: LegendTitle({'fontsize' : 12})})
plt.ylabel("Nombre de baits")
plt.xlabel("Nombre de baitModels")
plt.title("Nombre de baits retrouvés en fonction du nombre de baitModels")
plt.savefig("number_plot3{}.png".format(options))

# Proportion version
pEqNi = list(map(lambda t: t[0]/t[1], zip(pEqNi, pTotal)))
pDiNi = list(map(lambda t: t[0]/t[1], zip(pDiNi, pTotal)))
pEqIn = list(map(lambda t: t[0]/t[1], zip(pEqIn, pTotal)))
pDiIn = list(map(lambda t: t[0]/t[1], zip(pDiIn, pTotal)))
sum2 = [x + y for x, y in zip(pEqNi, pDiNi)]

plt.clf()
gpat = plt.bar(lengthBaitModels, pEqNi, width=1, color='tab:green')
opat = plt.bar(lengthBaitModels, pDiNi, width=1, color='tab:orange', bottom=pEqNi)
bpat = plt.bar(lengthBaitModels, pEqIn, width=1, color='tab:blue', bottom=sum2)
rpat = plt.bar(lengthBaitModels, pDiIn, width=1, color='tab:red',
        bottom=[x + y for x, y in zip(sum2, pEqIn)])
leg1 = ["bait incl. dans baitmodels", gpat, opat, "bait non incl. dans \
baitmodels", bpat, rpat]
plt.legend(leg1, leg2[1:], loc=4, handler_map={str: LegendTitle({'fontsize' : 12})})
plt.ylabel("Proportion")
plt.xlabel("Nombre de baitModels")
plt.title("Proportion de baits retrouvés en fonction du nombre de baitModels")
plt.savefig("proportion_plot3{}.png".format(options))

# ---------------------------------------------
