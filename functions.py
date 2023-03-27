import csv
from math import exp, trunc, ceil
from rich import box
from rich.table import Table
from rich.console import Console
import matplotlib.text as mtext
import matplotlib.pyplot as plt

class LegendTitle(object):
    """
    Used to generate plots (see StackOverflow answer at https://stackoverflow.com/a/38486135)
    """
    def __init__(self, text_props=None):
        self.text_props = text_props or {}
        super(LegendTitle, self).__init__()

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        title = mtext.Text(x0, y0, r'\underline{' + orig_handle + '}', usetex=True, **self.text_props)
        handlebox.add_artist(title)
        return title

def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    @params:
        number      - Required  : number to truncate (Float)
        decimals    - Optional  : positive number of decimals places (Int)
    @return truncated value of number
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return trunc(number)

    factor = 10.0 ** decimals
    return trunc(number * factor) / factor

def compare(bait, fusedBait):
    """
    Compute longest common contiguous subsequence to determine if the two
    sequences are the same and how much character they have in common (see
    https://stackoverflow.com/a/14044116)
    @params:
        bait        - Required  : sequence of the bait (Str)
        fusedBait   - Required  : sequence obtained with algorithm (Str)
    @return a Tuple (b, c) where b is a Boolean set to True if the two
    sequences are the same (False otherwise) and c is the number of characters
    that the two sequences have in common.
    """
    A, B = bait, fusedBait
    lenA, lenB = len(A), len(B)
    n = (lenA if lenA > lenB else lenB) + 1
    LCCS = [[0]*n for _ in range(n)]

    for _ in range(lenA, n-1):
        A += '-'
    for _ in range(lenB, n-1):
        B += '-'

    r, x, y = -1, -1, -1

    for i in range(1, lenA+1):
        for j in range(1, lenB+1):
            if A[i-1] == B[j-1]:
                LCCS[i][j] = LCCS[i-1][j-1]+1
            else:
                LCCS[i][j] = 0

            if LCCS[i][j] > r:
                x, y = i, j
                r = LCCS[i][j]

    r = 0 if r == -1 else r
    return (A == B), r

def scoreBM(stats):
    """
    Compute score for a baitModel
    @params:
        stats       - Required  : stats on the given baitModel (stats should be
                                  stored in a tuple, see function readStatsFile)
    @return the score of the baitModel
    """
    _, _, _, LS, _, GSC, GMC, GUM = stats
    return exp(LS/(10+10*GSC+100*GMC+1000*GUM))

def getMass(mono, sequence):
    """
    Compute mass of sequence
    @params:
        mono        - Required  : dictionary with single amino acid as key and
                                  the mass of the given amino acid as value (Dict)
        sequence    - Required  : the sequence of amino acids (Str)
    @return the mass of the sequence
    """
    mass = 0.0
    for aa in sequence:
        if aa in mono:
            mass += mono[aa]
    return mass

def simplifyBM(mono, originalBaitModels, baitModelsStats, massTable, tolerance,
               sensitivity):
    """
    Simplify baitModels
    @params:
        mono                - Required  : dictionary with single amino acid as key and
                                          the mass of the given amino acid as
                                          value (Dict)
        originalBaitModels  - Required  : the original baitModels (Str List)
        baitModelStats      - Required  : stats of all the baitModels [see readStatsFile]
                                          (Int Tuple)
        massTable           - Required  : mass table [see readMassTable] (Dict)
        tolerance           - Required  : error margin when looking up a mass
                                          in the mass table (Float)
        sensitivity         - Required  : smallest absolute amount of change to
                                          mass when looking up a mass in a the mass
                                          table [until sensitivity == tolerance] (Float)
    @return simplified baitModels [if the process can't be applied, e.g.
            a baitModel with no unknown mass, the original baitModel is
            returned] (Str List)
    """
    NUMS = "-.0123456789"
    baitModels = originalBaitModels.copy()

    for i in range(len(baitModels)):
        # If there is more than one unknown mass in the BM
        if baitModelsStats[i][7] > 1:
            firstunk, j, k = -1, 0, 0
            firstmass, masses, sequences = False, [], []
            seq, mas = "", ""

            # First read the different subsequences and masses
            for c in baitModels[i]:
                if j == 0 and c == '[':
                    firstmass = True
                j += 1
                if c == '[':
                    if seq != "":
                        sequences.append(seq)
                    seq = ""
                elif c == ']' and mas != "":
                    l, step, n = 0, 0, tolerance/sensitivity
                    mass, ncombi = truncate(float(mas), 2), -1
                    while 0 <= l < 3*n:
                        sign = -1 if l % 2 == 1 else 1
                        step += (1 if j % 2 == 1 else 0)
                        queryMass = truncate(mass + sign*sensitivity*step, 2)
                        l += 1
                        if str(abs(queryMass)) in massTable:
                            l = -abs(queryMass)
                            ncombi = len(massTable[str(abs(queryMass))])

                    # If the unknown mass is small enough
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

            # If a mass is unknown in the table mass, sum that mass with the
            # next (or previous) mass
            offset = 0 if firstmass else 1
            if firstunk == lenM-1:
                masses[firstunk-1] += masses[firstunk]
                masses[firstunk-1] += getMass(mono, sequences[firstunk-1+offset])
                masses[firstunk-1] = truncate(masses[firstunk-1], 2)
                masses[firstunk] = ''
                sequences[firstunk-1+offset] = ''
            else:
                masses[firstunk] += masses[firstunk+1]
                masses[firstunk] += getMass(mono, sequences[firstunk+offset])
                masses[firstunk] = truncate(masses[firstunk], 2)
                masses[firstunk+1] = ''
                sequences[firstunk+offset] = ''

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

            # Replace existing baitModel with the new version
            baitModels[i] = newBaitModel
    return baitModels

def make_table(color="", showline=False, showheader=True):
    """
    Utility function to build stats table (see printResults)
    """
    params = {"padding": (0,1,0,1), "pad_edge": False, "expand": True,
              "style": color, "header_style": "bold " + color,
              "show_edge": False, "show_lines": showline,
              "show_header": showheader, "box": box.ASCII_DOUBLE_HEAD}
    return Table(**params)

def add_column(table, text, justify="center", color=""):
    """
    Utility function to add column to stats table (see printResults)
    """
    table.add_column(text, justify=justify, style=color)

def fillResults(R, inBM, wholebaitmodel):
    """
    Utility function to fill results array according to booleans
    @params:
        R               - Required  : results array (Int List)
        inBM            - Required  : True if the bait is in the baitModels
                                      (Bool)
        wholebaitmodel  - Required  : True if a sequence without any mass that
                                      is different from the bait is in the
                                      baitModels (Bool)
    @return results array (Int List)
    """
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
    return R

def readStatsFile(fname, minNumBaits=0, maxNumBaits=float('inf'),
                  ignoreDuplicateBM=False):
    """
    Read stats file
    @params:
        fname             - Required  : Path to the stats file (Str)
        minNumBaits       - Optional  : Only load baits with at least minNumBaits
                                        baitModels (Int)
        maxNumBaits       - Optional  : Only load baits with at most maxNumBaits
                                        baitModels (Int)
        ignoreDuplicateBM - Optional  : Only load baits with at most maxNumBaits
                                        baitModels (Int)
    @return statistics on file and a dictionary with baits as key and all the
    data about each baits as value (Tuple)
    """
    numBait, totalInBM, baits = 0, 0, {}
    with open(fname, 'r') as file:
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
            offset = 0
            for i in range(i+1, i+1+nBaitModel):
                baitModel, SPC, SGscore, baitMass, LS, nMass, GSC, GMC, GUM = rows[i].split()
                if ignoreDuplicateBM and baitModel in baitModels:
                    offset += 1
                    continue
                baitModels.append(baitModel)
                st = [int(SPC), float(SGscore), float(baitMass)]
                st = st + list(map(int, [LS, nMass, GSC, GMC, GUM]))
                baitModelsStats.append(st)
            nBaitModel -= offset
            i += 1
            if minNumBaits <= nBaitModel <= maxNumBaits:
                numBait += 1
                totalInBM += 1 if bait in baitModels else 0
                baits[bait] = (baitModels, baitStats, baitModelsStats)
    return baits, nBait, nBaitOne, numBait, totalInBM, massDispCount, minCount, maxCount, moyCount

def readMassTable(fname):
    """
    Read mass table csv file
    @params:
        fname           - Required  : Path to the mass table csv file (Str)
    @return dictionary with mass as key and a list of all possible combinations
            of up to eight amino acids which mass is equal to the key as value
            (Dict)
    """
    massTable = {}
    with open(fname, 'r') as file:
        reader = csv.reader(file, delimiter=',')
        numrow = 0
        for row in reader:
            if numrow != 0:
                massTable[row[0]] = []
                for combi in row[1:len(row)]:
                    if combi != '':
                        massTable[row[0]].append(combi)
            numrow += 1
    return massTable

def printStats(verbose, trace, tolerance, sensitivity, numBait, nBait, nBaitOne,
               massDispCount, minCount, maxCount, moyCount, totalInBM):
    """
    Print file statistics
    """
    print("Verbose                                  : ", verbose)
    print("Trace                                    : ", trace, " Da")
    print("Tolerance                                : ", tolerance, " Da")
    print("Sensitivity                              : ", sensitivity, " Da")
    print("# of baits                               : ", numBait)
    print("# of baits in stats file                 : ", nBait)
    print("# of baits with at least two bait models : ", nBaitOne)
    print("# of baits with mass dispersion > 1.0 Da : ", massDispCount)
    print("Min # of baitModels in stats file        : ", minCount)
    print("Max # of baitModels in stats file        : ", maxCount)
    print("Avg # of baitModels in stats file        : ", truncate(moyCount, 2))
    print("# of bait sequence incl. in bait models  : ", totalInBM)

def getResultsTables(results, resultsPercent, fulltable=False):
    """
    Print table according to the results [see printResults]
    @params:
        results         - Required  : Results obtained with method
                                      [look for "# Some stats" in main file]
                                      (Int List)
        resultsPercent  - Required  : Results obtained with method (percent)
                                      [look for "# Some stats" in main file]
                                      (Float List)
        fulltable       - Optional  : Print table with all details or not (Bool)
    """
    mt, bt, fnd = "More than ", "Between ", " found"
    entries = {"Full reconstruction" : [("= bait", "= bait"),
                                        (mt+"21 AA"+fnd, mt+"80% AA"+fnd),
                                        (bt+"14 and 20 AA"+fnd, bt+"80% and 50% AA"+fnd),
                                        (bt+"7 and 13 AA"+fnd, bt+"50% and 30% AA"+fnd),
                                        ("≠ bait", "≠ bait")],
               "Partial reconstruction" : [(mt+"21 AA"+fnd, mt+"80% AA"+fnd),
                                           (bt+"14 and 20 AA"+fnd, bt+"80% and 50% AA"+fnd),
                                           (bt+"7 and 13 AA"+fnd, bt+"50% and 30% AA"+fnd),
                                           ("≠ bait", "≠ bait")]}
    table = Table(title="baitFusion detailed results", padding=(0,0,0,0),
                  box=box.ASCII_DOUBLE_HEAD)
    tablePercent = Table(title="baitFusion detailed results (percent)",
                  padding=(0,0,0,0), box=box.ASCII_DOUBLE_HEAD)

    columns, columnsP = [], []
    previousAAinEntry = False
    idx, contract, conidx = 0, False, -1
    for entry, subentries in entries.items():
        subtables, subtablesP = [], []
        add_column(table, entry)
        add_column(tablePercent, entry)
        for subentry, subentryP in subentries:
            AAinEntry = "AA" in subentry or "AA" in subentryP
            if not fulltable and AAinEntry:
                if not contract:
                    contract, conidx = True, idx
                idx += 1
                continue
            cl = "green" if (idx == 0) else "red"
            subtotal = sum(results[idx*3:3+idx*3])
            subtable = make_table(cl, True)
            subtable2, subtable3 = make_table(), make_table()
            # percent subtables
            subtotalP = sum(resultsPercent[idx*3:3+idx*3])
            subtableP = make_table(cl, True)
            subtable2P, subtable3P = make_table(), make_table()

            if fulltable and not AAinEntry and previousAAinEntry:
                add_column(subtable, "Others ("+subentry+")", color=cl)
                add_column(subtableP, "Others ("+subentryP+")", color=cl)
            else:
                add_column(subtable, subentry, color=cl)
                add_column(subtableP, subentryP, color=cl)
            add_column(subtable2, "bait incl. in \nbaitmodels")
            add_column(subtable2, "bait not incl.\nin baitmodels")
            add_column(subtable3, "baitmodel ≠ from\nbait without\nmodification")
            add_column(subtable3, "others")
            # create percent subtables columns
            add_column(subtable2P, "bait incl. in \nbaitmodels")
            add_column(subtable2P, "bait not incl.\nin baitmodels")
            add_column(subtable3P, "baitmodel ≠ from\nbait without\nmodification")
            add_column(subtable3P, "others")

            if not contract:
                subtable3.add_row(str(results[1+idx*3]), str(results[2+idx*3]))
                subtable2.add_row(str(results[idx*3]), subtable3)
                # fill percent subtables
                subtable3P.add_row(str(resultsPercent[1+idx*3]),
                                   str(resultsPercent[2+idx*3]))
                subtable2P.add_row(str(resultsPercent[idx*3]), subtable3P)
            else:
                subtotal += sum(results[conidx*3:idx*3])
                subtable3.add_row(str(sum(results[1+conidx*3::3])),
                                  str(sum(results[2+conidx*3::3])))
                subtable2.add_row(str(sum(results[conidx*3::3])), subtable3)
                # accumulate data if not fulltable (percent version)
                subtotalP += sum(resultsPercent[conidx*3:idx*3])
                subtable3P.add_row(str(sum(resultsPercent[1+conidx*3::3])),
                                   str(sum(resultsPercent[2+conidx*3::3])))
                subtable2P.add_row(str(sum(resultsPercent[conidx*3::3])),
                                   subtable3P)
                if contract and AAinEntry:
                    contract, conidx = False, -1
            subtable.add_row(subtable2)
            subtable.add_row(str(subtotal))
            # Add to subtable (percent version)
            subtableP.add_row(subtable2P)
            subtableP.add_row(str(subtotalP))

            subtables.append(subtable)
            subtablesP.append(subtableP)
            idx += 1
            previousAAinEntry = AAinEntry

        mergetable = make_table(showheader=False)
        mergetable.add_row(*subtables)
        columns.append(mergetable)

        mergetableP = make_table(showheader=False)
        mergetableP.add_row(*subtablesP)
        columnsP.append(mergetableP)
    table.add_row(*columns)
    tablePercent.add_row(*columnsP)

    return table, tablePercent

def printResults(solvedbaits, numBait, results, resultsPercent, fulltable=False):
    """
    Print results of the fusion of the baitModels
    @params:
        solvedbaits     - Required  : Number of solved baits (Int)
        numBait         - Required  : Total number of baits (Int)
        results         - Required  : Results obtained with method
                                      [look for "# Some stats" in main file]
                                      (Int List)
        resultsPercent  - Required  : Results obtained with method (percent)
                                      [look for "# Some stats" in main file]
                                      (Float List)
        fulltable       - Optional  : Print table with all details or not (Bool)
    """
    f = lambda a, b : a/b if b != 0 else 0
    fullRecEqBait = sum(results[0:3])
    fullRecNeqBait = sum(results[3:15])
    partialRecBait = sum(results[15:])
    total = fullRecEqBait + fullRecNeqBait + partialRecBait
    totalRec = total - partialRecBait
    print()
    print("# of baits found                         : {} / {} ({:.2f} %)".format(
          solvedbaits, numBait, f(solvedbaits, numBait)*100))
    print("# of fully reconstructed sequences (FRS) : {} / {} ({:.2f} %)".format(
          totalRec, total, f(totalRec, total)*100))
    print("# of baits found w/ respect to FRS       : {} / {} ({:.2f} %)".format(
          fullRecEqBait, totalRec, f(fullRecEqBait, totalRec)*100))

    table, tablePercent = getResultsTables(results, resultsPercent, fulltable)

    print()
    console = Console()
    if console.is_terminal:
        console.print(table, "\n")
        console.print(tablePercent)
    else:
        console = Console(width=400)
        console.print(table, "\n")
        console.print(tablePercent)

def writeResults(fname, csvheader, csvdata):
    """
    Store method output to csv file
    @params:
        fname               - Required  : csv file name
        csvheader           - Required  : header of csv file
        csvdata             - Required  : data to store in csv file
    """
    # open file in write mode
    with open(fname, 'w', encoding='UTF8', newline='') as f:
        # create csv writer
        writer = csv.writer(f)
        # write the header
        writer.writerow(csvheader)
        # write rows
        writer.writerows(csvdata)

def plotResults(options, lengthBaitModels, respBM):
    """
    Plot results
    @params:
        options             - Required  : Options of the algorithm [see PLOTS
                                          section in main file] (Str)
        lengthBaitModels    - Required  : All numbers of baitModels (Int List)
        respBM              - Required  : Results for baits per number of
                                          baitModels (Int List)
    """
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
    lab = ["baits retrouvés", "baits non retrouvés       "]
    leg2 = ["nombre de baits", ""] + lab + [""] + lab
    pEqNi = list(map(lambda L: L[4], respBM))
    pDiNi = list(map(lambda L: L[5], respBM))
    pEqIn = list(map(lambda L: L[6], respBM))
    pDiIn = list(map(lambda L: L[7], respBM))
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
    plt.savefig("number_plot2{}.png".format(options))

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
    plt.savefig("proportion_plot2{}.png".format(options))

def printProgressBar(solvedbaits, iteration, total, prefix = '', suffix = '', decimals = 2, length = 30, printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar (modified version of https://stackoverflow.com/a/34325723)
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    fillsol, fillit = '█', '▒'
    percent = ("{0:." + str(decimals) + "f}").format(100 * (solvedbaits / float(total)))
    solvedLength = int(length * solvedbaits // total)
    filledLength = int(length * iteration // total) - solvedLength
    bar = fillsol * solvedLength + fillit * filledLength + '-' * (length - solvedLength - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()
