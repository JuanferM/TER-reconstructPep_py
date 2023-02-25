import csv
from math import exp, trunc, ceil
from rich import box
from rich.table import Table
from rich.console import Console
import matplotlib.text as mtext
import matplotlib.pyplot as plt

# Used to generate plots
class LegendTitle(object):
    def __init__(self, text_props=None):
        self.text_props = text_props or {}
        super(LegendTitle, self).__init__()

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        title = mtext.Text(x0, y0, r'\underline{' + orig_handle + '}', usetex=True, **self.text_props)
        handlebox.add_artist(title)
        return title

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
def simplifyBM(originalBaitModels, baitModelsStats, massTable, uncertainty,
               NUMS):
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
            baitModels[i] = newBaitModel
    return baitModels

# Make stats table
def make_table(color="", showline=False, showheader=True):
    params = {"padding": (0,1,0,1), "pad_edge": False, "expand": True,
              "style": color, "header_style": "bold " + color,
              "show_edge": False, "show_lines": showline,
              "show_header": showheader, "box": box.ASCII_DOUBLE_HEAD}
    return Table(**params)

# Add stats table's column
def add_column(table, text, justify="center", color=""):
    table.add_column(text, justify=justify, style=color)

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
    return R

# Read stats file
def readStatsFile(fname, minNumBaits=0, maxNumBaits=float('inf')):
    print("Reading stats file...")
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
            for i in range(i+1, i+1+nBaitModel):
                baitModel, SPC, SGscore, baitMass, LS, nMass, GSC, GMC, GUM = rows[i].split()
                baitModels.append(baitModel)
                st = [int(SPC), float(SGscore), float(baitMass)]
                st = st + list(map(int, [LS, nMass, GSC, GMC, GUM]))
                baitModelsStats.append(st)
            i += 1
            if minNumBaits < nBaitModel <= maxNumBaits:
                numBait += 1
                totalInBM += 1 if bait in baitModels else 0
                baits[bait] = (baitModels, baitStats, baitModelsStats)
    print("Done\n")
    return nBait, nBaitOne, numBait, totalInBM, baits

# Read mass table
def readMassTable(fname):
    print("Reading mass table...")
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
    print("Done\n")
    return massTable

# Print file statistics
def printStats(verbose, trace, uncertainty, numBait, nBait, nBaitOne,
               massDispCount, minCount, maxCount, moyCount, totalInBM):
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
    print("# of bait sequence incl. in bait models  : ", totalInBM)

# Print results
def printResults(solvedbaits, numBait, charcount, totalchar, results,
                 fulltable=False):
    print("\nSolved baits\t\t : {} / {} ({:.2f} %)".format(solvedbaits,
                                                     numBait,
                                                     (solvedbaits/numBait)*100))
    print("# of matching characters : {} / {} ({:.2f} %)\n".format(charcount,
                                                     totalchar,
                                                     (charcount/totalchar)*100))

    entries = {"Reconstitution complète" : ["VP (= bait)", "FP (≠ bait)"],
               "Reconstitution incomplète": ["FN (= bait)", "80% du bait",
                                             "50% du bait", "30% du bait",
                                             "Reste (≠ bait)"]}
    colors = ["green", "dark_orange", "blue", "red"]
    table = Table(title="Statistiques linearBaitFusion.py", padding=(0,0,0,0),
                  box=box.ASCII_DOUBLE_HEAD)

    idx, lencol = 0, len(colors)
    symbolIn, symidx = False, -1
    columns = []
    for entry, subentries in entries.items():
        subtables = []
        add_column(table, entry)
        for subentry in subentries:
            if not fulltable and "%" in subentry:
                if not symbolIn:
                    symbolIn, symidx = True, idx
                idx += 1
                continue
            cl = colors[-1] if idx >= lencol else colors[idx]
            subtotal = sum(results[idx*3:3+idx*3])
            subtable = make_table(cl, True)
            subtable2, subtable3 = make_table(), make_table()

            add_column(subtable, subentry, color=cl)
            add_column(subtable2, "bait incl.\ndans les\nbaitmodels")
            add_column(subtable2, "bait non\nincl. dans les\nbaitmodels")
            add_column(subtable3, "baitmodel ≠\ndu bait sans\nmodification")
            add_column(subtable3, "reste")

            if not symbolIn:
                subtable3.add_row(str(results[1+idx*3]), str(results[2+idx*3]))
                subtable2.add_row(str(results[idx*3]), subtable3)
            else:
                subtotal += sum(results[symidx*3:idx*3])
                subtable3.add_row(str(sum(results[1+symidx*3::3])),
                                  str(sum(results[2+symidx*3::3])))
                subtable2.add_row(str(sum(results[symidx*3::3])), subtable3)
            subtable.add_row(subtable2)
            subtable.add_row(str(subtotal))

            subtables.append(subtable)
            idx += 1

        mergetable = make_table(showheader=False)
        mergetable.add_row(*subtables)
        columns.append(mergetable)
    table.add_row(*columns)

    console = Console()
    console.print(table)

# Plot resutls
def plotResults(options, lengthBaitModels, respBM):
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

