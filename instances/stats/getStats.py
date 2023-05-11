import sys
import pandas as pd
from ctypes import ArgumentError

if len(sys.argv) != 2:
    raise FileNotFoundError('File not found. Please provide the stats file name')
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

baits = {}
infilename = sys.argv[1]
baitPrintLength = 30
nPSM = 0
labCol = ["SPC", "SGS", "mass (Da)", "LS", "# G", "# GSC", "# GMC", "# GUM"]
nBait, nBaitOne, massDispCount, minCount, maxCount, moyCount = 0, 0, 0, 0, 0, 0

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
            if len(baitModel) > baitPrintLength:
                baitModel = baitModel[:baitPrintLength]+"..."
            baitModels.append(baitModel)
            st = [int(SPC), float(SGscore), float(baitMass)]
            st = st + list(map(int, [LS, nMass, GSC, GMC, GUM]))
            baitModelsStats.append(st)
        i += 1
        nPSM += len(baitModels)
        baits[bait] = (baitModels, baitStats, baitModelsStats)
print("Done\n")
# ---------------------------------------------

# ------------ PRINTING STATS FILE ------------
print("PSM = Peptide-Spectrum Matches")
print("SPC = Shared Peak Count")
print("SGS = SpecGlob Score")
print("LS  = Longest Stretch")
print("G   = Gaps")
print("GSC = Gaps with mass corresponding to a Single Combination")
print("GMC = Gaps with mass corresponding to Multiple Combinations")
print("GUM = Gaps with Unknown Mass\n")

print("Number of PSM                                    : ", nPSM)
print("Number of baits                                  : ", nBait)
print("Number of baits with at least one bait model     : ", nBaitOne)
print("Number of baits with mass dispersion > 1.0 Da    : ", massDispCount)
print("Min # of baitModels                              : ", minCount)
print("Max # of baitModels                              : ", maxCount)
print("Avg # of baitModels                              : ", moyCount, '\n')

bait = ""
while bait != "q":
    bait = input("Insert bait ('all' to get the stats for all baits, 'q' to exit) : ")
    if bait == "all":
        for k in baits.keys():
            print("Stats for ", k)
            baitModels, baitStats, baitModelsStats = baits[k]
            meanMass, sdMass = baitStats
            df = pd.DataFrame(baitModelsStats, columns = labCol, index = baitModels)
            print(df, '\n', '-'*(baitPrintLength+3))
            print("Mean mass (Da)               : ", "{:.2f}".format(meanMass))
            print("Mass standard deviation (Da) : ", "{:.2f}".format(sdMass), '\n')
    elif bait != "q":
        if bait in baits:
            baitModels, baitStats, baitModelsStats = baits[bait]
            meanMass, sdMass = baitStats
            df = pd.DataFrame(baitModelsStats, columns = labCol, index = baitModels)
            print(df, '\n', '-'*(baitPrintLength+3))
            print("Mean mass (Da)               : ", "{:.2f}".format(meanMass))
            print("Mass standard deviation (Da) : ", "{:.2f}".format(sdMass), '\n')
        else:
            print("Sorry, bait not found.\n")
# ---------------------------------------------
