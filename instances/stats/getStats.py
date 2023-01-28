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
nBait, minCount, maxCount, moyCount = 0, 0, 0, 0
labCol = ["mass (Da)", "LS", "# G", "# GSC", "# GMC", "# GUM"]

# ------------- READING STATS FILE ------------
print("Reading stats file...")
with open(infilename, 'r') as file:
    i, rows = 1, file.read().split('\n')
    nBait, minCount, maxCount, moyCount = rows[0].split()
    moyCount = float(moyCount)
    nBait, minCount, maxCount = map(int, (nBait, minCount, maxCount))
    while i < len(rows)-1:
        bait, nBaitModel, meanMass, sdMass = rows[i].split()
        nBaitModel = int(nBaitModel)
        baitStats = (float(meanMass), float(sdMass))
        baitModels, baitModelsStats = [], []
        for i in range(i+1, i+1+nBaitModel):
            baitModel, baitMass, stretch, nMass, corMass, ambMass, unkMass = rows[i].split()
            if len(baitModel) > baitPrintLength:
                baitModel = baitModel[:baitPrintLength]+"..."
            baitModels.append(baitModel)
            st = [float(baitMass)]
            st = st + list(map(int, [stretch, nMass, corMass, ambMass, unkMass]))
            baitModelsStats.append(st)
        i += 1
        baits[bait] = (baitModels, baitStats, baitModelsStats)
print("Done\n")
# ---------------------------------------------

# ------------ PRINTING STATS FILE ------------
print("LS  = LongestStretch")
print("G   = Gaps")
print("GSC = Gaps with mass corresponding to a Single Combination")
print("GMC = Gaps with mass corresponding to Multiple Combinations")
print("GUM = Gaps with Unknown Mass\n")

print("Number of baits     : ", nBait)
print("Min # of baitModels : ", minCount)
print("Max # of baitModels : ", maxCount)
print("Avg # of baitModels : ", moyCount, '\n')

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