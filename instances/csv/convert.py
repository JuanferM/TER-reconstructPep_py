import csv
import sys
import math
from ctypes import ArgumentError

if len(sys.argv) < 3 and len(sys.argv) > 4:
    raise FileNotFoundError('File not found. Please provide the csv filename\
 and the output instance name as arguments.')

uncertainty = 0.01
AA = "ARNDCEQGHILKMFPSTWYVCU"
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
baits, hitModified, massTable = {}, {}, {}
baitStats, baitModelsStats = {}, {}
infilename, outfilename, massfilename = sys.argv[1], sys.argv[2], ""
if len(sys.argv) == 4:
    massfilename = sys.argv[3]

def truncate(number, decimals=0):
    """
    Returns a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor

# ------------- READING CSV FILE --------------
print("Reading csv file...")
with open(infilename, 'r') as file:
    reader = csv.reader(file, delimiter=';')
    counter = 0
    for row in reader:
        if counter != 0:
            if row[0] not in baits:
                baits[row[0]] = []
            baits[row[0]].append(row[3])
            if row[0] not in hitModified:
                hitModified[row[0]] = []
            hitModified[row[0]].append(row[2])
        counter += 1
print("Done\n")
# ---------------------------------------------

# ------------- WRITE INSTANCE ----------------
print("Generating instance...")
with open(outfilename, 'w') as file:
    file.write(str(len(baits.keys()))+"\n")
    for k, v in baits.items():
        file.write(k+" "+str(len(baits[k]))+"\n")
        for baitModels in baits[k]:
            file.write(baitModels+"\n")
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

# ------------- COMPUTE STATS -----------------
print("Computing stats...")
minCount, maxCount, moyCount = sys.maxsize, 0, 0 # General info
for bait, L in baits.items():
    sz = len(L)
    minCount = min(minCount, sz)
    maxCount = max(maxCount, sz)
    moyCount += sz
    baitModelsStats[bait], baitModelsMass = [], []
    # Compute longest stretch, mass stats, number of offset/mass
    for baitModel in L:
        baitMass, currentMass = 0, ""
        corMass, ambMass, unkMass = 0, 0, 0
        longstretch, stretch, numOffset = -1, 0, 0
        for c in baitModel:
            if c in AA:
                stretch += 1
                baitMass += mono[c]
            else:
                if c == '[':
                    numOffset += 1
                elif c in NUMS:
                    currentMass += c
                elif c == ']' and massfilename != "":
                    i, mass, ncombi = 0, float(currentMass), -1
                    baitMass += mass
                    while i < 3 and i != -1:
                        mass = mass + (-1)*(i%2)*uncertainty*i
                        i += 1
                        if str(abs(mass)) in massTable:
                            i, ncombi = -1, len(massTable[str(abs(mass))])
                    currentMass = ""
                    if ncombi == -1:
                        unkMass += 1
                    elif ncombi == 1:
                        corMass += 1
                    else:
                        ambMass += 1
                longstretch = max(longstretch, stretch)
                stretch = 0
        baitModelsMass.append(baitMass)
        longstretch = max(longstretch, stretch)
        baitModelsStats[bait].append((baitModel, str(truncate(baitMass, 2)), str(longstretch),
                                      str(numOffset), str(corMass), str(ambMass), str(unkMass)))
    # Compute mass standard deviation and mean
    meanMass = sum(baitModelsMass)/sz
    quadmean = 0
    for i in range(sz):
        quadmean += (baitModelsMass[i] - meanMass) ** 2
    sdMass = math.sqrt(quadmean/sz)
    baitStats[bait] = (truncate(meanMass, 2), truncate(sdMass, 2))

moyCount /= len(baits)
# TODO Stats sur les hitModifieds
print("Done\n")
# ---------------------------------------------

print("Writing stats file...")
with open("stats_"+outfilename, 'w') as file:
    # write number of baits; min, max and mean number of baitModels
    file.write(str(len(baits.keys())) + " ")
    file.write(str(minCount) + " " + str(maxCount) + " ")
    file.write(str(moyCount) + "\n")

    # for each bait...
    for k, v in baitModelsStats.items():
        sz = str(len(baitModelsStats[k]))
        meanMass, sdMass = baitStats[k]
        # write bait and number of baitModels
        file.write(k+" "+sz+" "+str(meanMass)+" "+str(sdMass)+"\n")
        # for each baitModel
        for stats in baitModelsStats[k]:
            # write baitModel, bait mass, longest stretch, number of mass, ...
            file.write(" ".join(stats)+"\n")
print("Done\n")
