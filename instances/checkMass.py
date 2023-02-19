import sys
import math
from ctypes import ArgumentError

if len(sys.argv) != 2:
    raise FileNotFoundError('File not found. Please provide the stats file name')

infilename = sys.argv[1]
uncertainty = 1
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
AA = "ARNDCEQGHILKMFPSTWYVU"
NUMS = "0123456789."

print("uncertainty : ", uncertainty, " Da")

# ---------- CHECK BAITMODELS MASS ------------
print("\nChecking baitModels mass...")
good = True
with open(infilename, 'r') as file:
    rows = file.read().split('\n')
    i, nBait = 1, int(rows[0])
    while i < nBait-1:
        bait, nBaitModels = rows[i].split()
        nBaitModels = int(nBaitModels)
        pmass = -1
        for i in range(i+1, i+1+nBaitModels):
            mass, currentMass = 0, ""
            for c in rows[i]:
                if c == '[':
                    currentMass = ""
                elif c in AA:
                    mass += mono[c]
                elif c in NUMS:
                    currentMass += c
                elif c == ']':
                    mass += float(currentMass)

            if pmass != -1:
                good = abs(mass - pmass) <= uncertainty
            if not good:
                print("{:.2f}".format(abs(mass-pmass), 2), " > ", uncertainty)
                break
            pmass = mass
        i += 1
        if not good:
            raise ValueError("baitModels mass mismatch for bait " + bait)
print("Done")
# ---------------------------------------------
