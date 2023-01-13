import csv
import sys
from ctypes import ArgumentError

baits = {}

if len(sys.argv) != 3:
    raise FileNotFoundError('File not found. Please provide the csv filename\
 and the output instance name as arguments.')

print("Reading csv file...")
with open(sys.argv[1], 'r') as file:
    reader = csv.reader(file, delimiter=';')
    counter = 0
    for row in reader:
        if counter != 0:
            if row[0] not in baits:
                baits[row[0]] = []
            baits[row[0]].append(row[3])
        counter += 1

print("Done\nGenerating instance...")
with open(sys.argv[2], 'w') as file:
    file.write(str(len(baits.keys()))+"\n")
    for k, v in baits.items():
        file.write(k+" "+str(len(baits[k]))+"\n")
        for baitModels in baits[k]:
            file.write(baitModels+"\n")
print("Done")
