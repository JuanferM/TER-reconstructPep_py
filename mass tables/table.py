import math

code = ["", "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L",
        "K", "M", "F", "P", "S", "T", "W", "Y", "V", "O", "U"]

mono = {"A" : 71.037113805,
        "R" : 156.10111105,
        "N" : 114.04292747,
        "D" : 115.026943065,
        "C" : 103.009184505,
        "E" : 129.042593135,
        "Q" : 128.05857754,
        "G" : 57.021463735,
        "H" : 137.058911875,
        "I" : 113.084064015,
        "L" : 113.084064015,
        "K" : 128.09496305,
        "M" : 131.040484645,
        "F" : 147.068413945,
        "P" : 97.052763875,
        "S" : 87.032028435,
        "T" : 101.047678505,
        "W" : 186.07931298,
        "Y" : 163.063328575,
        "V" : 99.068413945,
        "O" : 237.147726925,
        "U" : 150.953633405
        }

comb = {}
mass = {}

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

def tocsv(dico, file):
    with open(file, 'w') as f:
        f.write(str(len(dico))+"\n")
        for k, v in dico.items():
            f.write(str(k)+",")
            for i in v:
                f.write(i+",")
            f.write("\n")

print("Computing...")
counter = 0
for i in code:
    for j in code:
        for k in code:
            for l in code:
                for m in code:
                    for n in code:
                        for o in code:
                            for p in code:
                                counter += 1
                                word = "".join(sorted(i+j+k+l+m+n+o+p))
                                if word != "" and word not in comb:
                                    comb[word] = True
                                    s = 0
                                    for c in word:
                                        if c != '':
                                            s += mono[c]
                                    s = truncate(s, 2)
                                    if s in mass:
                                        mass[s].append(word)
                                    else:
                                        mass[s] = [word]

print("Done.\nWriting to csv...")
tocsv(mass, "table.csv")
print("Done.")
