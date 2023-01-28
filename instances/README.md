# Collection of instances

*bait10000.txt* is a collection of instances generated from the csv file
*specOMS_test_for_glob_sample_10000_output_baitModel.csv* (found in the
directory `csv`).

The format of _bait10000.txt_ is as follows :
* the first line is the number of baits
* for each bait we have :
    - the name of the bait
    - the number of baitModels associated to the bait **n**
    - the following **n** lines are the baitModels associated to the bait

Here are the first 9 lines of the file that may help you undestand its structure :
```
10000
VCASIAQK 6
V[160.03]ASIA[128.06]K
[146.11]Q[55.97]SIAQK
V[231.06]SIAQK
V[47.04]QA[71.97]IAQK
VC[342.19]QK
VCASI[199.09]K
GGSGATIIMVVQR 5
```

# Script to check mass
checkMass.py takes the path to a collection and check if all baitModels
assigned to a bait have the same mass.
