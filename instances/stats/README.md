# Statistics

*stats_bait10000.txt* contains some statistics about the data in the csv
file (check instances/csv repository).

The format of *stats_bait10000.txt* is as follows :
* the first line contains the number of baits, the minimum number of
  baitModels, the maximum number of baitModels and the mean number of
  baitModels
* for each bait we have :
    - the name of the bait
    - the number of baitModels associated to the bait **n**
    - the mean of all baitModels masses
    - the standard deviation of all baitModels masses
    - the following **n** lines are the baitModels associated to the bait, the
      shared Peak Count of the baitModel, the specGlob score of the baitModel, the
      mass of the baitModel, the longest stretch in that baitModel, the number
      of gaps in that baitModel, the number of gaps in which the mass
      corresponds to a single combination, the number of gaps in which the mass
      corresponds to multiple combinations and the number of gaps in which the
      mass is unknown

The getStats.py script takes the path to the stats file and allows the user to
retrieve the statistics for a given bait (or all baits).
