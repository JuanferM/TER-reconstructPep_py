# Sequence alignment algorithms

Sequence alignment is an optimal arrangement of characters that
maximizes the number of matches and minimizes the number of mismatches
and gaps. Sequence alignment algorithms are often based on a scoring
scheme where there's a reward for matches and penalty for mismatches.
These algorithms are used for DNA/RNA and amino acids sequence alignment.

There are two types of alignment :
* Local alignment align regions having highest similarities. It aligns
substring of target with substring of query and it is suitable for more
divergent sequences.
* Global alignment tries to align entire sequence and align all letters
from query to target (suitable for closely related sequences).

* Local alignment algorithms
    - Smith-Waterman
* Global alignment algorithms
    - Needleman-Wunsch
