Fichier .csv avec quatre colonnes.

Les colonnes 1 et 2 sont des peptides : le bait et le hit associés par le logiciel SpecOMS car ils ont un nombre de pics en commun (SPC) supérieur ou égal au seuil (ici 7).

La colonne 3 est le hitModified renvoyé par le logiciel SpecGlob.

La colonne 4 est le baitModel ; pour le créer le hitModified est traité par un script maison afin que les modifications 1) soient correctement placées 2) interprétables soient transformées en séquences.

Chaque bait a un ou plusieurs hit(s). Par exemple le premier bait du fichier, VCASIAQK, est relié à 6 hits.

Ce fichier contient un échantillon de 10 000 baits, au total 56 412 PSMs (lien entre hit et bait).
