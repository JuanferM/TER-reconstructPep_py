# Fichier csv

Fichier .csv avec quatre colonnes.

Les colonnes 1 et 2 sont des peptides : le bait et le hit associés par le logiciel SpecOMS car ils ont un nombre de pics en commun (SPC) supérieur ou égal au seuil (ici 7).

La colonne 3 est le hitModified renvoyé par le logiciel SpecGlob.

La colonne 4 est le baitModel ; pour le créer le hitModified est traité par un script maison afin que les modifications 1) soient correctement placées 2) interprétables soient transformées en séquences.

Chaque bait a un ou plusieurs hit(s). Par exemple le premier bait du fichier, VCASIAQK, est relié à 6 hits.

Les fichiers contiennent :

- Les 10 000 premiers baits avec plus d'un hit, au total 56 412 PSMs (lien entre hit et bait)
- Les 100 000 premiers baits avec plus d'un hit, au total 703 354 PSMs
- L'intégralité des PSMs renvoyés par SpecOMS (même pour les baits qui n'ont qu'un seul hit), soient 455 404 baits et 5 164 896 PSMs

# Conversion du fichier csv
Le script convert.py permet de convertir ce fichier .csv en une instance dans
un autre format sans les hitModified (que l'on n'utilise pas a priori, voir le
README dans le répertoire instances). De plus, un fichier contenant des
statistiques sur les données (voir README dans le répertoire instances/stats)

Pour calculer toutes les statistiques, utilisez la commande :
```
python3 convert.py data.csv output table.csv clean.csv
```
Où *data.csv*  correspond au fichier .csv décrit ci-dessus;
   *output*    correspond au nom que l'on donne à l'instance de sortie;
   *table.csv* correspond à la table des masses à utiliser pour déterminer si une
        masse est inconnue ou équivalente à une ou plusieurs combinaisons
        d'acides aminés;
   *clean.csv* correspond au fichier .csv indiquant si un bait est "decoy" ou "target"
