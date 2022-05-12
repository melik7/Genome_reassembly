# Genome_reassembly


Cet outil permet de représenter et d’indexer les k-mers d’un fichier de séquençage. A partird’un fichier contenant l’ensemble des séquences (reads) issus du séquençage d’un génome et suite à une sélection des k-mers solides, l’outil représente via une séquence employant la notion de SPSS l’ensemble de k-mers retenus. La qualité de la SPSS est évaluée en comparant le nombre de faux positifs et négatifs entre la séquence SPSS et un génome de référence. Cette comparaison s’effectue de deux manières : une première comparaison est effectuée en utilisant la fonction in de python et la seconde en employant la méthode d’indexation de la Burrows-Wheeler (FM-index).
Pour un bon fonctionnement de l'outil, tous les arguments suivants doivent être précisés.

Arguments :

● -i [str] : fichier contenant les reads au format fasta
● -g [str] : fichier contenant le génome au format fasta
● -t [int] : Le seuil de solidité, tous les k-mers ayant un nombre d'occurrence supérieur à ce seuil sont conservés.
● -k [int] : la taille des k-mers

Exemple de commande :

python projet_algoo.py -i ecoli_sample_reads_1.fasta -g ecoli_genome_150k.fa -t 3 -k 31

En sortie, le programme affiche les informations suivantes :
- Temps pour compter et filtrer les k-mers canoniques solides
- Nombre de caractère total de la SPSS
- Nombre de séquences distinctes dans la SPSS
- Temps de construction de la SPSS
- Nombre de faux négatifs (FN) et le temps de calcul pour les compter sans ou avec
FM-index*
- Nombre de faux positifs (FP) et le temps de calcul pour les compter sans ou avec
FM-index
