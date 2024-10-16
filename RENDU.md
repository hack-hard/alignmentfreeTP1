# Eleve
- tristan hacquard
# Rapport de TD
Dans ce TD nous avons implementés une méthode pour calculer la similarité de jaccard.
Nous avons subdivisés les séquences en kmers de longeuer 20 et nous en avons sélectionné 10 000 par la méthode de minhash. Cependant, l'optimisation du tas en 16 bits n'as pas marché pour des raisons obscures et est données à titre indicative.

# Table de calculs
|                 | GCA_000069965.1 | GCA_000005845.2 | GCA_030271835.1 | GCA_000008865.2 | GCA_000013265.1 |
| --------------- | --------------- | --------------- | --------------- | --------------- | --------------- |
| GCA_000069965.1 |                 | 0.00257         | 0.0311          | 0.00232         | 0.00244         |
| GCA_000005845.2 |                 |                 | 0.00258         | 0.436           | 0.341           |
| GCA_030271835.1 |                 |                 |                 | 0.00232         | 0.00244         |
| GCA_000008865.2 |                 |                 |                 |                 | 0.307           |
| GCA_000013265.1 |                 |                 |                 |                 |                 |

Concernant l'alignement avec les Prokaryotes, l'algorithme continue de tourner et n'a pas donné de résultats.

|                  | GCA_029289425.3 | GCF_000001635.27 |
| ---------------- | --------------- | ---------------- |
| GCA_029289425.3  |                 | 1.0              |
| GCF_000001635.27 |                 |                  |

