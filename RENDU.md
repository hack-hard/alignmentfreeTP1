# Eleve
- tristan hacquard
# Rapport de TD
Dans ce TD nous avons implementés une méthode pour calculer une similarité entre séquence sans aligements.
Nous avons subdivisés les séquences en kmers et comparé la proportion de kmers qui étaient dans les 2 séquences.

# Table de calculs
|     | A   | B     | C      | D      | E      |
| --- | --- | ----- | ------ | ------ | ------ |
| A   |     | 0.001 | 0.028  | 0.002  | 0.002  |
| B   |     |       | 0.0018 | 0.435  | 0.334  |
| C   |     |       |        | 0.0016 | 0.0019 |
| D   |     |       |        |        | 0.299  |
| E   |     |       |        |        |        |
avec 
A =GCA_000069965.1
B =GCA_000005845.2 
C=GCA_030271835.1 
D =GCA_000008865.2 
E = GCA_000013265.1 

On observe que l'on a une homologie très forte entre B,D et E.


