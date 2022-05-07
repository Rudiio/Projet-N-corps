import matplotlib.pyplot as plt
import numpy as np


temps = []

for line in open(r'/home/elyas/Bureau/github/projet_N_corps_git/git/Projet-N-corps/src/donnees_temps.txt'):
    lines = [ i for i in line.split()]
    temps.append(lines[0])

print(temps)


taille = len(temps)
x = np.linspace(0,1,taille)

plt.figure(figsize=(10,10))
plt.tick_params(labelleft = False )
plt.plot(x, temps)
plt.title("Graphique reprensentant le temps passé dans les calculs de notre model à N-corps")
plt.ylabel('temps des calculs')
# plt.show()
plt.savefig( "/home/elyas/Bureau/github/projet_N_corps_git/git/Projet-N-corps/visuels/Graphique_durée_des_calculs.png")