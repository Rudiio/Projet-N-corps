import numpy as np
import matplotlib.pyplot as plt

#Lecture des données
E_BH = np.loadtxt("E_BH.txt")
E_NO = np.loadtxt("E_NO.txt")
E_N = np.loadtxt("E_N.txt")

###Affichage des données
##Barnes-Hut
plt.figure(figsize=(8,5))
plt.plot(E_BH[:,0],E_BH[:,1],label="Energie cinétique")
plt.plot(E_BH[:,0],-1*E_BH[:,2],label="Energie potentielle")
plt.plot(E_BH[:,0],E_BH[:,3],label="Energie mécanique")
plt.xlabel("temps")
plt.ylabel("Energie")
plt.title("Evolutions des énergies au cours du temps BH")
plt.legend()
plt.show()

##Barnes-Hut
plt.figure(figsize=(8,5))
plt.plot(E_NO[:,0],E_NO[:,1],label="Energie cinétique")
plt.plot(E_NO[:,0],-1*E_NO[:,2],label="Energie potentielle")
plt.plot(E_NO[:,0],E_NO[:,3],label="Energie mécanique")
plt.xlabel("temps")
plt.ylabel("Energie")
plt.title("Evolutions des énergies au cours du temps NO")
plt.legend()
plt.show()

##Barnes-Hut
plt.figure(figsize=(8,5))
plt.plot(E_N[:,0],E_N[:,1],label="Energie cinétique")
plt.plot(E_N[:,0],-1*E_N[:,2],label="Energie potentielle")
plt.plot(E_N[:,0],E_N[:,3],label="Energie mécanique")
plt.xlabel("temps")
plt.ylabel("Energie")
plt.title("Evolutions des énergies au cours du temps N")
plt.legend()
plt.show()