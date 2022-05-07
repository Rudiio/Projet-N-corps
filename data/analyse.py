import numpy as np
import matplotlib.pyplot as plt

#Lecture des données
BH_par = np.loadtxt("donnees_temps_BH.txt")
NE_par = np.loadtxt("donnees_temps_NE.txt")
N_par = np.loadtxt("donnees_temps_N.txt")

BH_seq = np.loadtxt("seq_donnees_temps_BH.txt")
NE_seq = np.loadtxt("seq_donnees_temps_NE.txt")
N_seq = np.loadtxt("seq_donnees_temps_N.txt")


###Affichage des données

##Comparaison entre méthodes

###Séquetiel
plt.plot(BH_seq[:,0],BH_seq[:,1],label="BH")
plt.plot(NE_seq[:,0],NE_seq[:,1],label="NE")
plt.plot(N_seq[:,0],N_seq[:,1],label="N")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en s")
plt.title("Durée moyenne de calcul par rapport au nombre de particules")
plt.legend()
plt.show()

###Parallèle
plt.plot(BH_par[:,0],BH_par[:,1],label="BH")
plt.plot(NE_par[:,0],NE_par[:,1],label="NE")
plt.plot(N_par[:,0],N_par[:,1],label="N")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en s")
plt.title("Durée moyenne de calcul par rapport au nombre de particules")
plt.legend()
plt.show()

##Comparaison parallèle/sequentiel
#Barnes-Hut
plt.plot(BH_seq[:,0],BH_seq[:,1],label="seq")
plt.plot(BH_par[:,0],BH_par[:,1],label="parallèle")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en s")
plt.title("Durée moyenne de calcul par rapport au nombre de particules BH")
plt.legend()
plt.show()

#Naïve optimisée
plt.plot(NE_seq[:,0],NE_seq[:,1],label="seq")
plt.plot(NE_par[:,0],NE_par[:,1],label="parallèle")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en s")
plt.title("Durée moyenne de calcul par rapport au nombre de particules NE")
plt.legend()
plt.show()

#Naïve
plt.plot(N_seq[:,0],N_seq[:,1],label="seq")
plt.plot(N_par[:,0],N_par[:,1],label="parallèle")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en s")
plt.title("Durée moyenne de calcul par rapport au nombre de particules N")
plt.legend()
plt.show()

