import numpy as np
import matplotlib.pyplot as plt

#Lecture des données
BH_par = np.loadtxt("donnees_temps_BH.txt")
NO_par = np.loadtxt("donnees_temps_NO.txt")
N_par = np.loadtxt("donnees_temps_N.txt")

BH_seq = np.loadtxt("seq_donnees_temps_BH.txt")
NO_seq = np.loadtxt("seq_donnees_temps_NO.txt")
N_seq = np.loadtxt("seq_donnees_temps_N.txt")

par_tree = np.loadtxt("BH_Tree_Construction_time.txt")
seq_tree = np.loadtxt("seq_BH_Tree_Construction_time.txt")

taille = NO_par.shape[0]
N= np.linspace(0,22000,50)
###Affichage des données

##Comparaison entre méthodes

##Arbre
plt.figure(figsize=(8,5))
plt.plot(par_tree[:,0],par_tree[:,1],"+-",label="parallèle")
plt.plot(seq_tree[:,0],seq_tree[:,1],"+-",label="sequentiel")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en secondes")
plt.title("Durée moyenne de construction de l'arbre")
plt.legend()
plt.show()

###Séquentiel
plt.figure(figsize=(8,5))
plt.plot(BH_seq[0:taille,0],BH_seq[0:taille,1],"+-",label="Barnes Hut")
plt.plot(NO_seq[:,0],NO_seq[:,1],"+-",label="Naïve Optimisée")
plt.plot(N_seq[:,0],N_seq[:,1],"+-",label="Naïve")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en secondes")
plt.title("Durée moyenne de calcul par rapport au nombre de particules séquentiel")
plt.legend()
plt.show()

###Parallèle
plt.figure(figsize=(8,5))
plt.plot(BH_par[0:taille,0],BH_par[0:taille,1],"+-",label="Barnes Hut")
plt.plot(NO_par[:,0],NO_par[:,1],"+-",label="Naïve Optimisée")
plt.plot(N_par[:,0],N_par[:,1],"+-",label="Naïve")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en secondes")
plt.title("Durée moyenne de calcul par rapport au nombre de particules parallèle")
plt.legend()
plt.show()

##Comparaison parallèle/sequentiel
#Barnes-Hut
plt.figure(figsize=(8,5))
plt.plot(BH_seq[:,0],BH_seq[:,1],"+-",label="sequentiel")
plt.plot(BH_par[:,0],BH_par[:,1],"+-",label="parallèle")
# plt.plot(N,N*np.log(N),"+-",label="$Nlog(N)$")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en secondes")
plt.title("Durée moyenne de calcul par rapport au nombre de particules BH")
plt.legend()
plt.show()

#Naïve optimisée
plt.figure(figsize=(8,5))
plt.plot(NO_seq[:,0],NO_seq[:,1],"+-",label="sequentiel")
plt.plot(NO_par[:,0],NO_par[:,1],"+-",label="parallèle")
# plt.plot(N,N*N,"+-",label="$N^2$")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en secondes")
plt.title("Durée moyenne de calcul par rapport au nombre de particules NO")
plt.legend()
plt.show()

#Naïve
plt.figure(figsize=(8,5))
plt.plot(N_seq[:,0],N_seq[:,1],"+-",label="sequentiel")
plt.plot(N_par[:,0],N_par[:,1],"+-",label="parallèle")
# plt.plot(N,N*N,"+-",label="$N^2$")
plt.xlabel("Nombre de particules")
plt.ylabel("Durée en secondes")
plt.title("Durée moyenne de calcul par rapport au nombre de particules N")
plt.legend()
plt.show()

