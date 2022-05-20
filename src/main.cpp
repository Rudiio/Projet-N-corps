/*
 * File:   main.cpp
 * Author: user
 *
 * Created on 5. Juli 2009, 22:15
 */

#include <cstdlib>
#include <iostream>

//------------------------------------------------------------------------------
#include "NBodyWnd.h"
 

//------------------------------------------------------------------------------
/*
 *
 */
int main(int argc, char** argv)
{
  try
  {
    //Choix du nombre de Particules
    int num=1000;
    std :: cout << "Veuillez choisir le nombre de particules\n";
    std ::cin >> num;

    //Choix de la méthode 
    int methode_calcul=1;
    std :: cout << "\nVeuillez choisir une méthode de calcul\n";
    std :: cout <<"1: Barnes-Hut | 2: Naïve optimisée | 3: Naïve\n";
    std :: cin >> methode_calcul;
    std :: cout << "\n";

    //Choix de l'initialisation
    int mode_init=1;
    std :: cout << "\nVeuillez choisir une méthode de calcul\n";
    std :: cout <<"0: Grosse galaxie | 1: Galaxie-atome | 2: 2 galaxies | 3:galaxie sphérique \n";
    std :: cin >> mode_init;
    std :: cout << "\n";

    NBodyWnd wndMain(700, "NBody Simulation (Barnes Hut algorithm)");

    //Le deuxième paramètre de Mainloop est le nombre d'itérations à effectuer
    //S'il vaut -1 le programme tourne tant que l'utilisateur n'a pas fermé le programme 
    wndMain.Init(num,methode_calcul,mode_init);
    wndMain.MainLoop(num,methode_calcul,-1);
  }
  catch(std::exception & exc)
  {
    std::cout << "Fatal error: " << exc.what() << std::endl;
  }
  catch(...)
  {
    std::cout << "Fatal error: unexpected exception" << std::endl;
  }
  return (EXIT_SUCCESS);
}

