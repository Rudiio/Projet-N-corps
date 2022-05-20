/*
 * File:   main.cpp
 * Author: user
 *
 * Created on 5. Juli 2009, 22:15
 */

#include <cstdlib>
#include <iostream>
#include <string>

//------------------------------------------------------------------------------
#include "NBodyWnd.h"
 

//------------------------------------------------------------------------------
/*
 *
 */
int main(int argc, char** argv)
{
  /* Permet de lancer une simulation en passant en argument le nombre de particules et
  la méthode à utiliser */
  
  try
  {
    //Choix du nombre de Particules
    int num=1000;
    num = std::stoi(argv[1]);

    //Choix de la méthode 
    int methode_calcul=1;
    methode_calcul=std::stoi(argv[2]);

    NBodyWnd wndMain(700, "NBody Simulation (Barnes Hut algorithm)");

    // Define simulation size
    //Le deuxième paramètre de Mainloop est le nombre d'itérations à effectuer
    //S'il vaut -1 le programme tourne tant que l'utilisateur n'a pas fermé le programme 
    wndMain.Init(num,methode_calcul);
    wndMain.MainLoop(num,methode_calcul,100);
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

