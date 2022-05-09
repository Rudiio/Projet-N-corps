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
  try
  {
    //Choix du nombre de Particules
    int num=1000;
    num = std::stoi(argv[1]);

    //Choix de la méthode 
    int methode_calcul=1;
    methode_calcul=std::stoi(argv[2]);
    
    //Nombre maximale d'itération
    int Nmax = 100;

    NBodyWnd wndMain(700, "NBody Simulation (Barnes Hut algorithm)");

    // Define simulation size
    //Modifier Init pour qu'elle prenne le nombre de particules à générer
    //implique des modifications d'autres fonctions 

    wndMain.Init(num,methode_calcul);
    wndMain.MainLoop(num,methode_calcul,Nmax);
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

