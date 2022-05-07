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
    std :: cout <<"1: Barnes-Hut 2: Naïve efficace 3: Naïve\n";
    std :: cin >> methode_calcul;
    std :: cout << "\n";

    NBodyWnd wndMain(700, "NBody Simulation (Barnes Hut algorithm)");

    // Define simulation size
    //Modifier Init pour qu'elle prenne le nombre de particules à générer
    //implique des modifications d'autres fonctions 

    wndMain.Init(num,methode_calcul);
    wndMain.MainLoop();
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

