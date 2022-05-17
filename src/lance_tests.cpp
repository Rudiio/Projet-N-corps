#include <iostream>
#include <cstdio>

using namespace std;


int main(int argc,char **argv)
{
    int i=1;
    // Première série de tests
    for(int n=100;n<2000;n+=100){ 
        
        char sbh[30];
        char sno[30];
        char sn[30];

        //Barnes Hut
        i=1;
        sprintf(sbh,"./main_mesure %d %d",n,1);
        i=system(sbh);

        if(i!=0)
        {
            exit(EXIT_FAILURE);
            cout << "erreur\n";
        }

        //naïve optimisée
        i=1;
        sprintf(sno,"./main_mesure %d %d",n,2);
        i=system(sno);

        if(i!=0)
        {
            exit(EXIT_FAILURE);
            cout << "erreur\n";
        }

        //naïve
        i=1;
        sprintf(sn,"./main_mesure %d %d",n,3);
        i=system(sn);

        if(i!=0)
        {
            exit(EXIT_FAILURE);
            
        }
    }

    //Deuxième série de tests
    for(int n=22000;n<=70000;n+=2000){ 
        
        char sbh[30];
        char sno[30];
        char sn[30];

        //Barnes Hut
        i=1;
        sprintf(sbh,"./main_mesure %d %d",n,1);
        i=system(sbh);

        if(i!=0)
        {
            exit(EXIT_FAILURE);
            cout << "erreur\n";
        }

        //naïve optimisée
        i=1;
        sprintf(sno,"./main_mesure %d %d",n,2);
        i=system(sno);

        if(i!=0)
        {
            exit(EXIT_FAILURE);
            cout << "erreur\n";
        }

        //naïve
        i=1;
        sprintf(sn,"./main_mesure %d %d",n,3);
        i=system(sn);

        if(i!=0)
        {
            exit(EXIT_FAILURE);
            
        }
    }

    return 0;
}