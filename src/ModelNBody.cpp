/*
 * File:   ModelMagneticPendulum.cpp
 * Author: user
 *
 * Created on 29. April 2009, 23:08
 */

#include "ModelNBody.h"

//--- Standard includes --------------------------------------------------

#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <omp.h>
#include <sys/time.h>
#include <time.h>
#include <fstream>
#include <string>

using namespace std;

struct timeval start;
struct timeval fin;

//------------------------------------------------------------------------
/** \brief Calcul des durées en secondes */
float time_diff(struct timeval *start, struct timeval *end)
{
    return (end->tv_sec - start->tv_sec) + 1e-6*(end->tv_usec - start->tv_usec);
}

//------------------------------------------------------------------------

/** \brief ModelNbody class constructor*/
ModelNBody::ModelNBody(int num,int methode_calcul,int mode_init) 

  :IModel("N-Body simulation (2D)")
  ,m_pInitial(NULL) //Tableau contenant les positions et vitesses initiales
  ,m_pAux(NULL)     //Tableau contenant les masses de chaque particule
  ,m_root(BHTreeNode(Vec2D(), Vec2D()))
  ,m_min()
  ,m_max()
  ,m_center()
  ,m_camDir()
  ,m_camPos()
  ,m_roi(1)                                                                           //Radial Orbit Instability
  ,m_timeStep(50)
  ,mass_sun(1.988435e30)                                                              //masse soleil M
  ,pc_in_m(3.08567758129e16)                                                          //taux de conversion de parsec en mètre L
  ,gamma_si(6.67428e-11)                                                              //Constante gravitationnelle M-1 L3 T-2
  ,gamma_1(gamma_si/(pc_in_m*pc_in_m*pc_in_m)*mass_sun*(365.25*86400)*(365.25*86400)) //dimension : L3 -> paramètre gravitationnel standard
  ,time_1(sqrt( (pc_in_m*pc_in_m*pc_in_m)/(gamma_si*mass_sun) ) / (365.25*86400))     //dimension : T1/2
  ,m_num(num)
  ,methode_calcul(methode_calcul)
  ,m_bVerbose(false)
  ,mesure_temps(0)
  ,mesure_construction_BH(0)
{
  BHTreeNode::s_gamma = gamma_1;

  gettimeofday(&start, NULL);

  //------------------------------------------------------------------------
  //Initialisation des particules
  InitCollision(num, mode_init);

  gettimeofday(&fin, NULL);
  std::cout << "  Le temps passé dans la fonction InitCollision est : " <<time_diff(&start, &fin) << "sec \n";
}

//------------------------------------------------------------------------
/** \brief Destructor de la classe ModelNBody*/

ModelNBody::~ModelNBody()
{
  delete m_pInitial;
  delete m_pAux;
}
//------------------------------------------------------------------------
float ModelNBody::Getmesuretemps()
{
  return mesure_temps;
}

//------------------------------------------------------------------------
float ModelNBody::Getmesureconstruction()
{
  return mesure_construction_BH;
}

//------------------------------------------------------------------------
void ModelNBody::SetROI(double roi)
{
  m_roi = roi;
}

//------------------------------------------------------------------------
double ModelNBody::GetSuggestedTimeStep() const
{
  return m_timeStep;
}

//------------------------------------------------------------------------
double ModelNBody::GetROI() const
{
  return m_roi;
}

//------------------------------------------------------------------------
Vec3D ModelNBody::GetCenterOfMass() const
{
  const Vec2D &cm2d = m_root.GetCenterOfMass();
  return Vec3D(cm2d.x, cm2d.y, 0);
}

//------------------------------------------------------------------------
const Vec3D& ModelNBody::GetCamDir() const
{
  return m_camDir;
}

//------------------------------------------------------------------------
const Vec3D& ModelNBody::GetCamPos() const
{
  return m_camPos;
}

//------------------------------------------------------------------------
double* ModelNBody::GetInitialState()
{
  return reinterpret_cast<double*>(m_pInitial);
}

//------------------------------------------------------------------------
double ModelNBody::GetMass(int i)
{
  return m_pAux[i].mass;
}

//------------------------------------------------------------------------
int ModelNBody::GetMethod()
{
  return methode_calcul;
}
//------------------------------------------------------------------------
/** \brief Calcule la vitesse orbitale de la particule 2 par rapport à la particule 1*/
void ModelNBody::GetOrbitalVelocity(const ParticleData &p1, const ParticleData &p2)
{
  double x1 = p1.m_pState->x,
         y1 = p1.m_pState->x,
         m1 = p1.m_pAuxState->mass;
  double x2 = p2.m_pState->x,
         y2 = p2.m_pState->y;

  // Calculate distance from the planet with index idx_main
  double r[2], dist;
  r[0] = x1 - x2;
  r[1] = y1 - y2;


  // distance in parsec
  dist = sqrt(r[0] * r[0] + r[1] * r[1]);

  // Based on the distance from the sun calculate the velocity needed to maintain a circular orbit
  double v = sqrt(gamma_1 * m1 / dist);

  // Calculate a suitable vector perpendicular to r for the velocity of the tracer
  double &vx = p2.m_pState->vx,
         &vy = p2.m_pState->vy;

  vx =( r[1] / dist) * v*100;
  vy =(-r[0] / dist) * v*100;
}

//------------------------------------------------------------------------
/** \brief Réinitialise un objet ModelNBodr*/
void ModelNBody::ResetDim(int num, double stepsize)
{
  m_num = num;
  SetDim(m_num*4);

  //Création du tableau de particules (position et vitesse)
  delete m_pInitial;  
  m_pInitial = new PODState[num];

  //Création du tableau de particules (masse)
  delete m_pAux;  
  m_pAux = new PODAuxState[num];

  //Delta t de la simulation
  m_timeStep = stepsize;

  // reset bounding box and center
  m_max.x = m_max.y = std::numeric_limits<double>::min();
  m_min.x = m_min.y = std::numeric_limits<double>::max();
  m_center = Vec2D(0,0);    // for storing the center of mass
}

//---------------------------------------------------------------------------------------------------- GALAXIE SPIRALE
/** \brief Initialise 2 black holes et des particules aléatoires autour des blackhole*/
void ModelNBody::InitCollision(int num, int mode)
{

  if(mode==0)
  {
    // Reset model size
    //Création des tableaux de particules +,le delta t
    ResetDim(num, m_timeStep);

    // initialize particles
    ParticleData blackHole;

    for (int i=0; i<m_num; ++i)
    {
      //Itération dans le tableau de particules
      PODState &st        = m_pInitial[i];
      PODAuxState &st_aux = m_pAux[i];

      if (i==0)
      {
        //On relie la particule balck hole 1 avec st[0] et st_aux[0]
        blackHole.m_pState = &st;
        blackHole.m_pAuxState = &st_aux;

        //Paramètres de black hole 1
        st.x  = st.y = 0;
        st.vx = st.vy = 0;      //vitesse nulle
        st_aux.mass = 4310;  // 4.31 Millionen Sonnenmassen
      }
      //---------------------------------------------------------------------------------------------------- POUR DÉCIDER À PARTIR DE QUELLE PARTICULE ON CHANGE DE GALAXIE
      else  //Initialisation des autres particules  random autour de black hole 1
      {
        const double rad = 100;
        double r = 0.1 + .8 * (rad * ((double)rand() / RAND_MAX))+1;
        double a = 2.0*M_PI*((double)rand() / RAND_MAX);
        st_aux.mass =1 + 2 * ((double)rand() / RAND_MAX);
        st.x = r*sin(a);
        st.y = r*cos(a);

        GetOrbitalVelocity(blackHole, ParticleData(&st, &st_aux));
      }

      // determine the size of the area including all particles
      m_max.x = std::max(m_max.x, st.x);
      m_max.y = std::max(m_max.y, st.y);
      m_min.x = std::min(m_min.x, st.x);
      m_min.y = std::min(m_min.y, st.y);
    }

    // The Barnes Hut algorithm needs square shaped quadrants.
    // calculate the height of the square including all particles (and a bit more space)
    double l = 1.05 * std::max(m_max.x - m_min.x,
                               m_max.y - m_min.y);

    m_roi = l * 1.5;

    // compute the center of the region including all particles
    Vec2D c(m_min.x + (m_max.x - m_min.x)/2.0,
            m_min.y + (m_max.y - m_min.y)/2.0);
    m_min.x = c.x - l/2.0;
    m_max.x = c.x + l/2.0;
    m_min.y = c.y - l/2.0;
    m_max.y = c.y + l/2.0;

    std::cout << "Initial particle distribution area\n";
    std::cout << "----------------------------------\n";
    std::cout << "Particle spread:\n";
    std::cout << "  xmin=" << m_min.x << ", ymin=" << m_min.y << "\n";
    std::cout << "  xmax=" << m_max.y << ", ymax=" << m_max.y << "\n";
    std::cout << "Bounding box:\n";
    std::cout << "  cx =" << c.x   << ", cy  =" << c.y   << "\n";
    std::cout << "  l  =" << l << "\n";
  }

  else if(mode==1)
  {
    // Reset model size
    //Création des tableaux de particules +,le delta t
    ResetDim(num, m_timeStep);

    // initialize particles
    ParticleData blackHole;

    float angle = 0;
    for (int i=0; i<m_num; ++i)
    {
      //Itération dans le tableau de particules
      PODState &st        = m_pInitial[i];
      PODAuxState &st_aux = m_pAux[i];

      angle += M_PI/4;

      if(i%8==0)
      {
        angle = 0;
      }

      if (i==0)
      {
        //On relie la particule balck hole 1 avec st[0] et st_aux[0]
        blackHole.m_pState = &st;
        blackHole.m_pAuxState = &st_aux;

        //Paramètres de black hole 1
        st.x  = st.y = 0;
        st.vx = st.vy = 0;      //vitesse nulle
        st_aux.mass = 431;  // 4.31 Millionen Sonnenmassen
      }
      //---------------------------------------------------------------------------------------------------- POUR DÉCIDER À PARTIR DE QUELLE PARTICULE ON CHANGE DE GALAXIE
      else  //Initialisation des autres particules  random autour de black hole 1
      {
        double r = i*0.0001+0.8*((double)rand()/RAND_MAX);
        st_aux.mass =0.03 + 20 * ((double)rand() / RAND_MAX);

        st.x = r*sin(angle);
        st.y = r*cos(angle);

        GetOrbitalVelocity(blackHole, ParticleData(&st, &st_aux));
      }

      // determine the size of the area including all particles
      m_max.x = std::max(m_max.x, st.x);
      m_max.y = std::max(m_max.y, st.y);
      m_min.x = std::min(m_min.x, st.x);
      m_min.y = std::min(m_min.y, st.y);
    }

    // The Barnes Hut algorithm needs square shaped quadrants.
    // calculate the height of the square including all particles (and a bit more space)
    double l = 1.05 * std::max(m_max.x - m_min.x,
                               m_max.y - m_min.y);

    m_roi = l * 1.5;

    // compute the center of the region including all particles
    Vec2D c(m_min.x + (m_max.x - m_min.x)/2.0,
            m_min.y + (m_max.y - m_min.y)/2.0);
    m_min.x = c.x - l/2.0;
    m_max.x = c.x + l/2.0;
    m_min.y = c.y - l/2.0;
    m_max.y = c.y + l/2.0;

    std::cout << "Initial particle distribution area\n";
    std::cout << "----------------------------------\n";
    std::cout << "Particle spread:\n";
    std::cout << "  xmin=" << m_min.x << ", ymin=" << m_min.y << "\n";
    std::cout << "  xmax=" << m_max.y << ", ymax=" << m_max.y << "\n";
    std::cout << "Bounding box:\n";
    std::cout << "  cx =" << c.x   << ", cy  =" << c.y   << "\n";
    std::cout << "  l  =" << l << "\n";
  }

  else if(mode==2)
  {
    // Reset model size
    //Création des tableaux de particules +,le delta t
    ResetDim(num, m_timeStep);

    // initialize particles
    ParticleData blackHole;
    ParticleData blackHole2;

    for (int i=0; i<m_num; ++i)
    {
      //Itération dans le tableau de particules
      PODState &st        = m_pInitial[i];
      PODAuxState &st_aux = m_pAux[i];

      if (i==0)
      {
        //On relie la particule balck hole 1 avec st[0] et st_aux[0]
        blackHole.m_pState = &st;
        blackHole.m_pAuxState = &st_aux;

        //Paramètres de black hole 1
        st.x  = st.y = 0;
        st.vx = st.vy = 0;      //vitesse nulle
        st_aux.mass = 431000;  // 4.31 Millionen Sonnenmassen
      }
      //---------------------------------------------------------------------------------------------------- POUR DÉCIDER À PARTIR DE QUELLE PARTICULE ON CHANGE DE GALAXIE
      else if (i<num/2)  //Initialisation des autres particules  random autour de black hole 1
      {
        const double rad = 10;
        double r = 0.1 + .8 * (rad * ((double)rand() / RAND_MAX))+1;
        double a = 2.0*M_PI*((double)rand() / RAND_MAX);
        st_aux.mass =0.03 + 20 * ((double)rand() / RAND_MAX);
        st.x = r*sin(a);
        st.y = r*cos(a);

        GetOrbitalVelocity(blackHole, ParticleData(&st, &st_aux));
      }
      else if (i==num/2)   //Initialisation de black hole 2
      {
        blackHole2.m_pState = &st;
        blackHole2.m_pAuxState = &st_aux;

        //Paramètres de black hole 2
        st.x = st.y = 10;
        st_aux.mass = 1000;

        //Vitesse orbitale de black hole 2 par rapport à black hole 1
        GetOrbitalVelocity(blackHole, blackHole2);

        //Vitesse orbitale non nulle de black hole 2 par rapport à black hole 1 
        blackHole2.m_pState->vx *= 0.9;
        blackHole2.m_pState->vy *= 0.9;

      }
      else
      {
        //Calculs des paramètres des particules random de black hole 2
        const double rad = 3;
        double r = 0.1 + .8 *  (rad * ((double)rand() / RAND_MAX));
        double a = 2.0*M_PI*((double)rand() / RAND_MAX);
        st_aux.mass = 3 + 20 * ((double)rand() / RAND_MAX);
        st.x = blackHole2.m_pState->x + r*sin(a);
        st.y = blackHole2.m_pState->y + r*cos(a);

        //Calcule de la vitesse orbitale des particules random autour de black hole 2
        GetOrbitalVelocity(blackHole2, ParticleData(&st, &st_aux));

        //On rajoute la vitesse de black hole 2 aux vitesses
        st.vx+=blackHole2.m_pState->vx;
        st.vy+=blackHole2.m_pState->vy;
      }

      // determine the size of the area including all particles
      m_max.x = std::max(m_max.x, st.x);
      m_max.y = std::max(m_max.y, st.y);
      m_min.x = std::min(m_min.x, st.x);
      m_min.y = std::min(m_min.y, st.y);
    }

    // The Barnes Hut algorithm needs square shaped quadrants.
    // calculate the height of the square including all particles (and a bit more space)
    double l = 1.05 * std::max(m_max.x - m_min.x,
                               m_max.y - m_min.y);

    m_roi = l * 1.5;

    // compute the center of the region including all particles
    Vec2D c(m_min.x + (m_max.x - m_min.x)/2.0,
            m_min.y + (m_max.y - m_min.y)/2.0);
    m_min.x = c.x - l/2.0;
    m_max.x = c.x + l/2.0;
    m_min.y = c.y - l/2.0;
    m_max.y = c.y + l/2.0;

    std::cout << "Initial particle distribution area\n";
    std::cout << "----------------------------------\n";
    std::cout << "Particle spread:\n";
    std::cout << "  xmin=" << m_min.x << ", ymin=" << m_min.y << "\n";
    std::cout << "  xmax=" << m_max.y << ", ymax=" << m_max.y << "\n";
    std::cout << "Bounding box:\n";
    std::cout << "  cx =" << c.x   << ", cy  =" << c.y   << "\n";
    std::cout << "  l  =" << l << "\n";
  }
  else if(mode==3){
    // Reset model size
    //Création des tableaux de particules +,le delta t
    ResetDim(num, m_timeStep);

    // initialize particles
    ParticleData blackHole;
    ParticleData blackHole2;

    for (int i=0; i<m_num; ++i)
    {
      //Itération dans le tableau de particules
      PODState &st        = m_pInitial[i];
      PODAuxState &st_aux = m_pAux[i];

      if (i==0)
      {
        //On relie la particule balck hole 1 avec st[0] et st_aux[0]
        blackHole.m_pState = &st;
        blackHole.m_pAuxState = &st_aux;

        //Paramètres de black hole 1
        st.x  = st.y = 0;
        st.vx = st.vy = 0;    //vitesse nulle
        st_aux.mass = 4310000;//4310000; //43100;//1000000; //431000;   // 4.31 Millionen Sonnenmassen
      }
      //---------------------------------------------------------------------------------------------------- POUR DÉCIDER À PARTIR DE QUELLE PARTICULE ON CHANGE DE GALAXIE
      else if (i<70000)  //Initialisation des autres particules  random autour de black hole 1
      {
        const double rad = 10;
        double r = 0.1 + .8 * (rad * ((double)rand() / RAND_MAX))+1;
        double a = 2.0*M_PI*((double)rand() / RAND_MAX);
        st_aux.mass =0.03 + 20 * ((double)rand() / RAND_MAX);
        st.x = r*sin(a);
        st.y = r*cos(a);

        GetOrbitalVelocity(blackHole, ParticleData(&st, &st_aux));
      }
      else if (i==70000)   //Initialisation de black hole 2
      {
        blackHole2.m_pState = &st;
        blackHole2.m_pAuxState = &st_aux;

        //Paramètres de black hole 2
        st.x = st.y = 10;
        st_aux.mass = 100000;

        //Vitesse orbitale de black hole 2 par rapport à black hole 1
        GetOrbitalVelocity(blackHole, blackHole2);

        //Vitesse orbitale non nulle de black hole 2 par rapport à black hole 1 
        blackHole2.m_pState->vx *= 0.9;
        blackHole2.m_pState->vy *= 0.9;

      }
      else
      {
        //Calculs des paramètres des particules random de black hole 2
        const double rad = 3;
        double r = 0.1 + .8 *  (rad * ((double)rand() / RAND_MAX));
        double a = 2.0*M_PI*((double)rand() / RAND_MAX);
        st_aux.mass = 0.03 + 20 * ((double)rand() / RAND_MAX);
        st.x = blackHole2.m_pState->x + r*sin(a);
        st.y = blackHole2.m_pState->y + r*cos(a);

        //Calcule de la vitesse orbitale des particules random autour de black hole 2
        GetOrbitalVelocity(blackHole2, ParticleData(&st, &st_aux));

        //On rajoute la vitesse de black hole 2 aux vitesses
        st.vx+=blackHole2.m_pState->vx;
        st.vy+=blackHole2.m_pState->vy;
      }

      // determine the size of the area including all particles
      m_max.x = std::max(m_max.x, st.x);
      m_max.y = std::max(m_max.y, st.y);
      m_min.x = std::min(m_min.x, st.x);
      m_min.y = std::min(m_min.y, st.y);
    }

    // The Barnes Hut algorithm needs square shaped quadrants.
    // calculate the height of the square including all particles (and a bit more space)
    double l = 1.05 * std::max(m_max.x - m_min.x,
                              m_max.y - m_min.y);

    m_roi = l * 1.5;

    // compute the center of the region including all particles
    Vec2D c(m_min.x + (m_max.x - m_min.x)/2.0,
            m_min.y + (m_max.y - m_min.y)/2.0);
    m_min.x = c.x - l/2.0;
    m_max.x = c.x + l/2.0;
    m_min.y = c.y - l/2.0;
    m_max.y = c.y + l/2.0;

    std::cout << "Initial particle distribution area\n";
    std::cout << "----------------------------------\n";
    std::cout << "Particle spread:\n";
    std::cout << "  xmin=" << m_min.x << ", ymin=" << m_min.y << "\n";
    std::cout << "  xmax=" << m_max.y << ", ymax=" << m_max.y << "\n";
    std::cout << "Bounding box:\n";
    std::cout << "  cx =" << c.x   << ", cy  =" << c.y   << "\n";
    std::cout << "  l  =" << l << "\n";
  }

}

//------------------------------------------------------------------------------
/** \brief A coder lorsque que l'on fera l'algorithme de Barnes Hut*/
void ModelNBody::CalcBHArea(const ParticleData &data)
{

  // reset bounding box
  m_max.x = m_max.y = std::numeric_limits<double>::min();
  m_min.x = m_min.y = std::numeric_limits<double>::max();

  //Détermine les positions de la boîte contenant toutes les particules
  #pragma omp parallel for
  for (int i=0; i<m_num; ++i)
  {
    PODState &s = data.m_pState[i];

    // determine the size of the area including all particles
    m_max.x = std::max(m_max.x, s.x);
    m_max.y = std::max(m_max.y, s.y);
    m_min.x = std::min(m_min.x, s.x);
    m_min.y = std::min(m_min.y, s.y);
  }

  // The Barnes Hut algorithm needs square shaped quadrants.
  // calculate the height of the square including all particles (and a bit more space)
  double l = 1.05 * std::max(m_max.x - m_min.x,
                             m_max.y - m_min.y);

  // compute the center of the region including all particles
  Vec2D c(m_min.x + (m_max.x - m_min.x)/2.0,
          m_min.y + (m_max.y - m_min.y)/2.0);
  m_min.x = c.x - l/2.0;
  m_max.x = c.x + l/2.0;
  m_min.y = c.y - l/2.0;
  m_max.y = c.y + l/2.0;

}

//------------------------------------------------------------------------------
/** \brief A recoder lorsqu'on fera l'algorithme de BH (Build the barnes hut tree by adding all particles that are inside
           the region of interest.)
*/
void ModelNBody::BuiltTree(const ParticleData &all)
{
  // Reset the quadtree, make sure only particles inside the roi
  // are handled. The renegade ones may live long and prosper
  // outside my simulation
  m_root.Reset(Vec2D(m_center.x - m_roi-20, m_center.y - m_roi-20),
               Vec2D(m_center.x + m_roi+20, m_center.y + m_roi+20));

  // build the quadtree
  int ct = 0;

  // #pragma omp parallel for
  for (int i=0; i<m_num; ++i)
  {
//    PODState *st = &(all.m_pState[i]);

    try
    {
      // extract data for a single particle
      ParticleData p(&(all.m_pState[i]),
                     &(all.m_pAuxState[i]));

      // insert the particle, but only if its inside the roi
      m_root.Insert(p,0);
      // cout <<"i= "<<  i <<"\n";
      ++ct;
    }

    catch(std::exception &exc)
    {
/*
      std::cout << exc.what() << "\n";
      std::cout << "Particle " << i << " (" << st->x << ", " << st->y << ") is outside the roi (skipped).\n";
      std::cout << "  roi size   =   " << m_roi << "\n";
      std::cout << "  roi center = (" << m_center.x << ", " << m_center.y << ")\n";
*/
    }
  }

  // compute masses and center of mass on all scales of the tree
  m_root.ComputeMassDistribution();

  if (m_bVerbose)
  {
    std::cout << "Tree Dump\n";
    std::cout << "---------\n";
    m_root.DumpNode(-1, 0);
    std::cout << "\n\n";
  }

  // update the center of mass
  m_center = m_root.GetCenterOfMass();
}

//------------------------------------------------------------------------
const PODAuxState* ModelNBody::GetAuxState() const
{
  return m_pAux;
}

//------------------------------------------------------------------------
BHTreeNode* ModelNBody::GetRootNode()
{
  return &m_root;
}

//------------------------------------------------------------------------
int ModelNBody::GetN() const
{
  return m_num;
}

//------------------------------------------------------------------------
double ModelNBody::GetTheta() const
{
  return m_root.GetTheta();
}

//------------------------------------------------------------------------
void ModelNBody::SetVerbose(bool bVerbose)
{
  m_bVerbose = bVerbose;
}

//------------------------------------------------------------------------
void ModelNBody::SetTheta(double theta)
{
  m_root.SetTheta(theta);
}

//------------------------------------------------------------------------
double ModelNBody::GetTimeUnit() const
{
  return time_1;
}

//------------------------------------------------------------------------
void ModelNBody::Fusion(double *a_state)
{
  /*Fusionne les particules qui sont à la mêmes position*/

  PODState *pState = reinterpret_cast<PODState*>(a_state);

  double dx=0;
  double dy=0;
  double r=0;

  for (int i=1; i<m_num; ++i){

    for(int j=1; j<m_num;j++){
      dx=(pState+i)->x - (pState+j)->x;    //distance en x
      dy=(pState+i)->y - (pState+j)->y;    //distance en y
      r=(double)sqrt(dx*dx + dy*dy);    //distance 

      //Masses
      double m= (m_pAux+i)->mass + (m_pAux+j)->mass;
      double m1= (m_pAux+i)->mass;  //particule actuelle
      double m2= (m_pAux+j)->mass;  //particule à supprimer

      //vitesses
      double v1x = (pState+i)->vx;
      double v2x = (pState+j)->vx;
      double v1y = (pState+i)->vy;
      double v2y = (pState+j)->vy;
      // cout << "r=" << r << "\n";
      if(r<0.1 && m_pAux[i].mass>0 && m_pAux[j].mass>0  && i!=j){
        m_pAux[i].mass = m;
        //En conservant la quantité de mouvement des particules à fusionner
        (pState+i)->vx = v1x *(m1*m1)/m + v2x*(m2*m2)/m;
        (pState+i)->vy = v1y *(m1*m1)/m + v2y*(m2*m2)/m;

        //En conservant l'énergie cinétique des particules à fusionner
        // (pState+i)->vx = sqrt(v1x*v1x * m1/m + v2x*v2x * m2/m);
        // (pState+i)->vy = sqrt(v1y*v1y * m1/m + v2y*v2y * m2/m);

        //suppression de la particule
        m_pAux[j].mass = 0;
        (pState+j)->vx=0;
        (pState+j)->vy=0;
        (pState+j)->x=0;
        (pState+j)->y=0;
      }

    }
  }

}

//------------------------------------------------------------------------
/** \brief Calcule l'accélération et la vitesse pour chaque particule*/

void ModelNBody::Eval(double *a_state, double a_time, double *a_deriv)
{
  // wrap the complete particle data together for easier treatment
  // in the following algorithms

  PODState *pState = reinterpret_cast<PODState*>(a_state);
  PODDeriv *pDeriv = reinterpret_cast<PODDeriv*>(a_deriv);
  ParticleData all(pState, m_pAux);

  CalcBHArea(all);

  //Construction de l'arbre et mesure de temps
  gettimeofday(&start,NULL);
  BuiltTree(all);
  gettimeofday(&fin,NULL);
  mesure_construction_BH  += time_diff(&start,&fin);

  //Tableaux pour le calcul efficace
  double acx[m_num];
  double acy[m_num];

  //initialisation des tableaux pour le calcul efficace
  #pragma omp parallel for
  for(int j=0; j<m_num; j++){
   acx[j]=0;
   acy[j]=0;
  }

  gettimeofday(&start, NULL);

  //Nombre de threads à utiliser en région parallèle
  // omp_set_num_threads(8);

  //Calcul des forces
  #pragma omp parallel for
  for (int i=0; i<m_num; ++i)   //i=1 de base
  {

    //Affichage du nombre de threads
    //std:: cout << "num threads =" <<omp_get_num_threads()<< endl;

    ParticleData p(&pState[i], &m_pAux[i]);

    Vec2D acc;

    //Calcul Barnes Hut
    if(methode_calcul==1)
      acc = m_root.CalcTreeForce(p);

    //Calcul efficace
    else if(methode_calcul==2)
      acc = m_root.CalcNaiveForceefficient(p, pState, m_pAux,m_num,i,acx,acy);

    //Calcul non efficace
    else
      acc = m_root.CalcNaiveForce(p, pState, m_pAux,m_num,i);  


    //Multiplication par la constante gravitationnelle
    acc.x *= gamma_si;
    acc.y *= gamma_si;

    //Stockage des forces et vitesses
    pDeriv[i].ax = acc.x;
    pDeriv[i].ay = acc.y;
    pDeriv[i].vx = pState[i].vx;
    pDeriv[i].vy = pState[i].vy;

  }

  // Particle 0 is calculated last, because the statistics
  // data relate to this particle. They would be overwritten
  // otherwise
  m_root.StatReset(); 

  ParticleData p(&pState[0], &m_pAux[0]);

  Vec2D acc;

  //Calcul Barnes-Hut
  if(methode_calcul==1){
    acc.x=0;
    acc.y=0;
    acc = m_root.CalcTreeForce(p);
  }

  //Version efficace
  else if(methode_calcul==2){
    acc.x =acx[0]/m_pAux[0].mass;
    acc.y =acy[0]/m_pAux[0].mass;
  }

  //Version non efficace
  else
    acc = m_root.CalcNaiveForce(p, pState, m_pAux,m_num,0);
  
  //Multiplication par la constante gravitationnelle
  acc.x *= gamma_si;
  acc.y *= gamma_si;

  //Stockage des forces et vitesses
  pDeriv[0].ax = acc.x;
  pDeriv[0].ay = acc.y;
  pDeriv[0].vx = pState[0].vx;
  pDeriv[0].vy = pState[0].vy;

  // Save vectors for camera orientations
  m_camDir.x = pState[0].x - pState[4000].x;
  m_camDir.y = pState[0].y - pState[4000].y;
  m_camPos.x = 0;
  m_camPos.y = 0;
  m_camPos.x = m_root.GetCenterOfMass().x;
  m_camPos.y = m_root.GetCenterOfMass().y;

  //Anti-warning
  acc.x = acx[0];
  acc.y = acy[0];

  gettimeofday(&fin, NULL);
  mesure_temps  += time_diff(&start, &fin);
}

//------------------------------------------------------------------------
bool ModelNBody::IsFinished(double *state)
{
  return false;
}

