/*
 * File:   NBodyWnd.cpp
 * Author: user
 *
 * Created on 8. Juli 2009, 21:54
 */
#include "NBodyWnd.h"

//--- Standard icludes ---------------------------------------------------------
#include <iostream>
#include <cmath>
#include <cassert>
#include <limits>
#include <omp.h>
#include <iostream>

using namespace std;
//------------------------------------------------------------------------------
#include "IntegratorEuler.h"
#include "IntegratorLeapFrog.h"

//------------------------------------------------------------------------------
/** \brief Constructeur NBodyWnd */
NBodyWnd::NBodyWnd(int sz, std::string caption)
  :SDLWindow(sz, sz, 30.0, caption)
  ,m_pModel(NULL)
  ,m_pSolver(NULL)
  ,m_camOrient(0)
  ,m_flags(dspBODIES | /*dspTREE |*/ dspAXIS | dspSTAT | dspVERBOSE)
  ,m_bDumpState(false)
  ,m_bDumpImage(true)
{}

//------------------------------------------------------------------------------
/** \brief Destructeur */
NBodyWnd::~NBodyWnd()
{
  delete m_pModel;
  delete m_pSolver;
  delete colorspread;
}

//------------------------------------------------------------------------------
void NBodyWnd::Init(int num, int methode_calcul,int mode_init)
{
  // Clean stuff
  delete m_pModel;
  delete m_pSolver;

  // Create the n-body model class
  m_pModel = new ModelNBody(num,methode_calcul,mode_init);
  std::cout <<"timestep =" << m_pModel->GetSuggestedTimeStep() << "\n";
  
  // Assign model to the solver and set the integration step width
  m_pSolver = new IntegratorLeapFrog(m_pModel, m_pModel->GetSuggestedTimeStep()); // Construction de l'intégrateur avec un pointeur vers le modèle
  m_pSolver->SetInitialState(m_pModel->GetInitialState());

  //--------------------------------------------------------------------------------------------------- POUR RÉGLER LES COULEURS
  // Create colorspread
  colorspread = new int[num];
  for(int i=0; i<num; i++)
  {
    if(m_pModel->GetMass(i)<2)
      colorspread[i] = 0;
    else if(m_pModel->GetMass(i)<5)
      colorspread[i] = 1;
    else
      colorspread[i] = 2;
  }

  if (m_bDumpState)
  {
    std::string sName = std::string("traces_barnes_hut_") + m_pSolver->GetID() + ".dat";
    m_outfile.open(sName.c_str());
  }

  // OpenGL initialization
  glClear(GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT);

  SetCameraOrientation(Vec3D(0, 1, 0));

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

//------------------------------------------------------------------------------
void NBodyWnd::Render()
{
  static int ct = 0;
  if (!(m_flags & dspPAUSE))
  {
    m_pSolver->SingleStep();

    ++ct;
    if (m_bDumpState && ct%1000000==0)
    {
      int num = m_pModel->GetDim()/4;
      PODState *state = reinterpret_cast<PODState*>(m_pSolver->GetState());
      m_outfile << m_pSolver->GetTime() << ", ";
      for (int i=0; i<num; ++i)
      {
        m_outfile << state[i].x << ", "
                  << state[i].y << ", ";
      }
      m_outfile << std::endl;
    }


    // if (m_bDumpImage && ct%2000==0)
    // {
    //   SaveToTGA();
    // }

    // the bottleneck is the output, reduce load a bit
//    if (ct%10000!=0)
//      return;
  }

//Passer en commentaire tout ce qui est en dessous pour supprimer l'affichage
  glClear(GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT);

  Vec3D orient;
  switch (m_camOrient)
  {
  case 0: orient.x = 0;
          orient.y = 1;
          orient.z = 0;
          break;

  case 1:
          orient = m_pModel->GetCamDir();
          break;
  }
  SetCameraOrientation(orient);

  if (m_flags & dspAXIS)
  {
    // Draw axis at position of the first particle
//    PODState *state = reinterpret_cast<PODState*>(m_pSolver->GetState());
    // const Vec3D &cm = m_pModel->GetCenterOfMass();
    DrawAxis(Vec2D(0, 0));
  }

  if (m_flags & dspTREE)
    DrawTree();

  if (m_flags & dspBODIES)
    DrawBodies();

  if (m_flags & dspROI)
    DrawROI();

  if (m_flags & dspSTAT)
    DrawStat();

  if (m_flags & dspHELP)
    DrawHelp();

  SDL_GL_SwapBuffers();
}

//------------------------------------------------------------------------------
void NBodyWnd::DrawBodies()
{
  assert(m_pSolver);

  PODState *state = reinterpret_cast<PODState*>(m_pSolver->GetState());
  //  const PODAuxState *state_aux  = m_pModel->GetAuxState();

  //glColor3f(1,1,1);
  glPointSize(2);
  glBegin(GL_POINTS);

  for(int i=0; i<m_pModel->GetN(); i++)
  {
    if(colorspread[i]==0) {
      glColor3f(0,0.8,0.8);
      //glPointSize(1);
    } else if(colorspread[i]==1) {
      glColor3f(0,01.7,0.4);
      //glPointSize(2);
    } else {
      glColor3f(0,0.6,0.6);
      //glPointSize(2);      
    }
    glVertex3f(state[i].x, state[i].y, 0.0f);
  }

  glEnd();
}

//------------------------------------------------------------------------------
void NBodyWnd::DrawStat()
{
  double x0 = 10, y0 = 20, dy = 20;
  int line = 0;

  // acquire pointer to the root node of the barnes hut tree from the model
  BHTreeNode *pRoot = m_pModel->GetRootNode();

  const Vec3D &camPos = GetCamPos(),
//              &camOrient = GetCamOrient(),
              &camLookAt = GetCamLookAt();
  glColor3f(1, 1, 1);
  TextOut(x0, y0 + dy * line++, "Number of bodies (outside tree): %d (%d)", pRoot->GetNum(), pRoot->GetNumRenegades());
  TextOut(x0, y0 + dy * line++, "Theta: %2.1f", pRoot->GetTheta());
  TextOut(x0, y0 + dy * line++, "FPS: %d", GetFPS());
  TextOut(x0, y0 + dy * line++, "Time: %2.1f y", m_pSolver->GetTime());
  TextOut(x0, y0 + dy * line++, "Camera: %2.2f, %2.2f, %2.2f", camPos.x, camPos.y, camPos.z);
  TextOut(x0, y0 + dy * line++, "LookAt: %2.2f, %2.2f, %2.2f", camLookAt.x, camLookAt.y, camLookAt.z);
  TextOut(x0, y0 + dy * line++, "Field of view: %2.2f pc", GetFOV());
  TextOut(x0, y0 + dy * line++, "Calculations: %d", pRoot->StatGetNumCalc());
  TextOut(x0, y0 + dy * line++, "Solver: %s", m_pSolver->GetID().c_str());
  if(m_pModel->GetMethod()==1)
    TextOut(x0, y0 + dy * line++, "Method: Barnes-Hut");
  else if(m_pModel->GetMethod()==2)
    TextOut(x0, y0 + dy * line++, "Method: Naive optimised");
  else 
    TextOut(x0, y0 + dy * line++, "Method: Naive");
}

//------------------------------------------------------------------------------
void NBodyWnd::DrawHelp()
{
  double x0 = 10, y0 = 20, dy = 20;
  int line = 0;
  Vec3D p;

  glColor3f(1, 1, 1);
  TextOut(x0, y0 + dy * line++, "Keyboard commands");
  TextOut(x0, y0 + dy * line++, "a     - display axis");
  TextOut(x0, y0 + dy * line++, "t     - display Barnes Hut tree");
  TextOut(x0, y0 + dy * line++, "s     - display statistic data");
  TextOut(x0, y0 + dy * line++, "c     - display center of mass");
  TextOut(x0, y0 + dy * line++, "b     - display particles");
  TextOut(x0, y0 + dy * line++, "h     - display help text");
  TextOut(x0, y0 + dy * line++, "0..1  - Set camera orientation");
  TextOut(x0, y0 + dy * line++, "pause - halt simulation");
}

//------------------------------------------------------------------------------
void NBodyWnd::DrawROI()
{
  const Vec3D &cm = m_pModel->GetCenterOfMass();

  double l = GetFOV()/20;
  glColor3f(1, 0, 0);

  glPushMatrix();
  glTranslatef(cm.x, cm.y, cm.z);
  // cross at the center of mass
  glBegin(GL_LINES);
    glVertex3f(-l,  0, 0);
    glVertex3f( l,  0, 0);
    glVertex3f( 0, -l, 0);
    glVertex3f( 0,  l, 0);
  glEnd();

  // region of interest
  l = m_pModel->GetROI()/2.0;
  glBegin(GL_LINE_STRIP);
    glVertex3f( -l,  l, 0);
    glVertex3f(  l,  l, 0);
    glVertex3f(  l, -l, 0);
    glVertex3f( -l, -l, 0);
    glVertex3f( -l,  l, 0);
  glEnd();

  glPopMatrix();
}

//------------------------------------------------------------------------------
void NBodyWnd::DrawTree()
{
  struct DrawTree
  {
    enum EWhat
    {
      COMPLETE,  ///< Display all nodes
      APPROX,    ///< Display only the nodes as they where used for the force calculation
    };

    DrawTree(BHTreeNode *pRoot, EWhat what, int flags, int sz)
      :displayFlags(flags)
    {
      DrawNode(pRoot, 0, what, sz);
    }

    void DrawNode(BHTreeNode *pNode, int level, EWhat what, int sz)
    {
      assert(pNode);

      double col = 1 - level*0.2;
      switch (what)
      {
        case COMPLETE:  glColor3f(col, 1, col); break;
        case APPROX:  glColor3f(0, 1, 0); break;
      }

      // Draw node rectangle
      if ( what==COMPLETE ||
          (what==APPROX  && !pNode->WasTooClose()) )
      {
        const Vec2D &min = pNode->GetMin(),
                    &max = pNode->GetMax();
        glBegin(GL_LINE_STRIP);
           glVertex3f(min.x, min.y, 0);
           glVertex3f(max.x, min.y, 0);
           glVertex3f(max.x, max.y, 0);
           glVertex3f(min.x, max.y, 0);
           glVertex3f(min.x, min.y, 0);
        glEnd();

        // Draw a cross at the center of mass is the corresponding flag is set
        if (displayFlags & dspCENTER_OF_MASS && !pNode->IsExternal())
        {
          double len = (double)sz/50 * std::max(1 - level*0.2, 0.1);
          glPointSize(4);
          glColor3f(col, 1, col);

          const Vec2D cm = pNode->GetCenterOfMass();
          glBegin(GL_LINES);
            glVertex3f(cm.x-len, cm.y, 0);
            glVertex3f(cm.x+len, cm.y, 0);
          glEnd();
          glBegin(GL_LINES);
            glVertex3f(cm.x, cm.y-len, 0);
            glVertex3f(cm.x, cm.y+len, 0);
          glEnd();
        }
      }

      // If the node was not subdivided in the force calculation
      // dont draw its subnodes either
      if (what!=COMPLETE && !pNode->WasTooClose())
        return;

      for (int i=0; i<4; ++i)
      {
        if (pNode->m_quadNode[i])
          DrawNode(pNode->m_quadNode[i], level+1, what, sz);
      }
    } // DrawTree::DrawNode

    int displayFlags;
  }; // struct DrawTree

  BHTreeNode *pRoot = m_pModel->GetRootNode();
  if ( (m_flags & dspTREE) && (m_flags & dspTREE_COMPLETE) )
    DrawTree DrawComplete(pRoot, DrawTree::COMPLETE, m_flags, GetFOV());
  else if ( (m_flags & dspTREE) && !(m_flags & dspTREE_COMPLETE) )
    DrawTree DrawFar(pRoot, DrawTree::APPROX, m_flags, GetFOV());
}

//------------------------------------------------------------------------------
void NBodyWnd::DrawNode(BHTreeNode *pNode, int level)
{
  assert(pNode);
  double len = 0.01 * std::max(1 - level*0.2, 0.1);

  double col = 1 - level*0.2;

//  glColor3f(col, 1, col);

  if (pNode->WasTooClose())
    glColor3f(1, col, col);
  else
    glColor3f(col, 1, col);

  const Vec2D &min = pNode->GetMin(),
              &max = pNode->GetMax();
  glBegin(GL_LINE_STRIP);
    glVertex3f(min.x, min.y, 0);
    glVertex3f(max.x, min.y, 0);
    glVertex3f(max.x, max.y, 0);
    glVertex3f(min.x, max.y, 0);
    glVertex3f(min.x, min.y, 0);
  glEnd();

  if (m_flags & dspCENTER_OF_MASS && !pNode->IsExternal())
  {
    Vec2D cm = pNode->GetCenterOfMass();

    glPointSize(4);
    glColor3f(col, 1, col);
    glBegin(GL_LINES);
      glVertex3f(cm.x-len, cm.y,     0);
      glVertex3f(cm.x+len, cm.y,     0);
      glVertex3f(cm.x,     cm.y-len, 0);
      glVertex3f(cm.x,     cm.y+len, 0);
    glEnd();
  }

  if (!pNode->WasTooClose())
    return;

  for (int i=0; i<4; ++i)
  {
    if (pNode->m_quadNode[i])
    {
      DrawNode(pNode->m_quadNode[i], level+1);
    }
  }
}

//------------------------------  double len = 0.01 * std::max(1 - level*0.2, 0.1);------------------------------------------------
void NBodyWnd::OnProcessEvents(uint8_t type)
{
  switch (type)
  {
    case SDL_KEYDOWN:
        switch (m_event.key.keysym.sym)
        {
          case SDLK_i:
               m_camOrient = 0;
               break;

          case SDLK_j:
               m_camOrient = 1;
               break;

          case SDLK_a:
                std::cout << "Display:  Toggling axis " << ((m_flags & dspAXIS) ? "off" : "on") << "\n";
                m_flags ^= dspAXIS;
                break;

          case  SDLK_b:
                std::cout << "Display:  Toggling bodies " << ((m_flags & dspBODIES) ? "off" : "on") << "\n";
                m_flags ^= dspBODIES;
                break;

          case  SDLK_t:
                {
                  if (!(m_flags & dspTREE))
                  {
                    m_flags ^= dspTREE;
                    std::cout << "Display:  Tree cells used in force calculation\n";
                  }
                  else if ( (m_flags & dspTREE) && !(m_flags & dspTREE_COMPLETE))
                  {
                    m_flags ^= dspTREE_COMPLETE;
                    std::cout << "Display:  Complete tree\n";
                  }
                  else if ( (m_flags & dspTREE) && (m_flags & dspTREE_COMPLETE))
                  {
                    m_flags &= ~(dspTREE | dspTREE_COMPLETE);
                    std::cout << "Display:  No tree\n";
                  }
                }
                break;

          case  SDLK_c:
                std::cout << "Display:  Center of mass " << ((m_flags & dspCENTER_OF_MASS) ? "off" : "on") << "\n";
                m_flags ^= dspCENTER_OF_MASS;
                break;

          case  SDLK_h:
                std::cout << "Display:  Help text" << ((m_flags & dspHELP) ? "off" : "on") << "\n";
                m_flags ^= dspHELP;
                m_flags &= ~dspSTAT;
                break;

          case  SDLK_PAUSE:
                std::cout << "Simulation:  pause " << ((m_flags & dspPAUSE) ? "off" : "on") << "\n";
                m_flags ^= dspPAUSE;
                break;

          case  SDLK_v:
                std::cout << "Simulation:  verbose mode " << ((m_flags & dspVERBOSE) ? "off" : "on") << "\n";
                m_flags ^= dspVERBOSE;
                if (m_pModel)
                  m_pModel->SetVerbose(m_flags & dspVERBOSE);
                break;

          case  SDLK_s:
                std::cout << "Display:  statistics " << ((m_flags & dspSTAT) ? "off" : "on") << "\n";
                m_flags ^= dspSTAT;
                break;

          case  SDLK_f:
                std::cout << "Display:  force indicator " << ((m_flags & dspARROWS) ? "off" : "on") << "\n";
                m_flags ^= dspARROWS;
                break;

          case SDLK_r:
               std::cout << "Display:  region of interest " << ((m_flags & dspROI) ? "off" : "on") << "\n";
               m_flags ^= dspROI;
               break;

          case SDLK_p:
               SaveToTGA();
               break;

          case SDLK_y:
               m_pModel->SetTheta(m_pModel->GetTheta() + 0.1);
               break;

          case SDLK_x:
               m_pModel->SetTheta(std::max(m_pModel->GetTheta() - 0.1, 0.1));
               break;

          case SDLK_k:
               ScaleAxis(0.9);
               SetCameraOrientation(Vec3D(0,1,0));
               break;

          case SDLK_l:
               ScaleAxis(1.1);
               SetCameraOrientation(Vec3D(0,1,0));
               break;

          case SDLK_ESCAPE:
              ExitMainLoop();
              break;
              
          default:
               break;

        }

        break;
  }
}

//------------------------------------------------------------------------------
float NBodyWnd::gettime()
{
  return m_pModel->Getmesuretemps();
}

//------------------------------------------------------------------------------
float NBodyWnd::getconstruction()
{
  return m_pModel->Getmesureconstruction();
}

//------------------------------------------------------------------------------
void NBodyWnd::CalcultateEnergy(int methode)
{
  //Position et vitesse actuelle
  PODState *pState = reinterpret_cast<PODState*>(m_pSolver->GetState());
  int num = m_pModel->GetN();
  double gamma = 6.67428e-11;
  Ecin = Epot = 0;
  
  for(int i=0;i<num;i++){
    double vx = pState[i].vx;
    double vy = pState[i].vy;

    //Clalcul de l'énergie cinétique
    Ecin += (double)1/2 * m_pModel->GetMass(i)* (vx*vx + vy*vy);


    for(int j=0;j<num;j++){
      if(i!=j){
        //Calcul de l'énergie potentielle
        double dx = pState[i].x - pState[j].x;  
        double dy = pState[i].y - pState[j].y; 
        double d = sqrt(dx*dx + dy*dy + 0.5);
        Epot += gamma * m_pModel->GetMass(i)*m_pModel->GetMass(j)/(d);
      }
    }
  }

  //Ecriture dans les fichiers
  //Barne-Hut
  if(methode==1){
    ofstream fichier("./data/E_BH.txt", ios_base::app ); 
    if(fichier.is_open()){
      fichier << m_pSolver->GetTime() << " " << Ecin << " " << Epot << " " << Ecin -Epot << endl;
      fichier.close();
    }
  }

  //Naïve optimisée
  else if(methode==2){
    std::ofstream fichier("./data/E_NO.txt", ios_base::app ); 
    if(fichier.is_open()){
      fichier << m_pSolver->GetTime() << " " << Ecin << " " << Epot << " " << Ecin -Epot << endl;
      fichier.close();
    }
  }
  //Naïve
  else if(methode==3){
    std::ofstream fichier("./data/E_N.txt", ios_base::app ); 
    if(fichier.is_open()){
      fichier << m_pSolver->GetTime() << " " << Ecin << " " << Epot << " " << Ecin -Epot << endl;
      fichier.close();
    }
  }
}