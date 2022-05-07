#include "BHTree.h"

//--- Standard includes --------------------------------------------------------
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>


#include <SDL/SDL.h>


#define softening 0.5
//------------------------------------------------------------------------------
// static variables
double BHTreeNode::s_theta = 0.9;
std::vector<ParticleData> BHTreeNode::s_renegades;
BHTreeNode::DebugStat BHTreeNode::s_stat = {0};
double BHTreeNode::s_gamma = 0;       // gravitational constant is set from the outside
double BHTreeNode::s_soft = 0.1*0.1;        // approx. 3 light year

//------------------------------------------------------------------------------

BHTreeNode::BHTreeNode(const Vec2D &min,
                       const Vec2D &max,
                       BHTreeNode *parent)
  :m_particle()
  ,m_mass(0)
  ,m_cm()
  ,m_min(min)
  ,m_max(max)
  ,m_center(min.x+(max.x-min.x)/2.0, min.y+(max.y-min.y)/2.0)
  ,m_parent(parent)
  ,m_num(0)
  ,m_bSubdivided(false)
{
  s_stat.m_nNumCalc =0;
  m_quadNode[0] = m_quadNode[1] = m_quadNode[2] = m_quadNode[3] = NULL;
}

//------------------------------------------------------------------------------
bool BHTreeNode::IsRoot() const
{
  return m_parent==NULL;
}

//------------------------------------------------------------------------------
bool BHTreeNode::IsExternal() const
{
  return  m_quadNode[0]==NULL &&
          m_quadNode[1]==NULL &&
          m_quadNode[2]==NULL &&
          m_quadNode[3]==NULL;
}

//------------------------------------------------------------------------------
bool BHTreeNode::WasTooClose() const
{
  return m_bSubdivided;
}

//------------------------------------------------------------------------------
const Vec2D& BHTreeNode::GetMin() const
{
  return m_min;
}

//------------------------------------------------------------------------------
const Vec2D& BHTreeNode::GetMax() const
{
  return m_max;
}

//------------------------------------------------------------------------------
const Vec2D& BHTreeNode::GetCenterOfMass() const
{
  return m_cm;
}

//------------------------------------------------------------------------------
double BHTreeNode::GetTheta() const
{
  return s_theta;
}

//------------------------------------------------------------------------------
void BHTreeNode::SetTheta(double theta)
{
  s_theta = theta;
}

//------------------------------------------------------------------------------
int BHTreeNode::StatGetNumCalc() const
{
  return s_stat.m_nNumCalc;
}

//------------------------------------------------------------------------------
/** \brief Returns the number of particles not assigned to any node. */
int BHTreeNode::GetNumRenegades() const
{
  return s_renegades.size();
}

//------------------------------------------------------------------------------
/** \brief Returns the number of particles inside this node. */
int BHTreeNode::GetNum() const
{
  return this->m_num;
}

//------------------------------------------------------------------------------
void BHTreeNode::StatReset()
{
  if (!IsRoot())
    throw std::runtime_error("Only the root node may reset statistics data.");

  s_stat.m_nNumCalc = 0;

  struct ResetSubdivideFlags
  {
    ResetSubdivideFlags(BHTreeNode *pRoot)
    {
      ResetFlag(pRoot);
    }

    void ResetFlag(BHTreeNode *pNode)
    {
      pNode->m_bSubdivided = false;
      for (int i=0; i<4;++i)
      {
        if (pNode->m_quadNode[i])
          ResetFlag(pNode->m_quadNode[i]);
      }
    }
  } ResetFlagNow(this);
}

//------------------------------------------------------------------------------
void BHTreeNode::Reset(const Vec2D &min,
                       const Vec2D &max)
{
  if (!IsRoot())
    throw std::runtime_error("Only the root node may reset the tree.");

  for (int i=0; i<4; ++i)
  {
    delete m_quadNode[i];
    m_quadNode[i] = NULL;
  }

  m_min = min;
  m_max = max;
  m_center = Vec2D(min.x + (max.x-min.x)/2.0,
                   min.y + (max.y-min.y)/2.0);
  m_num = 0;
  m_mass = 0;
  m_cm = Vec2D(0, 0);

  s_renegades.clear();
}

//------------------------------------------------------------------------------
BHTreeNode::EQuadrant BHTreeNode::GetQuadrant(double x, double y) const
{
  /* Retourne le quadrant de l'espace correspondant */

  //Ouest
  if(x <= m_center.x ){ 
    if(y >= m_center.y)   // Sud-Ouest
      return SW;
    else          // Nord-Ouest
      return NW;
  }
  //Est
  else if(x>m_center.x){
    if(y >= m_center.y)   //Sud-Est
      return SE;
    else        //Nord-Est
      return NE;
  }
  else
    return NONE;
}

//------------------------------------------------------------------------------
BHTreeNode* BHTreeNode::CreateQuadNode(EQuadrant eQuad)
{
  /*Crée le noeud correspondant au quadrant */

  if(eQuad == NW) // Création du noeud Nord Ouest
    return new BHTreeNode(m_min,m_center, this);
  else if(eQuad == SW)  // Création du noeud sud Ouest
    return new BHTreeNode(Vec2D(m_min.x,m_center.y), Vec2D(m_center.x,m_max.y), this);
  else if (eQuad == NE) // Création du noeud Nord Est
    return new BHTreeNode(Vec2D(m_center.x,m_min.y), Vec2D(m_max.x,m_center.y), this);
  else if(eQuad == SE)  // Création du noeud Nord Est
    return new BHTreeNode(m_center,m_max, this);
  else
    return NULL;
}


void BHTreeNode:: CreatequadTree()
{
  /*Crée tous les fils du noeud actuel*/

  m_quadNode[NE]=CreateQuadNode(NE);
  m_quadNode[NW]=CreateQuadNode(NW);
  m_quadNode[SW]=CreateQuadNode(SW);
  m_quadNode[SE]=CreateQuadNode(SE);

}

//------------------------------------------------------------------------------
/* \brief Calcule la masse et le centre de masse de chaque noeud de l'arbre */
void BHTreeNode::ComputeMassDistribution()
{
  if(m_num==1){
      m_mass = m_particle.m_pAuxState->mass;
      m_cm.x = m_particle.m_pState->x;
      m_cm.y = m_particle.m_pState->y;
}
  else {
    m_mass=0;
    m_cm.x=0;
    m_cm.y=0;
    for(int i=0; i<4; i++){
      if(m_quadNode[i]!=NULL){
        m_quadNode[i]->ComputeMassDistribution();
        m_mass += m_quadNode[i]->m_mass;
        m_cm.x += m_quadNode[i]->m_cm.x*m_quadNode[i]->m_mass;
        m_cm.y += m_quadNode[i]->m_cm.y*m_quadNode[i]->m_mass;
      }   
    }
    m_cm.x = m_cm.x/m_mass;
    m_cm.y = m_cm.y/m_mass; 
  }

}

//------------------------------------------------------------------------------
/** \brief Calculate the acceleration caused by gravitaion of p2 on p1. */
Vec2D BHTreeNode::CalcAcc(const ParticleData &p1, const ParticleData &p2) const
{
  Vec2D acc;
  
  double &x1(p1.m_pState->x);
  double &y1(p1.m_pState->y);
  double &x2(p2.m_pState->x);
  double &y2(p2.m_pState->y);

  double &m2(p2.m_pAuxState->mass);

  double r = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

  acc.x = (x2-x1)*m2/(r*r*r + softening);
  acc.y = (y2-y1)*m2/(r*r*r + softening);

  return acc;
}

//------------------------------------------------------------------------------
Vec2D BHTreeNode::CalcForce(const ParticleData &p1) const
{
  Vec2D acc;
  return acc;
}

Vec2D BHTreeNode::CalcNaiveForce(const ParticleData &p1, PODState *pState, PODAuxState *m_pAux,int p_num,int pos) const
{
  Vec2D acc;
  acc.x=0;
  acc.y=0;

  double dx=0;
  double dy=0;
  double r3=0;
  int i=0;
  if(m_pAux[pos].mass > 0.0){
    for(i = 0;i<pos;i++){
      if(m_pAux[i].mass > 0.0){
        s_stat.m_nNumCalc += 4;
        dx=pState[i].x - p1.m_pState->x;    //distance en x
        dy=pState[i].y - p1.m_pState->y;    //distance en y
        r3=std::pow(dx*dx + dy*dy,1.5);    //distance au cube

        acc.x += (m_pAux[i].mass * dx)/(r3+softening);
        acc.y += (m_pAux[i].mass *dy)/(r3+softening);
      }
        
    }

    for(i = pos+1;i<p_num;i++){
      if(m_pAux[i].mass > 0.0){
        s_stat.m_nNumCalc += 4;
        dx=p1.m_pState->x - pState[i].x;    //distance en x
        dy=p1.m_pState->y - pState[i].y;    //distance en y
        r3=std::pow(dx*dx + dy*dy,1.5);    //distance au cube
        

        acc.x += (m_pAux[i].mass * dx)/(r3+softening);
        acc.y += (m_pAux[i].mass *dy)/(r3+softening);
      }
    }
  }

  return acc;
}

Vec2D BHTreeNode::CalcNaiveForceefficient(const ParticleData &p1, PODState *pState, PODAuxState *m_pAux,int p_num,int pos,double *acx,double *acy) const
{
  Vec2D acc;
  acc.x=0;
  acc.y=0;

  double dx=0;
  double dy=0;
  double r3=0;
  int i=0;
  
  if(m_pAux[pos].mass > 0.0){

    for(i = pos;i<p_num;i++){
      if(m_pAux[i].mass > 0.0){
        //Nombre de calculs effectués
        s_stat.m_nNumCalc += 4;

        //Calculs des variations de distances
        dx=pState[i].x - p1.m_pState->x;    //distance en x
        dy=pState[i].y - p1.m_pState->y;    //distance en y
        r3=std::pow(dx*dx + dy*dy,1.5);    //distance au cube

        double fx = (m_pAux[i].mass*m_pAux[pos].mass * dx)/(r3+softening);
        double fy = (m_pAux[i].mass*m_pAux[pos].mass *dy)/(r3+softening);
        //force réciproque pour les autres particules
        acx[i]-=fx;
        acy[i]-=fy;

        //Addition des forces pour la particules p
        acx[pos]+=fx;
        acy[pos]+=fy;

      }
    }
  }

    acc.x = acx[pos]/m_pAux[pos].mass;
    acc.y = acy[pos]/m_pAux[pos].mass;

  return acc;
}

//------------------------------------------------------------------------------
/**  \brief Compute the force acting from this node and it's child
            to a particle p.
*/
Vec2D BHTreeNode::CalcTreeForce(const ParticleData &p1) const
{
  Vec2D acc;

  double &x1(p1.m_pState->x);
  double &y1(p1.m_pState->y);
  //double &x2(p2.m_pState->x);
  //double &y2(p2.m_pState->y);

  if(m_num==1) {

    acc = CalcAcc(p1, m_particle);
    s_stat.m_nNumCalc++;

  } else {

    double d = sqrt( (x1 - m_cm.x)*(x1 - m_cm.x) + (y1 - m_cm.y)*(y1 - m_cm.y) );
    double D = m_max.x - m_min.x;

    if(D/d < s_theta) {

      m_bSubdivided = false;
      acc.x = m_mass * (m_cm.x - x1) / (d*d*d + softening);
      acc.y = m_mass * (m_cm.y - y1) / (d*d*d + softening);

    } else {

      m_bSubdivided = true;
      Vec2D vect;
      for(int i=0; i<4; i++) {
          if(m_quadNode[i]!=NULL) {
            vect = m_quadNode[i]->CalcTreeForce(p1);
            acc.x += vect.x;
            acc.y += vect.y;
          }
      }

    }

  }

  return acc;
}

//------------------------------------------------------------------------------
void BHTreeNode::DumpNode(int quad, int level)
{
  std::string space;
  for (int i=0; i<level; ++i)
      space+= "  ";

  std::cout << space << "Quadrant " << quad << ": ";
  std::cout << space << "(num=" << m_num << "; ";
  std::cout << space << "mass=" << m_mass << ";";
  std::cout << space << "cx=" << m_cm.x << ";";
  std::cout << space << "cy=" << m_cm.y << ")\n";

  for (int i=0; i<4;++i)
  {
    if (m_quadNode[i])
    {
      m_quadNode[i]->DumpNode(i, level+1);
    }
  }
}

//------------------------------------------------------------------------------
/** \brief Ajoute une particule dans l'arbre */
void BHTreeNode::Insert(const ParticleData &newParticle,int level)
{
  /* Insère une nouvelle particule dans l'arbre */

  //On vérifie si la particule rentre dans l'arbre
  if(m_min.x <newParticle.m_pState->x && newParticle.m_pState->x < m_max.x && m_min.y <newParticle.m_pState->y && newParticle.m_pState->y < m_max.y ){
    
    //Il y a plus d'une particule dans le noeud 
    if(m_num > 1){
      EQuadrant c = GetQuadrant(newParticle.m_pState->x,newParticle.m_pState->y);
      //Création du fils
      if(!m_quadNode[c])
        m_quadNode[c]= CreateQuadNode(c);
      m_quadNode[c]->Insert(newParticle,level+1);
      m_num++;
    }

    //Il y a une particule dans l'arbres
    else if(m_num == 1){
      //Création des fils de l'arbres
      // CreatequadTree();
      
      //On replace l'ancienne particule
      EQuadrant c = GetQuadrant(m_particle.m_pState->x,m_particle.m_pState->y);
      //Création du fils
      m_quadNode[c] = CreateQuadNode(c);
      m_quadNode[c]->Insert(m_particle,level+1);

      //On place la nouvelle particule
      EQuadrant d = GetQuadrant(newParticle.m_pState->x,newParticle.m_pState->y);
      //Création du fils
      if(!m_quadNode[d])
        m_quadNode[d] = CreateQuadNode(d);
      m_quadNode[d]->Insert(newParticle,level+1);
      m_num=2;
    }

    //Le noeud est une feuille
    else if(m_num==0){
      m_particle = ParticleData(newParticle);
      m_num=1;
    }
  }
}

//------------------------------------------------------------------------------
BHTreeNode::~BHTreeNode()
{
  for (int i=0; i<4; ++i)
    delete m_quadNode[i];
}
