#include "IntegratorLeapFrog.h"

//--- Standard includes --------------------------------------------------------
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include "Types.h"
//------------------------------------------------------------------------------
IntegratorLeapFrog::IntegratorLeapFrog(IModel *pModel, double h)
  :IIntegrator(pModel, h)
  ,m_state(new double[m_dim])
  ,statehalf(new double[m_dim])
  ,statehalfacc(new double[m_dim])
{
  if (pModel==NULL)
    throw std::runtime_error("Model pointer may not be NULL.");

  
  std::stringstream ss;
  ss << "Euler (dt=" << m_h << ")";
  SetID(ss.str());
}

//------------------------------------------------------------------------------
IntegratorLeapFrog::~IntegratorLeapFrog()
{
  delete [] m_state;
  delete [] statehalf;
  delete [] statehalfacc;
}

//------------------------------------------------------------------------------
/** \brief Performs a single integration step. */
void IntegratorLeapFrog::SingleStep()
{

  PODState *pState = reinterpret_cast<PODState*>(m_state);
  PODState *halfState = reinterpret_cast<PODState*>(statehalf);
  PODDeriv *halfStateacc = reinterpret_cast<PODDeriv*>(statehalfacc);

  //Drift
  for(unsigned int i=0;i<m_dim/4;i++){
    halfState[i].x = pState[i].x + pState[i].vx * m_h/2;
    halfState[i].y = pState[i].y + pState[i].vy * m_h/2;
  }

  //Kick
  m_pModel->Eval(statehalf,m_time + m_h/2,statehalfacc);

  for(unsigned int i=0;i<m_dim/4;i++){
    pState[i].vx += m_h * halfStateacc[i].ax;
    pState[i].vy += m_h * halfStateacc[i].ay;
  }

  //Drift
  for(unsigned int i=0;i<m_dim/4;i++){
    pState[i].x = halfState[i].x + m_h/2 * pState[i].vx;
    pState[i].y = halfState[i].y + m_h/2 * pState[i].vy;
  }

  m_time += m_h;

}

//------------------------------------------------------------------------------
/** \brief Sets the initial state of the simulation. */
void IntegratorLeapFrog::SetInitialState(double *state)
{
  for (unsigned i=0; i<m_dim; ++i)
    m_state[i] = state[i];

  m_time = 0;
}

//------------------------------------------------------------------------------
double* IntegratorLeapFrog::GetState() const
{
  return m_state;
}
