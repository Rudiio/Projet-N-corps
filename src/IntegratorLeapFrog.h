#ifndef _INTEGRATOR_LEAPFROG_H
#define _INTEGRATOR_LEAPFROG_H

//------------------------------------------------------------------------------
#include <string>

//------------------------------------------------------------------------------
#include "IIntegrator.h"


//------------------------------------------------------------------------------
/** \brief Inplementation of the simple Euler integration scheme. */
class IntegratorLeapFrog : public IIntegrator
{
public:

  IntegratorLeapFrog(IModel *pModel, double h);
  virtual ~IntegratorLeapFrog();
  virtual void SingleStep();
  virtual void SetInitialState(double *state);
  virtual double* GetState() const;

private:

  double *m_state;
  double *statehalf;
  double *statehalfacc;

};


#endif // include guard