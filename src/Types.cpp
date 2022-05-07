#include "Types.h"

//--- Standard includes ------------------------------------------------------------------
#include <cassert>
#include <cstdlib>


//----------------------------------------------------------------------------------------
/** \brief Constructeur à vide de ParticleData */
ParticleData::ParticleData()
  :m_pState(NULL)
  ,m_pAuxState(NULL)
{}

//----------------------------------------------------------------------------------------
/** \brief Constructeur avec  2 valeurs  de ParticleData */
ParticleData::ParticleData(PODState *pState, PODAuxState *pAuxState)
  :m_pState(pState)
  ,m_pAuxState(pAuxState)
{
  assert(m_pState);
  assert(m_pAuxState);
}

//----------------------------------------------------------------------------------------
/** \brief Constructeur avec objet/struct de ParticleData */
ParticleData::ParticleData(const ParticleData &ref)
  :m_pState(ref.m_pState)
  ,m_pAuxState(ref.m_pAuxState)
{}

//----------------------------------------------------------------------------------------
/** \brief overload de = pour ParticleData*/
ParticleData& ParticleData::operator=(const ParticleData &ref)
{
  if (this!=&ref)
  {
    m_pState    = ref.m_pState;
    m_pAuxState = ref.m_pAuxState;
  }

  return *this;
}

//----------------------------------------------------------------------------------------
/** \brief Reset de ParticleData */
void ParticleData::Reset()
{
  m_pState    = NULL;
  m_pAuxState = NULL;
}

//----------------------------------------------------------------------------------------
/** \brief Vérifie si une ParticleData est vide */
bool ParticleData::IsNull() const
{
  return m_pState && m_pAuxState;
}
