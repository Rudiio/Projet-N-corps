#ifndef _TYPES_H
#define _TYPES_H


#pragma pack(push, 1)

//---------------------------------------------------------------------
/** \brief struct de la position et la vitesse d'une particule */
struct PODState
{
  double x;
  double y;
  double vx;
  double vy;
};

//---------------------------------------------------------------------
/** \brief struct de la masse d'une particule */
struct PODAuxState
{
  double mass;
};

//---------------------------------------------------------------------
/** \brief struct de la vitesse et de l'accélération*/
struct PODDeriv
{
  double vx;
  double vy;
  double ax;
  double ay;
};

#pragma pack(pop)

//---------------------------------------------------------------------
/** \brief Classe d'une particule */
struct ParticleData
{
  ParticleData();
  ParticleData(PODState *m_pState, PODAuxState *m_pAuxState);
  ParticleData(const ParticleData &ref);
  ParticleData& operator=(const ParticleData &ref);

  void Reset();
  bool IsNull() const;

  PODState *m_pState;
  PODAuxState *m_pAuxState;
};

#endif // _TYPES_H
