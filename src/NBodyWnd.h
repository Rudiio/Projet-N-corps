/*
 * File:   NBodyWnd.h
 * Author: user
 *
 * Created on 8. Juli 2009, 21:54
 */

#ifndef _NBODYWND_H
#define	_NBODYWND_H

#include <stdint.h>
#include <fstream>

//--- Framework ----------------------------------------------------------------
#include "SDLWnd.h"
#include "BHTree.h"
#include "ModelNBody.h"
#include "IIntegrator.h"


//------------------------------------------------------------------------------
/** \brief Main window of th n-body simulation. */
class NBodyWnd : public SDLWindow
{
public:

    NBodyWnd(int sz, std::string caption);
    virtual ~NBodyWnd();

    virtual void Render();
    virtual void OnProcessEvents(uint8_t type);
    virtual float gettime();
    virtual float getconstruction();
    void Init(int num,int methode_calcul);

private:

    enum EDisp
    {
      dspNONE           = 0,
      dspAXIS           = 1 << 0,
      dspBODIES         = 1 << 1,
      dspSTAT           = 1 << 2,
      dspTREE           = 1 << 3,
      dspTREE_COMPLETE  = 1 << 4,
      dspCENTER_OF_MASS = 1 << 5,
      dspPAUSE          = 1 << 6,
      dspVERBOSE        = 1 << 7,
      dspHELP           = 1 << 8,
      dspARROWS         = 1 << 9,
      dspROI            = 1 << 10
    };

    NBodyWnd(const NBodyWnd& orig);

    void DrawBodies();
    void DrawStat();
    void DrawTree();        //Marche pas
    void DrawHelp();
    void DrawROI();
    void DrawCenterOfMass();      //Marche pas
    void DrawNode(BHTreeNode *pNode, int level);  //Marche pas

    ModelNBody *m_pModel;     //Objet : Modèle de résolution
    IIntegrator *m_pSolver;   //Objet : Intégrateur

    int m_camOrient;    ///< Index of the camera orientation to use
    uint32_t m_flags;   ///< The display flags
    bool m_bDumpState;
    bool m_bDumpImage;
    std::ofstream m_outfile;

    int* colorspread;
};

#endif	/* _NBODYWND_H */

