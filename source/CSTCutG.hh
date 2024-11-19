#ifndef CSTCUTG_HH
#define CSTCUTG_HH

#include "TCutG.h"
class CSTCutG {
    public:
    CSTCutG();
    ~CSTCutG();

    int IsIn(double CsIEnergy, double SiEnergy);
    int IsIn(bool isWeird, double CsIEnergy, double SiEnergy);


    private:
    TCutG *PCut;
    TCutG *DCut;
    TCutG *TCut;
    TCutG *wPCut;
    TCutG *wDCut;
    TCutG *wTCut;
    


};




#endif