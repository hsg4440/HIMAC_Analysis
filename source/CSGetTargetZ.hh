#ifndef CSGETTARGETZ_HH
#define CSGETTARGETZ_HH

#include "TMath.h"
#include "TF1.h"


using namespace TMath;
namespace CSGetTargetZ {
    

    double getTargetZ(int wing, double slope, double constant);
    double getTargetSi(int wing, double slope, double constant);
    
    

}

#endif