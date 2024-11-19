#ifndef CSTRANSFORM_HH
#define CSTRANSFORM_HH

#include "TMath.h"
#include "TF1.h"



namespace CSTransform{

    struct point {
        double x;
        double y;
        double z;
    };

 
    double TPC_Point(bool wing, double zpos, double slope_xy, double const_xy, double slope_yz, double const_yz,int xyz);
    point TPC_Point(bool wing, double zpos, double slope_xy, double const_xy, double slope_yz, double const_yz);
    double Detector_Point(bool wing, double xpos, double ypos, double zpos,int xyz);
    point Detector_Point(bool wing, point pos);
}

#endif