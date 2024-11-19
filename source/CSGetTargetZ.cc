#include "CSGetTargetZ.hh"

namespace CSGetTargetZ{
    // ========= ATTPC-constant ===============
    const double fTPCHieghtToBeam = 89.; //[mm]
    const double fTimeoffsetL = 239.45; //[mm] //5.7-245.35 //5.6-239.45 //5.5-233.60 //5.4-227.75 //5.3-221.86 //5.8-251.21
    const double fTimeoffsetR = 238.08; //[mm] //5.7-243.91 //5.6-238.08 //5.5-232.25 //5.4-226.41 //5.3-220.58 //5.8-249.76
    const double fTPCADCThreshold = 50.; // [ADC]
    const double fPadHeight = 11.9; //[mm]
    const double fPadWeith = 1.9; //[mm]
    const double fPadGap = 0.1; //[mm]
    const double fhalfH = 89; //[mm]
    const double fDistance_TPC = 550; //[mm]
    const double fDistance_Si = 715.1; //[mm]
    // ========= ATTPC-constant ===============


    double getTargetZ(int wing, double slope, double constant){
        if(wing == 0){
            return(-(fDistance_TPC-(fPadHeight+fPadGap)*3.5)*slope + constant - fhalfH - fTimeoffsetR)/(Cos(50*Pi()/180));
        }
        else{
            return((fDistance_TPC+(fPadHeight+fPadGap)*3.5)*slope + constant - fhalfH - fTimeoffsetL)/(Cos(50*Pi()/180));
        }
        
    }
    
    double getTargetSi(int wing, double slope, double constant){
        
        TF1 *track = new TF1("track","(x-[1])/[0]");
        track->SetParameter(0,slope);
        track->SetParameter(1,constant);

        double result = 0;

        if(wing == 1){
            result = track->Eval(-(fDistance_Si-fDistance_TPC)-((fPadHeight+fPadGap)*3.5))-(fPadWeith+fPadGap)*15.5;
        }
        else{
            result = track->Eval((fDistance_Si-fDistance_TPC)+((fPadHeight+fPadGap)*3.5))-(fPadWeith+fPadGap)*15.5;
        }
        delete track;
        return result;
    }

    
}