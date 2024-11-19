#include "CSTransform.hh"



namespace CSTransform{
    // ========= ATTPC-constant ===============
    const double fPadHeight = 11.9; //[mm]
    const double fPadWidth = 1.9; //[mm]
    const double fPadGap = 0.1; //[mm]
    const double fhalfH = 89; //[mm]
    const double fDistance_TPC = 550; //[mm]

    const double fTimeoffsetL = 239.16; //[mm] //5.7-245.35 //5.6-239.45 //5.5-233.60 //5.4-227.75 //5.3-221.86 //5.8-251.21
    const double fTimeoffsetR = 238.16; //[mm] //5.7-243.91 //5.6-238.08 //5.5-232.25 //5.4-226.41 //5.3-220.58 //5.8-249.76

    
    // ========= ATTPC-constant ===============

    const double lengthYPadL = (fPadHeight+fPadGap)*3.5 + fDistance_TPC;
    const double lengthYPadR = -(fPadHeight+fPadGap)*3.5 + fDistance_TPC;
    const double lengthXPad = (fPadWidth+fPadGap)*16;
    const double lengthZPad = fhalfH;

    double Detector_Point(bool wing, double xpos, double ypos, double zpos, int xyz){
        point fSlopeToPoint = {0,0,0};

        double Zpoint = zpos;
        double Ypoint = ypos;
        double Xpoint = xpos;
       
        if(wing){
            double trans_x =  Xpoint - lengthXPad;
            double trans_y =  Ypoint - lengthYPadL;
            double trans_z =  Zpoint - lengthZPad ;

            //1st rot ( -50 degree along x axis )
            double rot_x = trans_x ;
            double rot_y = trans_y*cos(-50*TMath::DegToRad()) - trans_z*sin(-50*TMath::DegToRad());
            double rot_z = trans_y*sin(-50*TMath::DegToRad()) + trans_z*cos(-50*TMath::DegToRad());
        
            //2nd rot ( 90 degree along z axis )
            double rot_x2 = rot_x*cos(90*TMath::DegToRad()) - rot_y*sin(90*TMath::DegToRad());
            double rot_y2 = rot_x*sin(90*TMath::DegToRad()) + rot_y*cos(90*TMath::DegToRad());
            double rot_z2 = rot_z;

            point result_point;
            result_point.x = rot_x2;
            result_point.y = rot_y2;
            result_point.z = rot_z2;
            fSlopeToPoint = result_point;
        }
        else{
            double trans_x =  Xpoint - lengthXPad;
            double trans_y =  Ypoint + lengthYPadR;
            double trans_z =  Zpoint - lengthZPad;
            //1st rot ( 50 degree along x axis )
            double rot_x = trans_x ;
            double rot_y = trans_y*cos(50*TMath::DegToRad()) - trans_z*sin(50*TMath::DegToRad());
            double rot_z = trans_y*sin(50*TMath::DegToRad()) + trans_z*cos(50*TMath::DegToRad());
            
            //2nd rot ( 90 degree along z axis )
            double rot_x2 = rot_x*cos(90*TMath::DegToRad()) - rot_y*sin(90*TMath::DegToRad());
            double rot_y2 = rot_x*sin(90*TMath::DegToRad()) + rot_y*cos(90*TMath::DegToRad());
            double rot_z2 = rot_z;

            point result_point;
            result_point.x = rot_x2;
            result_point.y = rot_y2;
            result_point.z = rot_z2;
            fSlopeToPoint = result_point;

        
        }

        if(xyz==1){return fSlopeToPoint.x;}
        else if(xyz==2){return fSlopeToPoint.y;}
        else {return fSlopeToPoint.z;}
        
    }
    point Detector_Point(bool wing, point pos){

        point fSlopeToPoint = {0,0,0};

        double Zpoint = pos.z;
        double Ypoint = pos.y;
        double Xpoint = pos.x;
       
        if(wing){
            double trans_x =  Xpoint - lengthXPad;
            double trans_y =  Ypoint - lengthYPadL;
            double trans_z =  Zpoint - lengthZPad ;
            // double trans_z =  Zpoint - lengthZPad -fTimeoffsetL;

            //1st rot ( -50 degree along x axis )
            double rot_x = trans_x ;
            double rot_y = trans_y*cos(-50*TMath::DegToRad()) - trans_z*sin(-50*TMath::DegToRad());
            double rot_z = trans_y*sin(-50*TMath::DegToRad()) + trans_z*cos(-50*TMath::DegToRad());
        
            //2nd rot ( 90 degree along z axis )
            double rot_x2 = rot_x*cos(90*TMath::DegToRad()) - rot_y*sin(90*TMath::DegToRad());
            double rot_y2 = rot_x*sin(90*TMath::DegToRad()) + rot_y*cos(90*TMath::DegToRad());
            double rot_z2 = rot_z;

            point result_point;
            result_point.x = rot_x2;
            result_point.y = rot_y2;
            result_point.z = rot_z2;
            fSlopeToPoint = result_point;
        }
        else{
            double trans_x =  Xpoint - lengthXPad;
            double trans_y =  Ypoint + lengthYPadR;
            double trans_z =  Zpoint - lengthZPad;
            // double trans_z =  Zpoint - lengthZPad - fTimeoffsetR;
            //1st rot ( 50 degree along x axis )
            double rot_x = trans_x ;
            double rot_y = trans_y*cos(50*TMath::DegToRad()) - trans_z*sin(50*TMath::DegToRad());
            double rot_z = trans_y*sin(50*TMath::DegToRad()) + trans_z*cos(50*TMath::DegToRad());
            
            //2nd rot ( 90 degree along z axis )
            double rot_x2 = rot_x*cos(90*TMath::DegToRad()) - rot_y*sin(90*TMath::DegToRad());
            double rot_y2 = rot_x*sin(90*TMath::DegToRad()) + rot_y*cos(90*TMath::DegToRad());
            double rot_z2 = rot_z;

            point result_point;
            result_point.x = rot_x2;
            result_point.y = rot_y2;
            result_point.z = rot_z2;
            fSlopeToPoint = result_point;

        
        }

        return fSlopeToPoint;

    }
    double TPC_Point(bool wing,double zpos, double slope_xy, double const_xy, double slope_yz, double const_yz,int xyz){
        point fSlopeToPoint = {0,0,0};
        
        // double Ypoint = ypos;
        // double Xpoint = (Ypoint-const_xy)/slope_xy; // y = ax +b
        // double Zpoint = slope_yz*Ypoint + const_yz; // z = ay +b

        double Zpoint = zpos;
        double Ypoint = (Zpoint-const_yz)/slope_yz;
        double Xpoint = (Ypoint-const_xy)/slope_xy;
        


        if(wing){
            double trans_x =  Xpoint - lengthXPad;
            double trans_y =  Ypoint - lengthYPadL;
            double trans_z =  Zpoint - lengthZPad -fTimeoffsetL;

            //1st rot ( -50 degree along x axis )
            double rot_x = trans_x ;
            double rot_y = trans_y*cos(-50*TMath::DegToRad()) - trans_z*sin(-50*TMath::DegToRad());
            double rot_z = trans_y*sin(-50*TMath::DegToRad()) + trans_z*cos(-50*TMath::DegToRad());
        
            //2nd rot ( 90 degree along z axis )
            double rot_x2 = rot_x*cos(90*TMath::DegToRad()) - rot_y*sin(90*TMath::DegToRad());
            double rot_y2 = rot_x*sin(90*TMath::DegToRad()) + rot_y*cos(90*TMath::DegToRad());
            double rot_z2 = rot_z;

            point result_point;
            result_point.x = rot_x2;
            result_point.y = rot_y2;
            result_point.z = rot_z2;
            fSlopeToPoint = result_point;
        }
        else{
            double trans_x =  Xpoint - lengthXPad;
            double trans_y =  Ypoint + lengthYPadR;
            double trans_z =  Zpoint - lengthZPad - fTimeoffsetR;
            //1st rot ( 50 degree along x axis )
            double rot_x = trans_x ;
            double rot_y = trans_y*cos(50*TMath::DegToRad()) - trans_z*sin(50*TMath::DegToRad());
            double rot_z = trans_y*sin(50*TMath::DegToRad()) + trans_z*cos(50*TMath::DegToRad());
            
            //2nd rot ( 90 degree along z axis )
            double rot_x2 = rot_x*cos(90*TMath::DegToRad()) - rot_y*sin(90*TMath::DegToRad());
            double rot_y2 = rot_x*sin(90*TMath::DegToRad()) + rot_y*cos(90*TMath::DegToRad());
            double rot_z2 = rot_z;

            point result_point;
            result_point.x = rot_x2;
            result_point.y = rot_y2;
            result_point.z = rot_z2;
            fSlopeToPoint = result_point;

        
        }

        if(xyz==1){return fSlopeToPoint.x;}
        else if(xyz==2){return fSlopeToPoint.y;}
        else {return fSlopeToPoint.z;}
        
    }
    point TPC_Point(bool wing, double zpos, double slope_xy, double const_xy, double slope_yz, double const_yz){
        point fSlopeToPoint = {0,0,0};
        
        // double Ypoint = ypos;
        // double Xpoint = (Ypoint-const_xy)/slope_xy; // y = ax +b
        // double Zpoint = slope_yz*Ypoint + const_yz; // z = ay +b

        double Zpoint = zpos;
        double Ypoint = (Zpoint-const_yz)/slope_yz;
        double Xpoint = (Ypoint-const_xy)/slope_xy;
        


        if(wing){
            double trans_x =  Xpoint - lengthXPad;
            double trans_y =  Ypoint - lengthYPadL;
            // double trans_z =  Zpoint - lengthZPad;
            double trans_z =  Zpoint - lengthZPad -fTimeoffsetL;

            //1st rot ( -50 degree along x axis )
            double rot_x = trans_x ;
            double rot_y = trans_y*cos(-50*TMath::DegToRad()) - trans_z*sin(-50*TMath::DegToRad());
            double rot_z = trans_y*sin(-50*TMath::DegToRad()) + trans_z*cos(-50*TMath::DegToRad());
        
            //2nd rot ( 90 degree along z axis )
            double rot_x2 = rot_x*cos(90*TMath::DegToRad()) - rot_y*sin(90*TMath::DegToRad());
            double rot_y2 = rot_x*sin(90*TMath::DegToRad()) + rot_y*cos(90*TMath::DegToRad());
            double rot_z2 = rot_z;

            point result_point;
            result_point.x = rot_x2;
            result_point.y = rot_y2;
            result_point.z = rot_z2;
            fSlopeToPoint = result_point;
        }
        else{
            double trans_x =  Xpoint - lengthXPad;
            double trans_y =  Ypoint + lengthYPadR;
            // double trans_z =  Zpoint - lengthZPad ;
            double trans_z =  Zpoint - lengthZPad - fTimeoffsetR;
            //1st rot ( 50 degree along x axis )
            double rot_x = trans_x ;
            double rot_y = trans_y*cos(50*TMath::DegToRad()) - trans_z*sin(50*TMath::DegToRad());
            double rot_z = trans_y*sin(50*TMath::DegToRad()) + trans_z*cos(50*TMath::DegToRad());
            
            //2nd rot ( 90 degree along z axis )
            double rot_x2 = rot_x*cos(90*TMath::DegToRad()) - rot_y*sin(90*TMath::DegToRad());
            double rot_y2 = rot_x*sin(90*TMath::DegToRad()) + rot_y*cos(90*TMath::DegToRad());
            double rot_z2 = rot_z;

            point result_point;
            result_point.x = rot_x2;
            result_point.y = rot_y2;
            result_point.z = rot_z2;
            fSlopeToPoint = result_point;

        
        }

        return fSlopeToPoint;
    }

}