//basic header
#include <sstream>
#include <iostream>
#include <time.h>
#include <vector>
#include <algorithm>
//ROOT header
#include "TROOT.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TView.h"
#include "TRotation.h"

//personnal header
#include "CSTCutG.hh"
#include "CSTFileInput.hh"
#include "CSGetTargetZ.hh"
#include "CSTransform.hh"

using namespace std;

const double fTimeoffsetL = 239.16; //[mm] //5.7-245.35 //5.6-239.45 //5.5-233.60 //5.4-227.75 //5.3-221.86 //5.8-251.21
const double fTimeoffsetR = 238.16; //[mm] //5.7-243.91 //5.6-238.08 //5.5-232.25 //5.4-226.41 //5.3-220.58 //5.8-249.76



int main(){

    time_t start=time(NULL);
    ROOT::EnableImplicitMT(24);
    cout.precision(2);

    CSTCutG PidCut;
    // TString fileName = "/home/sghwang/workspace/SmallHIMACData_v5.6.root";
    TString fileName = "/home/sghwang/workspace/forHIMAC/SmallHIMACData.root";
    // TString fileName = "/home/sghwang/workspace/forHIMAC/SmallHIMACDatayc.root";

    CSTFileInput fileIn(fileName);
    TTree* treein = fileIn.getTree();

    int siHitNum;
    double csiData[2][2][2]; //[right,left][down,up][energy,sumADC]
    double siData[2][5][3][2]; //[right,left][hitnum][energy, ADC, position] energy=(omic,junc), ADC=(omic,junc), position(x,y)
    double trackData[2][2][8][32][4]; //[Right, Left][trkNum][yIdx][xIdx][x, y, z, adc]
    treein->SetBranchStatus("*",0);
    treein->SetBranchStatus("siHitNum",1);
    treein->SetBranchStatus("csiData",1);
    treein->SetBranchStatus("siData",1);
    treein->SetBranchAddress("siHitNum",&siHitNum);
    treein->SetBranchAddress("csiData",&csiData);
    treein->SetBranchAddress("siData",&siData);

    int TotalEventNum = treein->GetEntries();
    int valid_event = 0;
    int EventList[6] = {0};
    vector<int> hasEvent;
    vector<int> hasParticle;
    
    for(int event = 0; event < TotalEventNum ; event++){
        if(event%23146==0)cout <<fixed<< "    Event Cut Progress = " << event<< " / " << TotalEventNum << " ( " << (double)event/(double)TotalEventNum*100.<<" % )           \r" << flush;
        treein->GetEntry(event);

        if(siHitNum!=1)continue;

        int rightEvent = 0;
        int leftEvent  = 0;

        for(int wing=0; wing <2; wing++){
            double CsI_down_energy = csiData[wing][0][0];
            double CsI_up_energy   = csiData[wing][1][0];

            if(CsI_down_energy<0.1 && CsI_up_energy<0.1 )continue;
        
            for(int hitnum =0 ; hitnum < siHitNum ; hitnum++){
                double Si_energy = siData[wing][hitnum][0][1]+siData[wing][hitnum][0][0];
                if(Si_energy < 0.1)continue;

                double si_Y_position = siData[wing][hitnum][2][1];
                if(wing == 0){
                    if(si_Y_position>0)rightEvent = PidCut.IsIn(0,CsI_up_energy,Si_energy);
                    else               rightEvent = PidCut.IsIn(0,CsI_down_energy,Si_energy);
                }
                if(wing == 1){
                    if(si_Y_position>0)leftEvent = PidCut.IsIn(1,CsI_up_energy,Si_energy);
                    else               leftEvent = PidCut.IsIn(0,CsI_down_energy,Si_energy);
                }
            }
        }

        if(rightEvent!=0 && leftEvent!=0){
            hasEvent.push_back(event);
        }
        if(rightEvent==1 && leftEvent==1)                                        {
            EventList[0]++;
            hasParticle.push_back(1);
        }
        else if((rightEvent==2 && leftEvent==1)||(rightEvent==1 && leftEvent==2)){
            EventList[1]++;
            hasParticle.push_back(2);
        }
        else if((rightEvent==3 && leftEvent==1)||(rightEvent==1 && leftEvent==3)){
            EventList[2]++;
            hasParticle.push_back(3);
        }
        else if(rightEvent==2 && leftEvent==2)                                   {
            EventList[3]++;
            hasParticle.push_back(4);
        }
        else if((rightEvent==3 && leftEvent==2)||(rightEvent==2 && leftEvent==3)){
            EventList[4]++;
            hasParticle.push_back(5);
        }
        else if(rightEvent==3 && leftEvent==3)                                   {
            EventList[5]++;
            hasParticle.push_back(6);
        }
    }

    cout << "    Event Cut Progress = \033[3mDone\033[0m ( 100.0 % )                                                                        " << endl;
    cout << "    Cutted Both Event : " << hasEvent.size() << endl;
    cout << "    p-p : " << EventList[0] << "    p-d : " << EventList[1] << "    p-t : " <<EventList[2] << "    d-d : " << EventList[3] << "    d-t : " << EventList[4] << "    t-t : " << EventList[5] << endl;

    treein->SetBranchStatus("*",0);
    treein->SetBranchStatus("trackFitData",1);
    treein->SetBranchStatus("trkData", 1);
    treein->SetBranchAddress("trkData", &trackData);
    double trackFitData[2][2][2][3]; //[Right, Left][hitnum][xy,yz][slope,const]
    treein->SetBranchAddress("trackFitData",&trackFitData); 
    gStyle->SetOptStat(0);
  
    TCanvas *c2= new TCanvas("c2","c2",800,800);
    c2->cd(1);
    TView *view = TView::CreateView(1);
    view->RotateView(90,90);


    TGraph2D *trajectory_plane = new TGraph2D();
    trajectory_plane->SetPoint(0,-500,-500,-250);
    trajectory_plane->SetPoint(1, 500,-500,-250);
    trajectory_plane->SetPoint(2, 500, 500,-250);
    trajectory_plane->SetPoint(3,-500, 500,-250);
    trajectory_plane->SetPoint(4,-500, 500, 950);
    trajectory_plane->SetPoint(5, 500, 500, 950);
    trajectory_plane->SetPoint(6,-500,-500, 950);
    trajectory_plane->SetPoint(7, 500,-500, 950);
    trajectory_plane->SetLineColor(0);
    trajectory_plane->SetLineWidth(0);
    trajectory_plane->Draw();

    
    CSTransform::point det_point[8] = { {-1, -6,   0},
                                        {63, -6,   0},
                                        {63, 90,   0},
                                        {-1, 90,   0},
                                        {-1, 90, 200},
                                        {-1, -6, 200},
                                        {63, -6, 200},
                                        {63, 90, 200}};

    

    
    

    for(int event = 0; event <hasEvent.size(); event++){
        if(event%23==0)cout <<fixed<< "    Angle Calculation Progress = " << event<< " / " << hasEvent.size() << " ( " << (double)event/(double)hasEvent.size()*100.<<" % )           \r" << flush;

        TGraph2D *TPCright =new TGraph2D();
        TGraph2D *TPCleft = new TGraph2D();
        TGraph2D *TPCrightMax = new TGraph2D();
        TGraph2D *TPCleftMax = new TGraph2D();

        TPCrightMax -> SetMarkerColor(kRed);
        TPCleftMax -> SetMarkerColor(kRed);

        
        treein->GetEntry(hasEvent[event]);
        TGraph2D *trajectory_right= new TGraph2D();
        TGraph2D *trajectory_left= new TGraph2D();
        if(trackFitData[0][0][1][0]<-500 || trackFitData[0][0][1][1] <-500 || trackFitData[1][0][1][0] < -500 || trackFitData[1][0][1][1]<-500)continue;

        double right_slope_yz = trackFitData[0][0][1][0];
        double right_slope_xy = trackFitData[0][0][0][0];
        double right_const_yz = trackFitData[0][0][1][1];
        double right_const_xy = trackFitData[0][0][0][1];
        double left_slope_yz = trackFitData[1][0][1][0];
        double left_const_yz = trackFitData[1][0][1][1];
        double left_slope_xy = trackFitData[1][0][0][0];
        double left_const_xy = trackFitData[1][0][0][1];

        //drawing tragjectory to the target
        for(int i = -200; i < 600; i ++ ){
            double posZ = i * 1;
            CSTransform::point trans_right = CSTransform::TPC_Point(0,posZ,right_slope_xy,right_const_xy,right_slope_yz,right_const_yz+fTimeoffsetR);
            CSTransform::point trans_left =  CSTransform::TPC_Point(1,posZ,left_slope_xy,left_const_xy,left_slope_yz,left_const_yz+fTimeoffsetL);
            

            int NpointL = trajectory_left->GetN();
            int NpointR = trajectory_right->GetN();
            trajectory_right->SetPoint(NpointR,trans_right.x , trans_right.y , trans_right.z );
            trajectory_left->SetPoint (NpointL,trans_left.x ,  trans_left.y,trans_left.z);

            //[Right, Left][trkNum][yIdx][xIdx][x, y, z, adc]
            
        }

        int maxRowPerLayer[2][8];
        std::fill_n(&maxRowPerLayer[0][0], 16, -1);


        for(int i = 0; i<8; i++){
            CSTransform::point trans_point = CSTransform::Detector_Point(0,det_point[i]);
            double trans_x = trans_point.x;
            double trans_y = trans_point.y;
            double trans_z = trans_point.z;

            int dotN = TPCright->GetN();
            TPCright->SetPoint(dotN,trans_x,trans_y,trans_z);

            trans_point = CSTransform::Detector_Point(1,det_point[i]);
            trans_x = trans_point.x;
            trans_y = trans_point.y;
            trans_z = trans_point.z;

            dotN = TPCleft->GetN();
            TPCleft->SetPoint(dotN,trans_x,trans_y,trans_z);
            
        }
        TPCright->Draw("same P0");
        TPCleft-> Draw("same P0");

        for(int j = 0 ; j < 8 ; j ++){

            double tmpWRight = 0;
            double tmpWLeft = 0;
            for(int i = 0 ; i < 32 ; i++){
                CSTransform::point track_right = {trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2]};
                CSTransform::point transform_track_right = CSTransform::Detector_Point(0,track_right);
                CSTransform::point track_left = {trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2]};
                CSTransform::point transform_track_left = CSTransform::Detector_Point(1,track_left);

                int dotN = TPCright->GetN();
                TPCright->SetPoint(dotN, transform_track_right.x,  transform_track_right.y,  transform_track_right.z  );
                dotN = TPCleft->GetN();
                TPCleft->SetPoint(dotN,  transform_track_left.x,   transform_track_left.y,   transform_track_left.z  );

                double rightW = trackData[0][0][j][i][3];
                double leftW = trackData[1][0][j][i][3];

                if(tmpWRight < rightW){
                    maxRowPerLayer[0][j] = i;
                    tmpWRight = rightW;
                }
                if(tmpWLeft < leftW){
                    maxRowPerLayer[1][j] = i;
                    tmpWLeft = leftW;
                }
            }
        }

       
        
        for(int j=0; j<8; j++){

            int i= maxRowPerLayer[0][j]; 
            double rightx = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2],       1);
            double righty = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2],       2);
            double rightz = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2],       3);
            int dotN = TPCrightMax->GetN();
            TPCrightMax -> SetPoint(dotN, rightx, righty, rightz);

            i= maxRowPerLayer[1][j]; 
            double leftx = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2],       1);
            double lefty = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2],       2);
            double leftz = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2],       3);
            int dotN2 = TPCleftMax->GetN();
            TPCleftMax -> SetPoint(dotN2, leftx, lefty, leftz);
        }

        trajectory_plane->SetTitle(Form("Event #%i",hasEvent[event]));

        trajectory_left->SetLineWidth(2);
        trajectory_right->SetLineWidth(2);
        trajectory_left->Draw("same line");
        trajectory_right->Draw("same line");
        TPCrightMax -> Draw("same, p0");
        TPCleftMax -> Draw("same, p0");

        c2->SaveAs(Form("/home/sghwang/workspace/CS_analysis/pictures/pic_trajec%i.png",event));
        
        delete trajectory_left;
        delete trajectory_right;
        delete TPCright;
        delete TPCleft;
        delete TPCrightMax;
        delete TPCleftMax;
        
        valid_event++;
    }
    
    // c1->SaveAs("/home/sghwang/workspace/CS_analysis/pictures/pic_posZ.png");
    
    cout << "    Angle Calculation Progress = \033[3mDone\033[0m ( 100.0 % )                            " << endl;
    time_t end = time(NULL);
    cout << "    Valid Event : " << valid_event << endl;
    cout << "    Elapsed Time =  " << end-start << " Seconds"<<endl;
    cout << "---------------------------------------------------------------------------------------------------"<<endl;

    return 0;
}
