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





int main(){

    time_t start=time(NULL);
    ROOT::EnableImplicitMT(24);
    cout.precision(2);

    CSTCutG PidCut;
    TString fileName = "/home/sghwang/workspace/SmallHIMACData_v5.6.root";
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
    int EventList[6] = {0};
    vector<int> hasEvent;
    vector<int> hasParticle;
    
    for(int event = 0; event < TotalEventNum ; event++){
        if(event%2346==0)cout <<fixed<< "    Event Cut Progress = " << event<< " / " << TotalEventNum << " ( " << (double)event/(double)TotalEventNum*100.<<" % )           \r" << flush;
        treein->GetEntry(event);

        if(siHitNum!=1)continue;

        int rightEvent = 0;
        int leftEvent  = 0;

        for(int wing=0; wing <2; wing++){
            double CsI_down_energy = csiData[wing][0][0];
            double CsI_up_energy   = csiData[wing][1][0];

            if(CsI_down_energy<0.1 && CsI_up_energy<0.1 )continue;
        
            for(int hitnum =0 ; hitnum < siHitNum ; hitnum++){
                double Si_omic_energy = siData[wing][hitnum][0][0];
                double Si_junc_energy = siData[wing][hitnum][0][1];
                double Si_energy      = Si_junc_energy+Si_omic_energy;
                if(Si_energy < 0.1)continue;

                double si_Y_position  = siData[wing][hitnum][2][1];
                
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
    TCanvas *c1 = new TCanvas("c1","c1",2400,800);
    c1->Divide(2,1);
    TGraph *PlotCanvas = new TGraph();
    TGraph *PlotCanvas_xy = new TGraph();

    PlotCanvas->SetPoint(0,-10,-40);
    PlotCanvas->SetPoint(1, 40,-40);
    PlotCanvas->SetPoint(2, 40, 40);
    PlotCanvas->SetPoint(3,-10, 40);
    PlotCanvas->SetLineWidth(0);
    PlotCanvas->SetLineColor(0);

    PlotCanvas_xy->SetPoint(0,-100,-100);
    PlotCanvas_xy->SetPoint(1, 100,-100);
    PlotCanvas_xy->SetPoint(2, 100, 100);
    PlotCanvas_xy->SetPoint(3,-100, 100);
    PlotCanvas_xy->SetLineWidth(0);
    PlotCanvas_xy->SetLineColor(0);

    TH1D *targetZ = new TH1D("zpos","",35,-100,100);
    // TH1D *targetZ_right = new TH1D("zpos_r","zpos from each detector",35,-100,100);
    // TH1D *targetZ_left = new TH1D("zpos_l","zpos from left detector",35,-100,100);
    TH1D *targetZ_right = new TH1D("zpos_r","",35,-100,100);
    TH1D *targetZ_left = new TH1D("zpos_l","",35,-100,100);
    TF1 *fitgaus = new TF1("fitgaus","gaus");


    int valid_event  = 0;
    int valid_event2  = 0;
    for(int event = 0; event <hasEvent.size(); event++){
        if(event%23==0)cout <<fixed<< "    Track Calculation Progress = " << event<< " / " << hasEvent.size() << " ( " << (double)event/(double)hasEvent.size()*100.<<" % )           \r" << flush;
        
        treein->GetEntry(hasEvent[event]);
        if(trackFitData[0][0][1][0]<-500 || trackFitData[0][0][1][1] <-500 || trackFitData[1][0][1][0] < -500 || trackFitData[1][0][1][1]<-500)continue;

        double right_slope_yz = trackFitData[0][0][1][0];
        double right_slope_xy = trackFitData[0][0][0][0];
        double right_const_yz = trackFitData[0][0][1][1];
        double right_const_xy = trackFitData[0][0][0][1];
        double left_slope_yz = trackFitData[1][0][1][0];
        double left_const_yz = trackFitData[1][0][1][1];
        double left_slope_xy = trackFitData[1][0][0][0];
        double left_const_xy = trackFitData[1][0][0][1];

        double atz1 = CSGetTargetZ::getTargetZ(0,right_slope_yz,right_const_yz);
        double atz2 = CSGetTargetZ::getTargetZ(1,left_slope_yz,left_const_yz);

        targetZ->Fill(atz1-atz2);
        targetZ_right->Fill(atz1);
        targetZ_left->Fill(atz2);

        // c1->SaveAs(Form("/home/sghwang/workspace/CS_analysis/pictures/pic_temp%i.png",hasEvent[event]));
        valid_event++;

    }
    c1->cd(1);
    targetZ->Draw();
    targetZ->Fit(fitgaus,"Q");
    c1->cd(2);
    targetZ_right->Draw();
    targetZ_left->SetLineColor(kRed);
    targetZ_left->Draw("same");


    double Zsigma = fitgaus->GetParameter(2);
    double Zsigma_right = targetZ_right->GetStdDev();
    double Zsigma_left = targetZ_left->GetStdDev();

    TCanvas *c2= new TCanvas("c2","c2",800,800);
    c2->cd(1);
    gPad->SetGrid(1);
    TView *view = TView::CreateView(1);
    int irep;
    view->RotateView(90,90);
    // view->SetTheta(90);
    // view->SetPsi(90);
    c2->Update();
    cout << "    Track Calculation Progress = \033[3mDone\033[0m ( 100.0 % )                            " << endl;
    
    TGraph2D *trajectory_plane = new TGraph2D();
    trajectory_plane->SetPoint(0,-500,-500,-50);
    trajectory_plane->SetPoint(1, 500,-500,-50);
    trajectory_plane->SetPoint(2, 500, 500,-50);
    trajectory_plane->SetPoint(3,-500, 500,-50);
    trajectory_plane->SetPoint(4,-500, 500,950);
    trajectory_plane->SetPoint(5, 500, 500,950);
    trajectory_plane->SetPoint(6,-500,-500,950);
    trajectory_plane->SetPoint(7, 500,-500,950);
    trajectory_plane->SetLineColor(0);
    trajectory_plane->SetLineWidth(0);
    trajectory_plane->Draw();
    
    double det_point_x[8] = {-1,63,63,-1,-1,-1,63,63};
    double det_point_y[8] = {-6,-6,90,90,90,-6,-6,90};
    double det_point_z[8] = {0,0,0,0,200,200,200,200};
    
    
    


    for(int event = 0; event <hasEvent.size(); event++){
        if(event%23==0)cout <<fixed<< "    Angle Calculation Progress = " << event<< " / " << hasEvent.size() << " ( " << (double)event/(double)hasEvent.size()*100.<<" % )           \r" << flush;

        TGraph2D *TPCright = new TGraph2D();
        TGraph2D *TPCleft = new TGraph2D();
        TGraph2D *TPCrightMax = new TGraph2D();
        TGraph2D *TPCleftMax = new TGraph2D();
        TPCrightMax -> SetMarkerColor(kRed);
        TPCleftMax -> SetMarkerColor(kRed);

        
        for(int i = 0; i<8; i++){
            double trans_x = CSTransform::Detector_Point(0,det_point_x[i],det_point_y[i],det_point_z[i],1);
            double trans_y = CSTransform::Detector_Point(0,det_point_x[i],det_point_y[i],det_point_z[i],2);
            double trans_z = CSTransform::Detector_Point(0,det_point_x[i],det_point_y[i],det_point_z[i],3);

            int dotN = TPCright->GetN();
            TPCright->SetPoint(dotN,trans_x,trans_y,trans_z);
        
            trans_x = CSTransform::Detector_Point(1,det_point_x[i],det_point_y[i],det_point_z[i],1);
            trans_y = CSTransform::Detector_Point(1,det_point_x[i],det_point_y[i],det_point_z[i],2);
            trans_z = CSTransform::Detector_Point(1,det_point_x[i],det_point_y[i],det_point_z[i],3);
            dotN = TPCleft->GetN();
            TPCleft->SetPoint(dotN,trans_x,trans_y,trans_z);
            
        }
        

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

        double atz1 = CSGetTargetZ::getTargetZ(0,right_slope_yz,right_const_yz);
        double atz2 = CSGetTargetZ::getTargetZ(1,left_slope_yz,left_const_yz);

        if(abs(atz1-atz2)>Zsigma)continue;
        if(abs(atz1)>Zsigma_right)continue;
        if(abs(atz2)>Zsigma_left)continue;
    
        for(int i = -50; i < 600; i ++ ){
            double posZ = i * 1;
            int NpointL = trajectory_left->GetN();
            int NpointR = trajectory_right->GetN();
            double trans_X  = CSTransform::TPC_Point(0,posZ,right_slope_xy,right_const_xy,right_slope_yz,right_const_yz,1);
            double trans_Y  = CSTransform::TPC_Point(0,posZ,right_slope_xy,right_const_xy,right_slope_yz,right_const_yz,2);
            double trans_Z  = CSTransform::TPC_Point(0,posZ,right_slope_xy,right_const_xy,right_slope_yz,right_const_yz,3);
            double trans_X2 = CSTransform::TPC_Point(1,posZ,left_slope_xy,left_const_xy,left_slope_yz,left_const_yz,1);
            double trans_Y2 = CSTransform::TPC_Point(1,posZ,left_slope_xy,left_const_xy,left_slope_yz,left_const_yz,2);
            double trans_Z2 = CSTransform::TPC_Point(1,posZ,left_slope_xy,left_const_xy,left_slope_yz,left_const_yz,3);
            trajectory_right->SetPoint(NpointR,trans_X, trans_Y, trans_Z);
            trajectory_left->SetPoint( NpointL,trans_X2,trans_Y2,trans_Z2);

            //[Right, Left][trkNum][yIdx][xIdx][x, y, z, adc]
            
  
        }


        int maxRowPerLayer[2][8];
        std::fill_n(&maxRowPerLayer[0][0], 16, -1);

        for(int j = 0 ; j < 8 ; j ++){

            double tmpWRight = 0;
            double tmpWLeft = 0;
            for(int i = 0 ; i < 32 ; i++){
            

                double rightx = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2]-238.08,       1);
                double righty = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2]-238.08,       2);
                double rightz = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2]-238.08,3);
                double rightW = trackData[0][0][j][i][3];
                double leftx = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2]-239.45,       1);
                double lefty = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2]-239.45,       2);
                double leftz = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2]-239.45,3);
                double leftW = trackData[1][0][j][i][3];

                int dotN = TPCright->GetN();
                TPCright->SetPoint(dotN,rightx,righty,rightz);
                dotN = TPCleft->GetN();
                TPCleft->SetPoint(dotN,leftx,lefty,leftz);

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

        TPCrightMax -> Set(0);
        TPCleftMax -> Set(0);

        for(int j=0; j<8; j++){

            int i= maxRowPerLayer[0][j]; 
            double rightx = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2]-238.08,       1);
            double righty = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2]-238.08,       2);
            double rightz = CSTransform::Detector_Point(0,trackData[0][0][j][i][0],trackData[0][0][j][i][1],trackData[0][0][j][i][2]-238.08,3);
            int dotN = TPCrightMax->GetN();
            TPCrightMax -> SetPoint(dotN, rightx, righty, rightz);

            i= maxRowPerLayer[1][j]; 
            double leftx = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2]-239.45,       1);
            double lefty = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2]-239.45,       2);
            double leftz = CSTransform::Detector_Point( 1,trackData[1][0][j][i][0],trackData[1][0][j][i][1],trackData[1][0][j][i][2]-239.45,3);
            int dotN2 = TPCleftMax->GetN();
            TPCleftMax -> SetPoint(dotN2, leftx, lefty, leftz);
        }



        TPCright->Draw("same P0");
        TPCleft-> Draw("same P0");
        trajectory_plane->SetTitle(Form("Event #%i",hasEvent[event]));
        trajectory_left->SetLineWidth(2);
        trajectory_right->SetLineWidth(2);
        // trajectory_left->Draw("same p0");
        trajectory_left->Draw("same line");
        // trajectory_right->Draw("same p0");
        trajectory_right->Draw("same line");
        // trajectory_plane->GetYaxis()->SetRangeUser(-100,100);
        
        TPCrightMax -> Draw("same, p0");
        TPCleftMax -> Draw("same, p0");

        
        c2->SaveAs(Form("/home/sghwang/workspace/CS_analysis/pictures/pic_trajec%i.root",event));
        delete TPCright;
        delete TPCleft;
        delete TPCrightMax;
        delete TPCleftMax;

        delete trajectory_left;
        delete trajectory_right;
        // cout << hasEvent[event] << "    " << hasParticle[event] << endl;
        valid_event2++;


 
    }
    
    // c1->SaveAs("/home/sghwang/workspace/CS_analysis/pictures/pic_posZ.png");
    
    cout << "    Angle Calculation Progress = \033[3mDone\033[0m ( 100.0 % )                            " << endl;
    time_t end = time(NULL);
    cout << "    Elapsed Time =  " << end-start << " Seconds"<<endl;
    cout << "    Valid Event : " << valid_event2 << endl;
    cout << "---------------------------------------------------------------------------------------------------"<<endl;

    return 0;
}