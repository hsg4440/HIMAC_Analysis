//basic header
#include <sstream>
#include <iostream>
#include <time.h>
#include <vector>
//ROOT header
#include "TROOT.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

//personnal header
#include "CSTCutG.hh"
#include "CSTFileInput.hh"
#include "CSGetTargetZ.hh"

using namespace std;

struct EventData{
    int eventnum;
    int rightwinghitnum;
    int lefttwinghitnum;

    
};



int main(){

    time_t start=time(NULL);
    ROOT::EnableImplicitMT(24);
    cout.precision(3);

    CSTCutG PidCut;
    TString fileName = "/home/sghwang/workspace/SmallHIMACData_v5.6.root";
    CSTFileInput fileIn(fileName);
    TTree* treein = fileIn.getTree();

    int siHitNum;
    double csiData[2][2][2]; //[right,left][down,up][energy,sumADC]
    double siData[2][5][3][2]; //[right,left][hitnum][energy, ADC, position] energy=(omic,junc), ADC=(omic,junc), position(x,y)
    double trackFitData[2][2][2][3]; //[Right, Left][hitnum][xy,yz][slope,const]
    
    treein->SetBranchStatus("*",0);
    treein->SetBranchStatus("siHitNum",1);
    treein->SetBranchStatus("csiData",1);
    treein->SetBranchStatus("siData",1);
    treein->SetBranchStatus("trackFitData",1);
    treein->SetBranchAddress("siHitNum",&siHitNum);
    treein->SetBranchAddress("csiData",&csiData);
    treein->SetBranchAddress("siData",&siData);
    treein->SetBranchAddress("trackFitData",&trackFitData); 


 
    int TotalEventNum = treein->GetEntries();
    int BothEvent =0;
    int EventList[6] = {0};

    vector<int> hasEvent;
    vector<int> hasParticle;
    TCanvas *c1 =new TCanvas("c1","c1",800,800);
    c1->SetMargin(0.15,0.15,0.1,0.1);
    TH1D *ypos = new TH1D("ypos","ypos",   25,-70,70);
    TH1D *yposs = new TH1D("yposs","position_y at Si detector;pos [mm] ; count ",25,-70,70);
    
    int event_nocut = 0 ;
    int event_cut1  = 0 ;
    int event_cut2  = 0 ;
    EventData evtdata;
    std::vector<EventData> event_both;


    for(int event = 0; event < TotalEventNum ; event++){
        // if(event%12346==0)cout <<fixed<< "    Event Cut Progress = " << event<< " / " << TotalEventNum << " ( " << (double)event/(double)TotalEventNum*100.<<" % )           \r" << flush;
        
        treein->GetEntry(event);
        if(siHitNum<1)continue;
        

        bool rightEvent[2] = {0,0};
        bool leftEvent[2]  = {0,0};

        for(int wing=0; wing <2; wing++){
            double CsI_down_energy = csiData[wing][0][0];
            double CsI_up_energy   = csiData[wing][1][0];
            if(CsI_down_energy<0.1 && CsI_up_energy<0.1 )continue;

            for(int hitnum = 0 ; hitnum < siHitNum ; hitnum++){
                double Si_omic_energy = siData[wing][hitnum][0][0];
                double Si_junc_energy = siData[wing][hitnum][0][1];
                double Si_energy      = Si_junc_energy+Si_omic_energy;
                
                double si_Y_position  = siData[wing][hitnum][2][1];
                event_nocut++;
                if(Si_energy < 0.1){continue;}
                event_cut1++;
                if(trackFitData[wing][hitnum][0][0]<-500){continue;}
                event_cut2++;
                
        
                if(wing==1){
                    double pos_y = CSGetTargetZ::getTargetSi(wing,trackFitData[wing][hitnum][0][0],trackFitData[wing][hitnum][0][1]);
                    yposs->Fill(pos_y);
                    ypos->Fill(si_Y_position);
                    leftEvent[hitnum] = 1;
                }
                else if(wing==0){
                    double pos_y = CSGetTargetZ::getTargetSi(wing,trackFitData[wing][hitnum][0][0],trackFitData[wing][hitnum][0][1]);
                    yposs->Fill(pos_y);
                    ypos->Fill(si_Y_position);
                    rightEvent[hitnum] = 1;
                }
            



            }
        
        }
        cout << event_nocut << "       " << event_cut1 << "        " << event_cut2 << "      " << leftEvent[1] << endl;
        for(int rgthit = 0 ; rgthit < 2 ;  rgthit++ ){
            for(int lfthit = 0; lfthit <2 ; lfthit++){
                if(rightEvent[rgthit] && leftEvent[lfthit]){
                    
                    
                    evtdata = {event, rgthit, lfthit};
                    // cout << evtdata.eventnum << "   " << evtdata.rightwinghitnum << "    " << evtdata.lefttwinghitnum << endl;
                    event_both.push_back(evtdata);
                }
            }
        }

    }

    
    cout << "    Event Cut Progress = \033[3mDone\033[0m ( 100.0 % )                                                                        " << endl;
    cout << "    Cut Event Info : " << event_nocut << "   " << event_cut1 << "   " << event_cut2  << endl;
    cout << "    Cutted Both Event : " << event_both.size() << endl;
    cout << "    p-p : " << EventList[0] << "    p-d : " << EventList[1] << "    p-t : " <<EventList[2] << "    d-d : " << EventList[3] << "    d-t : " << EventList[4] << "    t-t : " << EventList[5] << endl;
    time_t end = time(NULL);
    cout << "    Elapsed Time =  " << end-start << " Seconds"<<endl;
    cout << "---------------------------------------------------------------------------------------------------"<<endl;
    gStyle->SetOptStat(0);
    ypos->Draw();
    yposs->Draw("same");
    
    
    
    ypos->SetLineColor(kRed);
   
    // c1->SaveAs("../pictures/ypos.png");

    return 0;
}