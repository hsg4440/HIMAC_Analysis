#include<iostream>
#include<sstream>
#include<string>

//ROOT library
#include "TCanvas.h"
#include "TPad.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TH1D.h"

using namespace TMath;
using namespace std;
int particle_pid(bool isWeird, double si_energy, double CsI_energy);
    
void cs_analysis(){
    time_t start = time(NULL);
    //setting multithreading
    ROOT::EnableImplicitMT(24);

    //read file 
    cout << "---------------------------------------------------------------------------------------------------" << endl;
    
    TString fileName = "/home/sghwang/workspace/SmallHIMACData_v5.6.root";
    TFile* filein = new TFile(fileName,"read");
    TTree* treein = (TTree*) filein->Get("data");

    cout << "   Target File Path = " << fileName << endl;

    double csiData[2][2][2]; //[right,left][down,up][energy,sumADC]
    treein->SetBranchAddress("csiData",&csiData);
    double siData[2][5][3][2]; //[right,left][hitnum][energy, ADC, position] energy=(omic,junc), ADC=(omic,junc), position(x,y)
    treein->SetBranchAddress("siData",&siData);
    int siHitNum;
    int event_Counter = 0;
    int event_Counter_left = 0;
    int event_Counter_right  = 0;
    int two_counter = 0;

    int pp_event = 0;
    int pd_event = 0;
    int pt_event = 0;
    int dd_event = 0;
    int dt_event = 0;
    int tt_event = 0;

    int temp_event1 = 0;
    int temp_event2 = 0;
    int temp_count =0;
    treein->SetBranchAddress("siHitNum",&siHitNum);


    //initialize parameters
    TCanvas* c1 = new TCanvas("c1","total",800,800);
    TCanvas* c2 = new TCanvas("c2","4Graph",1600,1600);
    TH2D* Si_CsI = new TH2D("Si_CsI",
                            "Total;CsI Energy[MeV]; Si Energy[MeV]",
                            150,
                            0.1,
                            180,
                            150,
                            0.1,
                            12);
    TH2D* Si_CsI_Part[2][2];
    
    Si_CsI_Part[0][0]= new TH2D("Si_CsI_Right_Down","Right_Down;CsI Energy[MeV]; Si Energy[MeV]",150,1,180,150,0.1,12);
    Si_CsI_Part[0][1]= new TH2D("Si_CsI_Right_Up"  ,"Right_Up;CsI Energy[MeV]; Si Energy[MeV]"  ,150,1,180,150,0.1,12);
    Si_CsI_Part[1][0]= new TH2D("Si_CsI_Left_Down" ,"Left_Down;CsI Energy[MeV]; Si Energy[MeV]" ,150,1,180,150,0.1,12);
    Si_CsI_Part[1][1]= new TH2D("Si_CsI_Left_Up"   ,"Left_Up;CsI Energy[MeV]; Si Energy[MeV]"   ,150,1,180,150,0.1,12);
    int total_eventnum = treein ->GetEntries();
    cout << "   Total Event Number = " << total_eventnum << endl;
    cout.precision(4);

    //containing energy data in "for loop"

    for(int event=0; event<total_eventnum; event++){
        if(event%1337==0)cout<<"   Event Progress = "<<event<<" / "<<total_eventnum << " ( "<< (double)event/(double)total_eventnum*100.<<" % )"<<"         \r"<<flush;
        
        treein->GetEntry(event);
        if(siHitNum!=1)continue;
        if(siData[0][0][0][0]+siData[0][0][0][1] >0 && siData[1][0][0][0]+siData[1][0][0][1] > 0){
            two_counter++;
        }
        temp_event1 = 0;
        temp_event2 = 0;


        //loop for each wing
        for(int wing=0; wing < 2 ; wing ++){
            

            double CsI_down_energy = csiData[wing][0][0];
            double CsI_up_energy = csiData[wing][1][0];
            if(CsI_down_energy<0.1 && CsI_up_energy<0.1)continue;
            //loop for each hitnum
        
            for(int hitnum = 0 ; hitnum < siHitNum ; hitnum++) {
                
                double Si_omic_energy = siData[wing][hitnum][0][0];
                double Si_junc_energy = siData[wing][hitnum][0][1];
                double Si_energy = Si_omic_energy+Si_junc_energy;
                
                if(Si_energy <0.1)continue;
                //to determine which CsI contacted with particle
                double si_Y_position =siData[wing][hitnum][2][1];
                // if(wing==1 && si_Y_position>0)Si_energy*=1.3;
                
                //Filling histogram
                if(si_Y_position < 0.){
                    Si_CsI->Fill(CsI_down_energy,Si_energy);
                    
                    Si_CsI_Part[wing][0]->Fill(CsI_down_energy,Si_energy);

                }
                else{
                    Si_CsI->Fill(CsI_up_energy,Si_energy);
                    Si_CsI_Part[wing][1]->Fill(CsI_up_energy,Si_energy);
                }       
                if(wing ==0 ){
                    event_Counter_right++;
                    if(si_Y_position>0)temp_event1=particle_pid(0,Si_energy,CsI_up_energy);
                    else temp_event1=particle_pid(0,Si_energy,CsI_down_energy);
                }
                else if(wing ==1){
                    event_Counter_left++;
                    if(si_Y_position>0)temp_event2=particle_pid(1,Si_energy,CsI_up_energy);
                    else temp_event2=particle_pid(0,Si_energy,CsI_down_energy);
                }
                event_Counter++;         

                

            }
        }
        if(temp_event1!=0 && temp_event2!=0)temp_count++;
        if(temp_event1==1 && temp_event2==1)pp_event++;
        if((temp_event1==1 && temp_event2==2)||(temp_event1 ==2 && temp_event2==1))pd_event++;
        if((temp_event1==1 && temp_event2==3)||(temp_event1 ==3 && temp_event2==1))pt_event++;
        if(temp_event1==2 && temp_event2==2)dd_event++;
        if((temp_event1==3 && temp_event2==2)||(temp_event1 ==2 && temp_event2==3))dt_event++;
        if(temp_event1==3 && temp_event2==3)tt_event++;

   
    }
    
    



    cout << "   Event Progress = Done ( 100.0 % )                                                                        " << endl;
    time_t end = time(NULL);
    cout << "   Elapsed Time =  " << end-start << " Seconds"<<endl;
    cout << "   Hitted event : " << event_Counter - two_counter << "    Left event : " << event_Counter_left << "   Right event : " << event_Counter_right << endl;
    cout << "   Hitted from both event : " << two_counter << endl;
    cout << endl;
    cout << "   Cutted both event count : " <<temp_count << endl;
    cout << "   p-p : " << pp_event << "    p-d : " << pd_event << "    p-t : " <<pt_event << "   d-d : " << dd_event << "  d-t : " << dt_event << "    t-t : " << tt_event << endl;




    gStyle->SetOptStat(0);
    gErrorIgnoreLevel=3000;
    c1->cd();
    
    //drawing histogram
    Si_CsI->Draw("colz");
    
    // c1->SetLogz();
    

    c1->SaveAs("temp.png");
    c2->Divide(2,2);
    c2->cd(1);
    Si_CsI_Part[1][1]->Draw("colz");
    c2->cd(2);
    Si_CsI_Part[0][1]->Draw("colz");
    c2->cd(3);
    Si_CsI_Part[1][0]->Draw("colz");
    c2->cd(4);
    Si_CsI_Part[0][0]->Draw("colz");
    
    c2->SaveAs("temp2.png");
    
    cout << "---------------------------------------------------------------------------------------------------" << endl;
    
        

} 

int particle_pid(bool isWeird, double si_energy, double CsI_energy){
    TCutG* P_cut = new TCutG("P_Cut",17);
    P_cut -> SetPoint(0,1,4.8);
    P_cut -> SetPoint(1,20,4.25);
    P_cut -> SetPoint(2,40,3.5);
    P_cut -> SetPoint(3,60,3);
    P_cut -> SetPoint(4,80,2.7);
    P_cut -> SetPoint(5,100,2.4);
    P_cut -> SetPoint(6,115,2.25);
    P_cut -> SetPoint(7,135,2.13);
    P_cut -> SetPoint(8,145,2.7);
    P_cut -> SetPoint(9,135,3.4);
    P_cut -> SetPoint(10,115,3.75);
    P_cut -> SetPoint(11,100,4);
    P_cut -> SetPoint(12,80,4.3);
    P_cut -> SetPoint(13,60,4.7);
    P_cut -> SetPoint(14,40,5.5);
    P_cut -> SetPoint(15,20,6.2);
    P_cut -> SetPoint(16,1,7);

    TCutG* D_cut = new TCutG("D_Cut",19);
    D_cut -> SetPoint(0,1,7.5);
    D_cut -> SetPoint(1,20,6.7);
    D_cut -> SetPoint(2,40,5.9);
    D_cut -> SetPoint(3,60,5.15);
    D_cut -> SetPoint(4,80,4.6);
    D_cut -> SetPoint(5,100,4.15);
    D_cut -> SetPoint(6,115,3.85);
    D_cut -> SetPoint(7,135,3.4);
    D_cut -> SetPoint(8,145,3.3);
    D_cut -> SetPoint(9,155,3.4);
    D_cut -> SetPoint(10,145,4.2);
    D_cut -> SetPoint(11,130,4.65);
    D_cut -> SetPoint(12,110,5.15);
    D_cut -> SetPoint(13,100,5.45);
    D_cut -> SetPoint(14,80,6);
    D_cut -> SetPoint(15,60,6.75);
    D_cut -> SetPoint(16,40,7.5);
    D_cut -> SetPoint(17,20,8.55);
    D_cut -> SetPoint(18,1,9.5);

    TCutG* T_cut = new TCutG("T_Cut",19);
    T_cut -> SetPoint(0,1,9.7);
    T_cut -> SetPoint(1,20,8.7);
    T_cut -> SetPoint(2,40,7.7);
    T_cut -> SetPoint(3,60,6.9);
    T_cut -> SetPoint(4,80,6.2);
    T_cut -> SetPoint(5,100,5.75);
    T_cut -> SetPoint(6,115,5.25);
    T_cut -> SetPoint(7,135,4.75);
    T_cut -> SetPoint(8,145,4.5);
    T_cut -> SetPoint(9,150,4.6);
    T_cut -> SetPoint(10,145,5.6);
    T_cut -> SetPoint(11,130,6.5);
    T_cut -> SetPoint(12,110,7);
    T_cut -> SetPoint(13,100,7.25);
    T_cut -> SetPoint(14,80,8);
    T_cut -> SetPoint(15,60,8.8);
    T_cut -> SetPoint(16,40,9.6);
    T_cut -> SetPoint(17,20,10.75);
    T_cut -> SetPoint(18,1,12);

    TCutG* wP_Cut = new TCutG("wP_Cut",17);
    wP_Cut -> SetPoint(0,1,3.6);
    wP_Cut -> SetPoint(1,20,3);
    wP_Cut -> SetPoint(2,40,2.5);
    wP_Cut -> SetPoint(3,60,2.2);
    wP_Cut -> SetPoint(4,80,2);
    wP_Cut -> SetPoint(5,100,1.8);
    wP_Cut -> SetPoint(6,115,1.6);
    wP_Cut -> SetPoint(7,130,1.6);
    wP_Cut -> SetPoint(8,140,2);
    wP_Cut -> SetPoint(9,130,2.6);
    wP_Cut -> SetPoint(10,115,2.8);
    wP_Cut -> SetPoint(11,100,2.9);
    wP_Cut -> SetPoint(12,80,3.2);
    wP_Cut -> SetPoint(13,60,3.7);
    wP_Cut -> SetPoint(14,40,4.2);
    wP_Cut -> SetPoint(15,20,4.8);
    wP_Cut -> SetPoint(16,1,5.4);

    TCutG* wD_Cut = new TCutG("wD_Cut");
    wD_Cut -> SetPoint(0,1,5.6);
    wD_Cut -> SetPoint(1,20,5.1);
    wD_Cut -> SetPoint(2,40,4.5);
    wD_Cut -> SetPoint(3,60,4.0);
    wD_Cut -> SetPoint(4,80,3.3);
    wD_Cut -> SetPoint(5,100,3.0);
    wD_Cut -> SetPoint(6,115,2.9);
    wD_Cut -> SetPoint(7,130,2.6);
    wD_Cut -> SetPoint(8,140,2.7);
    wD_Cut -> SetPoint(9,150,2.9);
    wD_Cut -> SetPoint(10,140,3.2);
    wD_Cut -> SetPoint(11,130,3.5);
    wD_Cut -> SetPoint(12,115,3.8);
    wD_Cut -> SetPoint(13,100,4.1);
    wD_Cut -> SetPoint(14,80,4.6);
    wD_Cut -> SetPoint(15,60,5.2);
    wD_Cut -> SetPoint(16,40,5.8);
    wD_Cut -> SetPoint(17,20,6.4);
    wD_Cut -> SetPoint(18,1,7.2);
    
    TCutG* wT_Cut = new TCutG("wT_Cut");
    wT_Cut ->SetPoint(0,1,7.3);
    wT_Cut ->SetPoint(1,20,6.5);
    wT_Cut ->SetPoint(2,40,5.9);
    wT_Cut ->SetPoint(3,60,5.3);
    wT_Cut ->SetPoint(4,80,4.7);
    wT_Cut ->SetPoint(5,100,4.3);
    wT_Cut ->SetPoint(6,115,3.9);
    wT_Cut ->SetPoint(7,130,3.6);
    wT_Cut ->SetPoint(8,140,3.4);
    wT_Cut ->SetPoint(9,150,4);
    wT_Cut ->SetPoint(10,140,4.6);
    wT_Cut ->SetPoint(11,130,5);
    wT_Cut ->SetPoint(12,115,5.3);
    wT_Cut ->SetPoint(13,100,5.7);
    wT_Cut ->SetPoint(14,80,6.1);
    wT_Cut ->SetPoint(15,60,6.7);
    wT_Cut ->SetPoint(16,40,7.3);
    wT_Cut ->SetPoint(17,20,8.3);
    wT_Cut ->SetPoint(18,1,9.2);

    if(isWeird ==1){
        if(wP_Cut->IsInside(CsI_energy,si_energy))return 1;
        else if(wD_Cut->IsInside(CsI_energy,si_energy))return 2;
        else if(wT_Cut->IsInside(CsI_energy,si_energy))return 3;
        else return 0; 
    }
    else{
        if(P_cut->IsInside(CsI_energy,si_energy))return 1;
        else if(D_cut->IsInside(CsI_energy,si_energy))return 2;
        else if(T_cut->IsInside(CsI_energy,si_energy))return 3;
        else return 0; 
    }
}