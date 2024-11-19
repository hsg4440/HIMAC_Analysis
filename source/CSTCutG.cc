#include "CSTCutG.hh"

CSTCutG::CSTCutG(){
    PCut = new TCutG("pCut");
    PCut -> SetPoint(0,1,4.8);
    PCut -> SetPoint(1,20,4.25);
    PCut -> SetPoint(2,40,3.5);
    PCut -> SetPoint(3,60,3);
    PCut -> SetPoint(4,80,2.7);
    PCut -> SetPoint(5,100,2.4);
    PCut -> SetPoint(6,115,2.25);
    PCut -> SetPoint(7,135,2.13);
    PCut -> SetPoint(8,145,2.7);
    PCut -> SetPoint(9,135,3.4);
    PCut -> SetPoint(10,115,3.75);
    PCut -> SetPoint(11,100,4);
    PCut -> SetPoint(12,80,4.3);
    PCut -> SetPoint(13,60,4.7);
    PCut -> SetPoint(14,40,5.5);
    PCut -> SetPoint(15,20,6.2);
    PCut -> SetPoint(16,1,7);
    DCut = new TCutG("dCut");
    DCut -> SetPoint(0,1,7.5);
    DCut -> SetPoint(1,20,6.7);
    DCut -> SetPoint(2,40,5.9);
    DCut -> SetPoint(3,60,5.15);
    DCut -> SetPoint(4,80,4.6);
    DCut -> SetPoint(5,100,4.15);
    DCut -> SetPoint(6,115,3.85);
    DCut -> SetPoint(7,135,3.4);
    DCut -> SetPoint(8,145,3.3);
    DCut -> SetPoint(9,155,3.4);
    DCut -> SetPoint(10,145,4.2);
    DCut -> SetPoint(11,130,4.65);
    DCut -> SetPoint(12,110,5.15);
    DCut -> SetPoint(13,100,5.45);
    DCut -> SetPoint(14,80,6);
    DCut -> SetPoint(15,60,6.75);
    DCut -> SetPoint(16,40,7.5);
    DCut -> SetPoint(17,20,8.55);
    DCut -> SetPoint(18,1,9.5);
    TCut = new TCutG("tCut");
    TCut -> SetPoint(0,1,9.7);
    TCut -> SetPoint(1,20,8.7);
    TCut -> SetPoint(2,40,7.7);
    TCut -> SetPoint(3,60,6.9);
    TCut -> SetPoint(4,80,6.2);
    TCut -> SetPoint(5,100,5.75);
    TCut -> SetPoint(6,115,5.25);
    TCut -> SetPoint(7,135,4.75);
    TCut -> SetPoint(8,145,4.5);
    TCut -> SetPoint(9,150,4.6);
    TCut -> SetPoint(10,145,5.6);
    TCut -> SetPoint(11,130,6.5);
    TCut -> SetPoint(12,110,7);
    TCut -> SetPoint(13,100,7.25);
    TCut -> SetPoint(14,80,8);
    TCut -> SetPoint(15,60,8.8);
    TCut -> SetPoint(16,40,9.6);
    TCut -> SetPoint(17,20,10.75);
    TCut -> SetPoint(18,1,12);
    wPCut = new TCutG("wpCut");
    wPCut -> SetPoint(0,1,3.6);
    wPCut -> SetPoint(1,20,3);
    wPCut -> SetPoint(2,40,2.5);
    wPCut -> SetPoint(3,60,2.2);
    wPCut -> SetPoint(4,80,2);
    wPCut -> SetPoint(5,100,1.8);
    wPCut -> SetPoint(6,115,1.6);
    wPCut -> SetPoint(7,130,1.6);
    wPCut -> SetPoint(8,140,2);
    wPCut -> SetPoint(9,130,2.6);
    wPCut -> SetPoint(10,115,2.8);
    wPCut -> SetPoint(11,100,2.9);
    wPCut -> SetPoint(12,80,3.2);
    wPCut -> SetPoint(13,60,3.7);
    wPCut -> SetPoint(14,40,4.2);
    wPCut -> SetPoint(15,20,4.8);
    wPCut -> SetPoint(16,1,5.4);
    wDCut = new TCutG("wdCut");
    wDCut -> SetPoint(0,1,5.6);
    wDCut -> SetPoint(1,20,5.1);
    wDCut -> SetPoint(2,40,4.5);
    wDCut -> SetPoint(3,60,4.0);
    wDCut -> SetPoint(4,80,3.3);
    wDCut -> SetPoint(5,100,3.0);
    wDCut -> SetPoint(6,115,2.9);
    wDCut -> SetPoint(7,130,2.6);
    wDCut -> SetPoint(8,140,2.7);
    wDCut -> SetPoint(9,150,2.9);
    wDCut -> SetPoint(10,140,3.2);
    wDCut -> SetPoint(11,130,3.5);
    wDCut -> SetPoint(12,115,3.8);
    wDCut -> SetPoint(13,100,4.1);
    wDCut -> SetPoint(14,80,4.6);
    wDCut -> SetPoint(15,60,5.2);
    wDCut -> SetPoint(16,40,5.8);
    wDCut -> SetPoint(17,20,6.4);
    wDCut -> SetPoint(18,1,7.2);
    wTCut = new TCutG("wtCut");
    wTCut ->SetPoint(0,1,7.3);
    wTCut ->SetPoint(1,20,6.5);
    wTCut ->SetPoint(2,40,5.9);
    wTCut ->SetPoint(3,60,5.3);
    wTCut ->SetPoint(4,80,4.7);
    wTCut ->SetPoint(5,100,4.3);
    wTCut ->SetPoint(6,115,3.9);
    wTCut ->SetPoint(7,130,3.6);
    wTCut ->SetPoint(8,140,3.4);
    wTCut ->SetPoint(9,150,4);
    wTCut ->SetPoint(10,140,4.6);
    wTCut ->SetPoint(11,130,5);
    wTCut ->SetPoint(12,115,5.3);
    wTCut ->SetPoint(13,100,5.7);
    wTCut ->SetPoint(14,80,6.1);
    wTCut ->SetPoint(15,60,6.7);
    wTCut ->SetPoint(16,40,7.3);
    wTCut ->SetPoint(17,20,8.3);
    wTCut ->SetPoint(18,1,9.2);
}
CSTCutG::~CSTCutG(){
    delete PCut;
    delete DCut;
    delete TCut;
    delete wPCut;
    delete wDCut;
    delete wTCut;
}

int CSTCutG::IsIn(double CsIEnergy, double SiEnergy){
    if(PCut->IsInside(CsIEnergy, SiEnergy))return 1;
    else if(DCut->IsInside(CsIEnergy, SiEnergy))return 2;
    else if(TCut->IsInside(CsIEnergy, SiEnergy))return 3;
    else return 0;
}
int CSTCutG::IsIn(bool isWeird, double CsIEnergy, double SiEnergy){
    if(isWeird){
        if(wPCut->IsInside(CsIEnergy, SiEnergy))return 1;
        else if(wDCut->IsInside(CsIEnergy, SiEnergy))return 2;
        else if(wTCut->IsInside(CsIEnergy, SiEnergy))return 3;
        else return 0;
    }
    else{
        if(PCut->IsInside(CsIEnergy, SiEnergy))return 1;
        else if(DCut->IsInside(CsIEnergy, SiEnergy))return 2;
        else if(TCut->IsInside(CsIEnergy, SiEnergy))return 3;
        else return 0;
    }
}



