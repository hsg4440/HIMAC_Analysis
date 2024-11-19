#include "CSTFileInput.hh"

CSTFileInput::CSTFileInput(const TString& fileName){

    gErrorIgnoreLevel =kFatal;
    cout << "---------------------------------------------------------------------------------------------------" << endl;
    cout << "    Target File Path = \033[33m" <<  fileName << "\033[0m"<< endl;
    fileIn = new TFile(fileName,"read");
    
    if(fileIn && !fileIn->IsZombie()){
        treeIn = (TTree*)fileIn->Get("data");
    }
    else{
        treeIn = nullptr;
        cout <<"   \033[1;31m File is not properly loaded. Aborting...\033[0m" <<endl;
        cout << "---------------------------------------------------------------------------------------------------" << endl;
        exit(EXIT_FAILURE);
    }
    if(treeIn){
        cout <<"    Total Event Number : " << treeIn->GetEntries() << endl;
    }
    else{
        cout <<"   \033[1;31m Tree is not properly loaded. Aborting...\033[0m" <<endl;
        cout << "---------------------------------------------------------------------------------------------------" << endl;
        exit(EXIT_FAILURE);
    }



    cout << "---------------------------------------------------------------------------------------------------" << endl;
}
CSTFileInput::~CSTFileInput(){
    if(fileIn){
        fileIn->Close();
        delete fileIn;
    }
}

TTree* CSTFileInput::getTree() const{
    return treeIn;
}

