#ifndef CSTFILEINPUT_HH
#define CSTFILEINPUT_HH

#include <iostream>
#include <string>


#include "TTree.h"
#include "TFile.h"
#include "TError.h"
using namespace std;
class CSTFileInput {
    public:
     CSTFileInput(const TString& fileName);
     ~CSTFileInput();
     
     TTree* getTree() const;
    
    private:
     TFile* fileIn;
     TTree* treeIn;
    
    
};

#endif