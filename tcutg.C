void tcutg(){

    double x_pos[10] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.};
    double y_pos[10] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.};

    TH2D* field = new TH2D("temp_field",
                           "Field;x;y",
                           19,
                           -0.5,
                           9.5,
                           19,
                           -0.5,
                           9.5);
    
    TCanvas *c1 = new TCanvas("c1","c1",800,800);

    for(int i = 0; i <10; i++){
        for(int j = 0; j<10; j++){
            field->Fill(i,j,i*j+1);
        }
    }

    TCutG *cut = new TCutG("cut");
    cut->SetPoint(0,2.5,2.5);
    cut->SetPoint(1,6.5,2.5);
    cut->SetPoint(2,6.5,6.5);
    cut->SetPoint(3,2.5,6.5);
    cut->SetPoint(4,2.5,2.5);


//출력 (별 다섯개)
    for(int j = 9 ; j >= 0 ; j--){
        for(int i = 0 ; i < 10; i++){
            cout << cut->IsInside(i,j) << " " ;
        }
        cout << endl;
    }



    field->Draw("colz");
    cut->SetLineColor(kRed);
    cut->SetLineWidth(2);
    cut->Draw("same");
    gStyle->SetOptStat(0);


    c1->SaveAs("field.png");
}