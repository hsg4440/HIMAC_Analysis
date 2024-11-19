
const double fPadHeight = 11.9; //[mm]
const double fPadWeith = 1.9; //[mm]
const double fPadGap = 0.1; //[mm]
const double fDistance = 550; //[mm]

const double centerXPad = 32.;
const double centerZPad = 89.;
const double lengthYPad = fPadHeight*3.5 + fPadGap*3.5 + fDistance;
 
struct point
{
    double x;
    double y;
    double z;
};

point test(double x, double y, double z){ // right tpc
    
    const double radian40 = 40.*TMath::Pi()/180.;
    const double cos40 = TMath::Cos(radian40);
    const double sin40 = TMath::Sin(radian40);

    point hit;
    point temp;

    temp.x = x-centerXPad;
    temp.y = y-lengthYPad;
    temp.z = z-centerZPad;
    // 1st Rotate (90 degree x-axis)
    hit.x = temp.x;
    hit.y = temp.z;
    hit.z = -temp.y;

    temp = hit;
    // 2nd Rotate  (90 degree to z-axis)
    hit.x = temp.y;
    hit.y = -1.*temp.x;
    hit.z = temp.z;

    temp = hit;
    // 3th Rotate (40 degree to y-axis)
    hit.x = cos40*temp.x -sin40*temp.z;
    hit.y = temp.y;
    hit.z = sin40*temp.x +cos40*temp.z;

    return hit;
}

void TransformTest()
{
    auto graph = new TGraph2D();
    auto base = new TGraph2D();
    base -> SetPoint(0, -500, -500, -500);
    base -> SetPoint(1, 300, 300, 300);


    for(int i=0; i<3; i++){
        point hits = test(20., i*10.+10., 50.);

        cout << hits.x << " " << hits.y << " " << hits.z << endl;
        graph -> SetPoint(i, hits.x, hits.y, hits.z);
    }

    auto c1 = new TCanvas();

    graph->SetMarkerStyle(20);
    base -> Draw();
    graph -> Draw("pcol,same");

    // auto arr = trans()
}