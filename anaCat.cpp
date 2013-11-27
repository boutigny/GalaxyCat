#define Risa_cxx

#include <iostream>
#include <iomanip>

#include "math.h"

#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"

#include "RisaCat.h"

using namespace std;

int main(int argc, char* argv[]) {

// The following is necessary on Ubuntu to avoid the "no dictionary for class ... " warnings
//
    TApplication* rootapp = new TApplication("Example",&argc,argv);

    TFile *f = new TFile("../Catalogs/RW_truth_190.root");
    TTree *tree = (TTree*)f->Get("Risa");
    Risa *r = new Risa(tree);

    if (r->fChain == 0) return 99;

    Long64_t nentries = r->fChain->GetEntriesFast();

    Long64_t nbytes=0, nb=0;
    double pi = acos(-1.0);
    cout << "pi = " << pi << endl;
    double z_1, z_2, ra_1, ra_2, dec_1, dec_2;
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        Long64_t ientry = r->LoadTree(jentry);
        if (ientry < 0) break;

        nb = r->fChain->GetEntry(jentry); 
        z_1 = r->z ; ra_1 = r->ra*pi/180. ; dec_1 = r->dec*pi/180.;
        cout << z_1 << endl;
        for (Long64_t kentry=jentry+1; kentry<nentries; kentry++) {
            Long64_t ikentry = r->LoadTree(kentry);
            if (ikentry < 0) break;
//            if (kentry == jentry) continue;
            nb = r->fChain->GetEntry(kentry);
            z_2 = r->z ; ra_2 = r->ra*pi/180. ; dec_2 = r->dec*pi/180.;

            double theta;  // angular distance
            theta = atan(sqrt(cos(dec_2)*cos(dec_2)*sin(ra_2-ra_1)*sin(ra_2-ra_1) + 
                        (cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))*(cos(dec_1)*sin(dec_2)-sin(dec_1)*cos(dec_2)*cos(ra_2-ra_1))) /
                        (sin(dec_1)*sin(dec_2) + cos(dec_1)*cos(dec_2)*cos(ra_2-ra_1)));

//            cout << "Angular distance : " << theta*180./pi << endl;

            if (kentry > 10000) break;
        }
        if (jentry > 10000) break;
    }
        


    return 0;
}

