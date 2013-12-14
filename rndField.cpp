#include <iostream>
#include "math.h"

#include "healpix_base.h"
#include "healpix_map.h"
#include "planck_rng.h"
#include "pointing.h"
#include "error_handling.h"

#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TH2F.h"


using namespace std;

class rndField {
    public:
        rndField(int);
        double *getRaDec();
        void genField(int, int);
        double* getMap();
    private:
        Healpix_Base _base;
        int _nside;
        double* _map;

};

rndField::rndField(int nside) : _nside(nside) {
    
    _base = Healpix_Base(_nside,RING);
}

void rndField::genField(int cell, int values) {

    _map = new double[2*values];

    pointing p = _base.pix2ang(cell);
    double pi = acos(-1.0);
    double th_c = p.theta;
    double ph_c = p.phi;

    double cth_min = cos(th_c-7.0*pi/180.0);
    double cth_max = cos(th_c+7.0*pi/180.0);
    if (cth_min > cth_max) {
        int tmp = cth_max;
        cth_max = cth_min;
        cth_min = tmp;
    }
    double ph_min = ph_c-7.0*pi/180.0;
    double ph_max = ph_c+7.0*pi/180.00;

    planck_rng rng;

    int count = 0;
    while (count < values) {
       double phi = (ph_max-ph_min)*rng.rand_uni()+ph_min;
       double theta = acos( (cth_max-cth_min)*rng.rand_uni()+cth_min );

        pointing pt(theta,phi);
        int pix = _base.ang2pix(pt);
        if (pix != cell) {
            continue;
        } 
        _map[count] = phi;
        _map[values+count+1] = pi/2.0 - theta;
        count++;
        cout << count << " " << theta*180.0/pi << " - " << phi*180.0/pi << " "  << pix <<  endl;
    }
}

double* rndField::getMap() {
    return _map;
}

int main(int argc, char* argv[]) {
    
    double pi = acos(-1.0);
    planck_rng rng;
//    Healpix_Map<double> map = Healpix_Map<double>();

    int nsides = 3;
    rndField* field = new rndField(nsides);


    Healpix_Base base = Healpix_Base(3,RING);

    double cth_min = cos(53.0*pi/180.0);
    double cth_max = cos(66.0*pi/180.0);
    double ph_min = 156.0*pi/180.0;
    double ph_max = 170.0*pi/180.00;

    int count = 0;

    TApplication* rootapp = new TApplication("Example",&argc,argv);
    TH2* hdd = new TH2F("hdd", "Random galaxy phi / theta distribution", 100, 53, 66, 100, 156, 170);

    int cell = 190;
    int nvalues = 10;
    field->genField(cell,nvalues); 

    double* map = field->getMap();
    for (int i=0; i<nvalues; i++) {
        cout << i << " " << map[i]*180.0/pi << " " << map[nvalues+i+1]*180.0/pi << " " << endl;
    } 

    return 0;
}


