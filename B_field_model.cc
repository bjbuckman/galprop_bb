
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * B_field_model.cc *                            galprop package * 4/14/2000
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<iostream>
#include <cmath>
#include <string>    //AWS20050624

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Magnetic field in Tesla (=1e4 Gauss); r,z in kpc

double B_field_model(double r,double z,int model)
{
    float Bo, rscale, zscale;
    float b_field=0.0;
    float ro=8.5;      // Sun galactocentric distance

    if (model==1) {
        Bo=6.0e-10;
        rscale=20.;
        zscale=5.;
        b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);
    }

// exactly the model used in cal3prop:
    if (model==2) {
        Bo=6.0e-10;
        rscale=20.;
        zscale=5.;
        b_field=Bo *exp(-r/rscale) * exp(-fabs(z)/zscale);
    }

    if (model==3) {
        Bo=5.0e-10;
        rscale=10.0;
        zscale=2.0;
        b_field=Bo*exp(-(r-ro)/rscale)*exp(-fabs(z)/zscale);
    }

    if ((model>=5) && (model<=10)) {
        double BoMW = 5.0e-10; //T
        double rscaleMW = 10.0; //kpc
        double zscaleMW = 2.0; // kpc
        double b_fieldMW = 0.0; // T

        double Bo82 = 5.0e-9; // T
        double rscale82 = 0.2; //kpc
        double zscale82 = 0.05; //kpc
        double bz_power = -2.0;
        double b_field82 = 0.0; // T

        double tix = (model-5.0)/5.0;

        b_fieldMW = BoMW*exp(-(r-ro)/rscaleMW)*exp(-fabs(z)/zscaleMW);
        b_field82 = Bo82*exp(-(r)/rscale82)*pow((fabs(z)+zscale82)/zscale82, bz_power);
        b_field = b_fieldMW*(1.-tix) + b_field82*tix;
    }

    if ((model>=20) && (model<=24)) { //M82 Simple
        int tix = model-20;
        double Bo=1.0e-9; // Telsa
        if (tix == 0) Bo = 1.0e-10;
        if (tix == 1) Bo = 5.0e-10;
        if (tix == 2) Bo = 1.0e-9;
        if (tix == 3) Bo = 5.0e-9;
        if (tix == 4) Bo = 1.0e-8;
        double rscale=0.2; //kpc
        double zscale=0.05; //kpc
        double bz_power=-2.;
        if (r <= rscale && fabs(z) <= zscale) {
            b_field = Bo;
        } else if (r <= rscale && fabs(z) > zscale) {
            b_field = Bo*pow( fabs(z)/zscale, bz_power);
        } else if (fabs(z) <= zscale && r > rscale) {
            b_field = Bo*exp(-(r-rscale)/rscale);
        } else {
            b_field = Bo*exp(-(r-rscale)/rscale)*pow( fabs(z)/zscale, bz_power);
        }

    }

    if ((model>=25) && (model<=32)) { //M82 Simple
        int tix = model-25;
        double Bo=1.0e-9; // Telsa
        double bz_power=-2.;
        if (tix == 0) bz_power = 0.;
        if (tix == 1) bz_power = -0.5;
        if (tix == 2) bz_power = -1.;
        if (tix == 3) bz_power = -2.;
        if (tix == 4) bz_power = -3.;
        if (tix == 5) bz_power = -4.;
        if (tix == 6) bz_power = -5.;
        if (tix == 7) bz_power = -6.;
        double rscale=0.2; //kpc
        double zscale=0.05; //kpc

        if (r <= rscale && fabs(z) <= zscale) {
            b_field = Bo;
        } else if (r <= rscale && fabs(z) > zscale) {
            b_field = Bo*pow( fabs(z)/zscale, bz_power);
        } else if (fabs(z) <= zscale && r > rscale) {
            b_field = Bo*exp(-(r-rscale)/rscale);
        } else {
            b_field = Bo*exp(-(r-rscale)/rscale)*pow( fabs(z)/zscale, bz_power);
        }

    }

    if (model> 100) {
        Bo=model*1.0e-10;
        rscale=20.;
        zscale=5.;
        b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);
    }

// Bo rscale zscale encoded in 9-digit number: BBBrrrzzz in units of 0.1
// e.g. 123456789 : Bo= 12.3E-10 Tesla  rscale=45.6 kpc zscale=78.9 kpc

    if (model > 1000) {
        Bo=           (model/1000000)                 * 0.1 *1.0e-10;
        rscale=(model-(model/1000000)*1000000 )/1000  * 0.1         ;
        zscale=(model%1000)                           * 0.1         ;
        b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);
    }

    //cout<<"B_field (r,z) = ("<<r<<","<<z<<") "<<b_field<<endl;
    return b_field;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// 3D  interface

double B_field_model(double x, double y, double z, int model)
{
    return B_field_model(sqrt(x*x+y*y),z,model);
}
