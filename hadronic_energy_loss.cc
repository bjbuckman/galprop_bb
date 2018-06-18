using namespace std;
#include <cmath>
#include <iostream>
#include "constants.h"

// hadronic energy losses for nuclei
// AWS20150412


/*
A          nucleus mass number
Etot       total energy of nucleus (for consistency with nucleon_loss routine in energy_loss.cc) in MeV
nhcm3      hydrogen gas density cm-3, HI, H_2 and HII
he_to_h    He-to-H ratio by number
method     choice of method
debug      0=silent, 1=print debug information
result     dE/dt in  eV s-1 for nucleus  (for consistency with nucleon_loss routine in energy_loss.cc)
*/

///////////////////////////////////////////////

double hadronic_energy_loss(int A, double Etot, double nhcm3, double he_to_h,int method, int debug)

{
    double result = 0.;

    double Ekin = Etot - A * m_proton; // kinetic energy of nucleus

    double Ekin_threshold= 282.; //threshold for pion production, MeV/nucleus, roughly

    if(Ekin<Ekin_threshold) return result;

    // simplest possible approximation for check

    if(method==0) {

        // 44 mb A^0.7 from ?
        // 30 mb from  Mannheim&Schlickeiser A&A 286, 983 (1994), A power from  Krakau and Schlickeiser, ApJ 802, 114 (2015)

        double total_inelastic_cross_section = 30.e-27; // cm^2

        double inelasticity = 0.5; //  Mannheim&Schlickeiser A&A 286, 983 (1994), near equation (4.7)

        // dE/dt in eV s^-1  units :  E * cm^2 * cm^-3 * cm s^-1 * 1e6
        // ignoring velocity factor for comparison with Krakau and Schlickeiser (2015)!

        result = Ekin * total_inelastic_cross_section * inelasticity * nhcm3 * C * pow(A, -0.75) * 1e6; // eV s-1
    }

// methods 1-3:

// Based on Krakau and Schlickeiser, ApJ 802, 114 (2015), claimed validity Etot= 1 to 1e8 GeV
// equation (34) for protons on ISM composition (quoted factor= 1.26)
// as check: equation (8) follows from equation (7) using factors 1.26 and 3 from equation (3).
// The formulae assume v=c so not valid near threshold! although v is included in the original equation (3).
// The same comment applies to Mannheim&Schlickeiser A&A 286, 983 (1994)
// Projectile A-dependence from equation (3) : limited validity
// equation (3),  needs to be checked, since written in terms of Lorentz factor gamma, with A power -0.47. Putting gamma=E/(Am) yields -0.75 power.
// NB Ep is total energy in these papers
// Theshold given as 1220 MeV, i.e. Ekin=1220-938=282 MeV.


    if(method==1 && A==1) {
        // only use if A=1 i.e. only for protons
        result = 3.85e-7 * pow(Etot/1e3, 1.28) * pow(Etot/1e3 + 200., -0.2) * nhcm3 ;
    }

    if(method==2 && A<=4) {
        //  only use for A=1-4 i.e.only for protons, deuterium  and 3He, 4He
        result = 3.85e-7 * pow(Etot/1e3, 1.28) * pow(Etot/1e3 + 200., -0.2) * nhcm3 ;
        result*= pow(A, -0.75);
    }


    if(method==3) {
        // use for all A (caution!)
        result = 3.85e-7 * pow(Etot/1e3, 1.28) * pow(Etot/1e3 + 200., -0.2) * nhcm3 ; // equation (34) is in GeV, constant=3.85e-16, convert input from MeV, output to eV.
        result*= pow(A, -0.75);
    }


    if(method==10) {
        // same as method 3 but multiply by 10, for test purposes to see the effects clearly
        result = 3.85e-7 * pow(Etot/1e3, 1.28) * pow(Etot/1e3 + 200., -0.2) * nhcm3 ;
        result*= pow(A, -0.75);
        result*=10.;
    }

    if(debug==1) cout<<"hadronic energy loss: A="<<A<<" Etot="<<Etot<<" nhcm3="<<nhcm3<<" he_to_h="<<he_to_h<<" method="<<method<<" dE/dt="<<result<<" eV s-1"
                         <<" loss time="<<Ekin*1e6/result<<" s "<<Ekin*1e6/result/3.15e7<<" y"<<endl;

    return result;
}

/////////////////////////////////////////////////////// test program
/*
int main()
{

  int A=1;
  double Ekin=1e3;
  double nhcm3=1;
  double he_to_h=0.1;
  int method;
  int debug=1;

  for (A=1;A<=26;A++)
  for (Ekin=100;Ekin<1e7;Ekin*=2)
  {
   for(method=0;method<=3;method++)
  {
    double Etot= Ekin+A*m_proton;
    double dEdt= hadronic_energy_loss(A, Etot, nhcm3, he_to_h, method,debug);

    cout<<"method="<<method<<" A="<<A<<" Ekin="<<Ekin<<" Etot="<<Etot<<" dEdt="<<dEdt<<" eV s^-1"<<" loss time="<<Ekin*1e6/dEdt<<" s "<<Ekin*1e6/dEdt/3.15e7<<" y"     <<endl;
   }
   cout<<endl;
  }


}
*/
