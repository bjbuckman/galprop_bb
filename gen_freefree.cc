
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_freefree.cc *                          galprop package * 3/08/2017
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate free free emissivity and absorption
/*
emissivity (erg cm^-3 sr^-1 Hz^-1 s^-1)
*/
using namespace std;

#include"galprop_classes.h"
#include"galprop_internal.h"

#include "kappa_free_free.h"
#include <cassert>
#include <string>
#include <cstring>
#include <ErrorLogger.h>


int Galprop::gen_freefree()
{
    cout<<"gen_freefree_emiss"<<endl;
    cout<<"generating freefree    emissivity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

    int stat=0;
    double Te=galdef.HII_Te;
    double clumping_factor=galdef.HII_clumping_factor;//AWS20110704
    double emiss_free_free;
    int debug=0;

    for(int inu=0; inu<galaxy.n_nu_synchgrid; inu++) {

        ////////////////////////////////////////// 2D

        if(2 == galaxy.n_spatial_dimensions) {
            for(int ir=0; ir < galaxy.n_rgrid; ir++) {
                for(int iz=0; iz < galaxy.n_zgrid; iz++) {
                    emiss_free_free = 0;
                    galaxy.free_free_absorb.d2[ir][iz].s[inu] = kappa_free_free(galaxy.nu_synch[inu], galaxy.n_HII.d2[ir][iz].s[0], Te, emiss_free_free, 0, debug)*clumping_factor;
                    // galaxy.free_free_absorb.d2[ir][iz].s[inu]*= clumping_factor;
                    // galaxy.free_free_absorb.d2[ir][iz].s[inu] = kappa_free_free_2D( galaxy.nu_synch[inu], ir, iz, Te, clumping_factor, emiss_free_free, 0, debug);
                    galaxy.free_free_emiss.d2[ir][iz].s[inu] = emiss_free_free*clumping_factor;
                    // galaxy.free_free_emiss.d2[ir][iz].s[inu]*= clumping_factor;
                }
            }
        }

        ////////////////////////////////////////// 3D

        if(3 == galaxy.n_spatial_dimensions) {
            for(int ix = 0; ix < galaxy.n_xgrid; ++ix) {
                for(int iy = 0; iy < galaxy.n_ygrid; ++iy) {
                    for(int iz = 0; iz < galaxy.n_zgrid; ++iz) {
                        emiss_free_free = 0;
                        galaxy.free_free_absorb.d3[ix][iy][iz].s[inu] = kappa_free_free(galaxy.nu_synch[inu], galaxy.n_HII.d3[ix][iy][iz].s[0], Te, emiss_free_free, 0, debug)*clumping_factor;
                        // galaxy.free_free_absorb.d3[ix][iy][iz].s[inu]*= clumping_factor;
                        // galaxy.free_free_absorb.d3[ix][iy][iz].s[inu] = kappa_free_free(galaxy.nu_synch[inu], galaxy.n_HII[ix][iy][iz].s[0], Te, emiss_free_free, 0, debug);
                        galaxy.free_free_emiss.d3[ix][iy][iz].s[inu] = emiss_free_free*clumping_factor;
                        // galaxy.free_free_emiss.d3[ix][iy][iz].s[inu]*= clumping_factor;
                    }
                }
            }
        }
    }

    if(galdef.verbose==-801) { // AWS20050727
        cout<<"   freefree    emissivity "<<endl;
        galaxy.free_free_emiss.print();
        cout<<"   freefree    absorption "<<endl;
        galaxy.free_free_absorb.print();
    }

    cout<<" <<<< gen_free_free_emiss"<<endl;
    return stat;
}
