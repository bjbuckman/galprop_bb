#ifndef GALPROP_H
#define GALPROP_H
#include <string>
#include <vector>
#include "GCR_data.h"
#include "Distribution.h"
#include "Particle.h"
#include "Spectrum.h"
#include "Galaxy.h"
#include "Configure.h"
#include "Galdef.h"
#include "los_integration.h"
#include "DistributionFunction.h"
#include "GasCubeData.h"

using namespace std;

class Galprop
{

public:

    Galprop();
    ~Galprop();

    int Run ( const std::string& galdefPath,
              const std::string& fitsPath,
              const std::string& outputPath,
              const std::string& outputPrefix,
              const std::string& runNumber );//int argc, char** argv);

    int AvXCO ( const std::string& galdefPath,
                const std::string& fitsPath,
                const std::string& outputPath,
                const std::string& outputPrefix,
                const std::string& runNumber );//int argc, char** argv);

    int create_gcr();
    int create_transport_arrays ( Particle& );
    int create_galaxy();
    int create_SNR();
    int cr_luminosity();
    int propagate_particles();
    int e_KN_loss ( Particle& );

    double D_pp ( double,double,double,double );
    int    D_xx ( Particle&,int,int,int,int,int,int );
    static double fu ( double );

// DM routines IMOS20050912
    int gen_DM_source ( Particle& );
    int gen_DM_emiss();
    double DM_profile ( double, double, double );
    double DM_profile_av ( double,double,double,double,double );
    double DM_profile_av ( double,double,double,double,double,double,double );
    int gen_DM_skymap();
    int store_DM_emiss();
    int store_DM_skymap();
    int gen_DM_skymap_pixel ( const double l, const double b, vector<double> &DM );

    double IC_anisotropy_factor ( double,double,double,double,int );//IMOS20060420
    void decayed_cross_sections ( int iz,int ia,int jz,int ja, double *Ekin,int np,double *sigma );

    int nuclei_normalize();
    int electrons_normalize();
    int positrons_normalize(); //AWS20131023

    int propel ( Particle& );
    int propel_diagnostics (); //AWS20110113
    int propel_diagnostics ( Particle& particle,
                             Distribution& alpha1_r,
                             Distribution& alpha1_z,
                             Distribution& alpha1_p,
                             Distribution& alpha2_r,
                             Distribution& alpha2_z,
                             Distribution& alpha2_p,
                             Distribution& alpha3_r,
                             Distribution& alpha3_z,
                             Distribution& alpha3_p,
                             Distribution& total_source_function,
                             double dt );
    int propel_diagnostics ( Particle& particle,
                             Distribution& alpha1_x,
                             Distribution& alpha1_y,
                             Distribution& alpha1_z,
                             Distribution& alpha1_p,
                             Distribution& alpha2_x,
                             Distribution& alpha2_y,
                             Distribution& alpha2_z,
                             Distribution& alpha2_p,
                             Distribution& alpha3_x,
                             Distribution& alpha3_y,
                             Distribution& alpha3_z,
                             Distribution& alpha3_p,
                             Distribution& total_source_function,
                             double dt );

    void protri ( Particle& particle,
                  Distribution& alpha1_x,
                  Distribution& alpha1_y,
                  Distribution& alpha1_z,
                  Distribution& alpha1_p,
                  Distribution& alpha2_x,
                  Distribution& alpha2_y,
                  Distribution& alpha2_z,
                  Distribution& alpha2_p,
                  Distribution& alpha3_x,
                  Distribution& alpha3_y,
                  Distribution& alpha3_z,
                  Distribution& alpha3_p,
                  Distribution& Nx1_,
                  Distribution& Ny1_,
                  Distribution& Nz1_,
                  Distribution& Np1_,
                  Distribution& Nx2_,
                  Distribution& Ny2_,
                  Distribution& Nz2_,
                  Distribution& Np2_,
                  Distribution& Nx3_,
                  Distribution& Ny3_,
                  Distribution& Nz3_,
                  Distribution& Np3_,
                  Distribution& total_source_function,
                  double dt,
                  int nrept_outer,
                  double f_use );

    double source_distribution(const double x, const double y, const double z, int srcModel, const std::vector<double> &parameters);
    int source_SNR_event ( Particle &particle,double t );
    int source_SNR_event_vec ( Particle &particle,double t,
                               float *total_source_function_x );

    int test_Particle();
    int test_source_SNR_event();
    int test_isotope_cs();
    int test_cfactor();

    int print_BC();
    float isrf_energy_density ( float rr, float zz );
    int HIR_iRing ( double RR );

    int read_gcr();
    int read_isrf ( const int version );
//  int read_HIR();
//  int read_COR();

    int read_gas_cubes ( char* ); //ECC20140824
    //double interp_HI_cube(double, double, double);//ECC20140824
    //double interp_CO_cube(double, double, double, double);//ECC20140824
    int get_gas_cube_coords(int, int, int, int, int, int);//ECC20140824

    int read_gas_maps ( char* ); //IMOS20080114
    int gas_iRing ( double ); //IMOS20080114
    double fX_CO ( double ) const;  //IMOS20080114

    //Classes for use with LOSintegrators
    class GasFunction : public SM::LOSfunction<double>
    {
    private:
        GasFunction();
        enum TYPE { HI, H2, CO, HII } ftype;
        const double frInd;
        const Galprop &fgp; //To have access to fX_CO
    public:
        // The functions can be optionally multiplied with radius to the index
        // rPLindex
        GasFunction(const std::string &type, double rPLindex, Galprop& gp);
        virtual double operator () (double x, double y, double z) const;
    };
    class GasEmissFunction : public SM::LOSfunction<std::valarray<double> >
    {
    private:
        GasEmissFunction();
        enum TYPE { BREMSS, PION } ftype;
        const GasFunction fgf;
        DistributionFunction *fdf;
        const Galprop &fgp;
    public:
        GasEmissFunction(const std::string &type, const std::string &gas_type, Galprop& gp);
        ~GasEmissFunction();
        virtual std::valarray<double> operator () (double x, double y, double z) const;
    };

    int store_gcr();
    int store_gcr_full();
    int store_gcr_source_functions ( Particle &particle );//AWS2010031

    void store_IC_emiss ( const std::string& type ); // TAP20090312
    int store_IC_skymap ( const std::string& type );  //AWS20090415
    int store_IC_skymap_comp ( const std::string& type );
    int store_bremss_emiss();
    int store_bremss_ionized_skymap();
    int store_bremss_skymap();
    int store_pi0_decay_emiss();
    int store_pi0_decay_skymap();
    int store_pi0_decay_H2R_skymap(); //AWS20041215
    int store_pi0_decay_HIR_skymap(); //AWS20041215
    int store_pi0_decay_HII_skymap(); //IMOS20080114*
    int store_bremss_H2R_skymap();    //AWS20041215
    int store_bremss_HIR_skymap();    //AWS20041215
    int store_bremss_HII_skymap();    //IMOS20080114*
    int store_synch_emiss();          //AWS20080314
    int store_synch_skymap();
    int store_ionization_rate();
    int store_skymap ( float array[], long naxes[4], const std::string name, double crval[4], double cdelt[4] );
    int store_mapcube_skymap ( float array[],
                               double energy[],
                               const int nComponents,
                               const int nEnergies,
                               const std::string& name,
                               const bool MeV ) const;

    int gen_secondary_source ( Particle& );
    int gen_isrf_energy_density();
    bool gen_IC_emiss(); // TAP20090311
    int gen_IC_skymap();
    //int gen_IC_skymap_pixel(const double l,
    //			    const double b,
    //			    vector< vector<double> >& iso_IC,
    //			    const Distribution& electrons,
    //			    const valarray<double>& electronE,
    //			    const valarray< Skymap<double> >& isrf,
    //			    const valarray<double>& targetE,
    //			    const valarray<double>& cosZeta,
    //			    const valarray<double>& crossSection,
    //			    const unsigned int cosThetaBins,
    //			    vector< vector<double> >& aniso_IC,
    //			    const Skymap<float>& nSkymap,
    //			    const unsigned int nSkymapPixels);

    int gen_IC_skymap_pixel(const double l, const double b,
                            vector< vector<double> >& iso_IC);

    //int gen_IC_skymap_pixel(const double l, const double b,
    //			    const Distribution & electrons,
    // 			    const int ielectrons,
    //			    vector< vector<double> >& iso_IC,
    //			    vector< vector<double> >& aniso_IC,
    //			    const double Etarget[], const double factor);

    //int gen_IC_iso_skymap_pixel(const double l, const double b, vector< vector<double> >& iso_IC);

    int gen_bremss_emiss();
    int gen_bremss_ionized_skymap();
    int gen_bremss_ionized_skymap_pixel ( const double l, const double b, vector<double> &bremss_ionized );
    int gen_bremss_skymap();
    int gen_bremss_skymap_pixel ( const double l, const double b, vector<double> &n_H_av, vector<double> &n_HI_av, vector<double> &n_H2_av, vector<double> &n_CO_av, vector< vector<double> > &emiss_av, vector< vector<double> > &emiss_HII, vector< vector<double> > &emiss_HI, vector< vector<double> > &emiss_H2, const int n_Ring );
    int gen_pi0_decay_emiss();
    int gen_pi0_decay_skymap();
    int gen_pi0_decay_skymap_pixel ( const double l, const double b, vector<double> &n_H_av, vector<double> &n_HI_av, vector<double> &n_H2_av, vector<double> &n_CO_av, vector< vector<double> > &emiss_av, vector< vector<double> > &emiss_HII, vector< vector<double> > &emiss_HI, vector< vector<double> > &emiss_H2, const int n_Ring );
    int gen_synch_emiss();
    int gen_freefree(); //BJB20170308
    int store_freefree(); //BJB20170308

    int gen_synch_skymap();
    int gen_synch_skymap_pixel ( const double l, const double b, vector<double> &synch );
    int gen_synch_skymap_pixel ( const double l, const double b, vector<double> &synch,  vector<double> &Q,  vector<double> &U ); //AWS20100709
    int gen_synch_skymap_pixel  (const double l, const double b, vector<double> &synch, vector<double> &Q,  vector<double> &U, vector<double>&free_free ); //AWS20131004

    int gen_secondary_positron_source ( Particle &particle );
    int gen_tertiary_antiproton_source ( Particle &particle );
    int gen_secondary_antiproton_source ( Particle &particle );
    int gen_secondary_proton_source ( Particle &particle );
    int gen_secondary_deuterium_fusion_source( Particle &particle ); //AWS20140502
    int gen_knock_on_electron_source ( Particle &particle ); //IMOS20060504
    double knock_on_cross_section ( double,double,int );     //IMOS20060504
    int gen_ionization_rate();
    int gen_luminosity(); //AWS20100121
    int test_suite ();

//////////////////////////////////////////////////
//  Hacking away at a comparison routine for cosmic ray data, need these
//  variables.  This should be encapsulated in its own class at one point.
    GCR_data gcr_data;
    double calculate_GCR_chisq ( bool*, bool*, double*, double*, double*, double* ); // The arrays represent what we want to fit and the modulation potential for the plot

//////////////////////////////////////////////////

    Spectrum S;
    Distribution Dist;
    Particle P;
    GasCubeData HICubeData;
    GasCubeData HIICubeData;
    GasCubeData COCubeData;
    GasCubeData COCubeData_rlb;
    GasCubeData HICubeData_rlb;
    GasCubeData HIICubeData_rlb;

    double source_norm_CO;
    double source_norm_param;
    bool spiral_fraction_particle_init;
    Particle spiral_fraction_particle;

    int n_species;
    int isrf_energy_density_i_comp; // required for cfactor

///////////////////////////////////
    gp::Configure configure;
    Galdef galdef;

    Particle* gcr; // all species
    Galaxy galaxy;
};
#endif
