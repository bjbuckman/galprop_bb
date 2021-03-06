//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_bremss_emiss.cc *                       galprop package * 4/14/2000
// * Modified to output 3D BJB20160615
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"
#include "fitsio.h"

#include <ErrorLogger.h>

int Galprop::store_bremss_emiss()
{

    INFO("Entry");

    int stat;

    stat=0;
    fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status, ii, jj;
    long  fpixel = 1, naxis = 4, nelements, exposure;
    long naxes[4]  ;

    if(galaxy.n_spatial_dimensions==2) {
        naxes[0]=galaxy.n_rgrid;
        naxes[1]=galaxy.n_zgrid;
        naxes[2]=galaxy.n_E_gammagrid;
        naxes[3]=1;
    } //2D

    if(galaxy.n_spatial_dimensions==3) {
        naxes[0]=galaxy.n_xgrid;
        naxes[1]=galaxy.n_ygrid;
        naxes[2]=galaxy.n_zgrid;
        naxes[3]=galaxy.n_E_gammagrid;
    } //3D

    //cout<<galaxy.n_spatial_dimensions<<" "<<naxes[0]<<" "<<naxes[1]<<"  "<<naxes[2]<<endl;

    nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];

    valarray<float> array(0., nelements);

    //char outfile[100];
    //strcpy(outfile,"!"); /* create new file or overwrite existing one */
    //strcat(outfile,configure.fFITSDataDirectory);
    //strcat(outfile,"bremss_emiss_");
    //strcat(outfile,galdef.galdef_ID);
    //strcat(outfile,".gz");
    const std::string outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "bremss_emiss_" + galdef.galdef_ID + ".gz";
    //cout<<" storing bremss_emiss in file "<<outfile<<endl;

    status = 0;         /* initialize status before calling fitsio routines */
    fits_create_file(&fptr, outfile.c_str(), &status);   /* create new file or overwrite existing one */

    /* Create the primary array image (32-bit float pixels */
    fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);

    /* Write a keyword; must pass the ADDRESS of the value */
    exposure = 1500;
    fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status);

    int i=0;

    if(galaxy.n_spatial_dimensions==2) {
        for (int ip=0; ip<naxes[2]; ip++)
            for (int iz=0; iz<naxes[1]; iz++)
                for (int ir=0; ir<naxes[0]; ir++) {
                    array[i]=galaxy.bremss_emiss.d2[ir][iz].s[ip];
                    array[i]*=(galaxy.n_HI.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0]);
                    array[i]*=pow(galaxy.E_gamma[ip],2);
                    i++;
                } //ir,iz,ip
    } //2D

    if(galaxy.n_spatial_dimensions==3) {
        for (int ip=0; ip<naxes[3]; ip++)
            for (int iz=0; iz<naxes[2]; iz++)
                for (int iy=0; iy<naxes[1]; iy++)
                    for (int ix=0; ix<naxes[0]; ix++) {
                        array[i]=galaxy.bremss_emiss.d3[ix][iy][iz].s[ip];
                        array[i]*=(galaxy.n_HI.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0]);
                        array[i]*=pow(galaxy.E_gamma[ip],2);
                        i++;
                    } //ix,iy,iz
    } //3D

    /* Write the array of floats to the image */
    fits_write_img(fptr, TFLOAT, fpixel, nelements, &array[0], &status);

    // write basic FITS keywords
    double crval1,crval2,crval3,crval4;
    double cdelt1,cdelt2,cdelt3,cdelt4;
    int units=1;

    double hi_fac=galdef.HI_norm;
    double hii_fac=galdef.HII_norm;

    if(galaxy.n_spatial_dimensions==2) {
        crval1=galaxy.r_min;
        crval2=galaxy.z_min;
        crval3=log10(galaxy.E_gamma_min);
        crval4=1;

        cdelt1=galaxy.dr;
        cdelt2=galaxy.dz;
        cdelt3=log10(galaxy.E_gamma_factor);
        cdelt4=1;

        fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of radial dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of Z dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of log10(energy grid/MeV)", &status);
        fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Nothing", &status);

        fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of radial dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of Z dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of log10(energy grid/MeV)", &status);
        fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Nothing", &status);
        fits_update_key_null(fptr, "UNITS ", "MeV^2 cm^-3 sr^-1 s^-1 MeV^-1", &status);
        fits_update_key_null(fptr, "HI_MAP", galdef.HICube_filename, &status);
        fits_update_key(fptr, TDOUBLE, "HI_FAC", &hi_fac,"Factor to multiply HI_map", &status);
        fits_update_key_null(fptr, "HIIMAP", galdef.HIICube_filename, &status);
        fits_update_key(fptr, TDOUBLE, "HIIFAC", &hii_fac,"Factor to multiply HIImap", &status);
    } //2D

    if(galaxy.n_spatial_dimensions==3) {
        crval1=galaxy.x_min;
        crval2=galaxy.y_min;
        crval3=galaxy.z_min;
        crval4=log10(galaxy.E_gamma_min);

        cdelt1=galaxy.dx;
        cdelt2=galaxy.dy;
        cdelt3=galaxy.dz;
        cdelt4=log10(galaxy.E_gamma_factor);

        fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of X dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of Y dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of Z dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of log10(energy grid/MeV)", &status);

        fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of X dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of Y dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of Z dimension (kpc)", &status);
        fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of log10(energy grid/MeV)", &status);
        fits_update_key_null(fptr, "UNITS ", "MeV^2 nH_atom^-1 cm^-3 sr^-1 s^-1 MeV^-1", &status);
        fits_update_key_null(fptr, "HI_MAP", galdef.HICube_filename, &status);
        fits_update_key(fptr, TDOUBLE, "HI_FAC", &hi_fac,"Factor to multiply HI_map", &status);
        fits_update_key_null(fptr, "HIIMAP", galdef.HIICube_filename, &status);
        fits_update_key(fptr, TDOUBLE, "HIIFAC", &hii_fac,"Factor to multiply HIImap", &status);
    } //3D

    /*
    // write keywords describing nuclei
    char keyword[20];
    char comment[40];

    for(int i_nucleus=0;i_nucleus<n_species;i_nucleus++){
    sprintf(keyword,"NUCZ%03d",       i_nucleus+1 ); // e.g. NUCZ012
    sprintf(comment,"Z of nucleus %d",i_nucleus+1 );
    cout<<keyword<<" "<<gcr[i_nucleus].Z<<endl;
    fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].Z,comment, &status);

    sprintf(keyword,"NUCA%03d",       i_nucleus+1 ); // e.g. NUCA012
    sprintf(comment,"A of nucleus %d",i_nucleus+1 );
    cout<<keyword<<" "<<gcr[i_nucleus].A<<endl;
    fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].A,comment, &status);


    }
    */

    fits_close_file(fptr, &status);            /* close the file */

    fits_report_error(stderr, status);  /* print out any error messages */
    //delete[] array;       //Gulli20070810
    //return( status );

    //cout<<" <<<< store_bremss_emiss"<<endl;
    //return stat;

    INFO("Exit");

    return status;

}
