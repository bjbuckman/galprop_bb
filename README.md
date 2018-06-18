# galprop_diff_r2504



The second commit into this repo is the galprop version r2504, and the third commit is my modified version, so you can use the diff tool to see the parts I added (also can search "ECC" in a given file typically):

https://github.com/erccarls/galprop_diff_r2504/commit/9502ee75ba9d3585ccded0c6039d758c6834f0eb


    Note that the variable grid spacing does not work since the transport arrays were not built for this and it was more work than it was worth to finish it.  This means disregard changes in galaxy.cc, galaxy.h

    There are then a handful of additional galdef keywords.

### Source distribtion parameters ###

spiral_fraction:
    (aka f_H2) A float between 0 and 1 which specifies the fraction of CR sources distributed according to our star formation model.  The remaining fraction is distributed according to the usual parameterized source dist specification.  The source calculation routines are in create_transport_arrays.cc since we have to integrate over the full galaxy it was much easier than putting them elsewhere.

kennicutt_index:
        This should be called the Schmidt index which specifies the volumetric power-law dependence of SFR and gas density.

kennicutt_threshold:
        Likewise, this specidies the minimal gas density to initiate star formation (in units of cm^-3).

nuc_g1_inner etc..:
        Parameters to change the CMZ injection spectrum (currently hardcoded as r_2D < .5 kpc).  I can't remember how normalizations are handled so be careful if using this.



### Gas Model parameters ###

        nH2_model/nHI_model: 1: use galprop default gas maps for propagation and for gamma-rays.  2: use the new PEB models for propagation and gamma-rays (but does not change source distribution).  The filenames must be specified below.

        COCube_filename: Provide the input CO gas map in cartesian coordinates (x,y,z).  This gas map should be in the galprop/FITS folder containing galprop data.  This is also used for source distributions with non-zero fh2 (which is the "spiral_fraction") keyword

        HICube_filename: same as above, but for HI.  I can provide these files, and you can also make your own additions to the gas/source model by modifying the fits files.  The spatial info is read from the fits headers.

        COCube_rlb and HICube_rlb: filenames of the input gas map in r,l,b coordinates.  These are used if "nH2_model"=2 "nHI_model"=2 which tells the gamma-ray generation and galactic gas model (for propagation not sources) to use the new 3D models.


### Wind parameters ###

        convection: 1-3 are galprop models.  4: turn on radial 3D winds in GC as defined in our paper (1603.06584).
        v0_conv: specified wind velocity in km/s

### Debugging parameters ###

        renorm_off: if 1, don't renormalize pixels to the gas map.
uniform_emiss:
        if 1 make CR's uniform throughout galaxy when calculating skymaps.  used for debugging gas maps
        single_component: if 1 or 2 only outputs HI or H2 contribution to the pi0 emission.

        There are also a few routines for calculating propagation Green's functions for dark matter injected antiprotons, but I don't remember the status of this.
