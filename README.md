# Exo-FMS_column_semigrey

Major Update History
 - May 2021 - initial models
 - Dec 2021 - major overhaul
 - Jun 2022 - Improvements + Bezier short char.
 - Aug 2023 - Added AA and VIM methods, cleanup - add 4SDA method to AA_E, AA_L and VIM
 - Jul 2024 - More methods + updates to methods, split sw and lw method

Elspeth K.H. Lee - Dec 2021

This is one part of a series of codes that build upon different two-stream approaches and schemes, primarily useful for the GCM modeling community.
This is the semi-grey (double-grey) version, where one band is in the visible, representing incident radiation from a host star, and one band in the IR, representing the internal radiation propagating inside the planetary atmosphere.

Some useful references for semi-grey RT modelling are: \
Guillot (2010) \
Heng et al. (2011) \
Rauscher & Menou (2012) \
Parmentier et al. (2014, 2015) \
Lee et al. (2021)

The advantage of semi-grey scheme is that they are fast and produce relatively realistic T-p profiles if parameters are chosen carefully.
Larger deviations to real gas typically occur at low pressure, where the cooling efficiencies are generally too low in the semi-grey framework, leading to highly isothermal upper atmospheres.

To compile enter 'make' in the src directory. To remove compiled code enter 'make clean'. \
Some compiler options for gfortran, nvfortran and intel are provided in the makefile. \
You will need to link a lapack installation, if this is not possible then don't compile the SH2 and SH4 methods and comment out their code in FMS_RC.f90

For the shortwave method we have:
1. Direct beam only method
2. Adding method
3. Spherical harmonic adding method (SDA)
4. Toon89 multiple scattering method
5. Spherical harmonic two-stream multiple scattering method (SH2)
6. Spherical harmonic four-stream multiple scattering method (SH4)
7. Two-stream DISORT version

Adding-Doubling (AD), SH2, SH4, Feautrier, DFE and disort are in development

For the longwave method we have:
1. Absorption approximation with exponential in tau function (AA_E)
2. Absorption approximation with linear in tau function (AA_L)
3. Short characteristics (sc) with linear interpolants 
4. Toon89 multiple scattering method
5. Variational iteration method (VIM)
6. Two-stream DISORT version

Adding-Doubling (AD), SH2, SH4, Feautrier, DFE, sc_para and disort are in development

You can also see the header comments in the source code for some additional information.

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques
as they would perform inside a GCM setting. This is also useful to see differences in each method and their various approximations.

For the longwave radiation, we use the weighted essentially non oscillatory order 4 interpolation method (WENO4) from Janett et al. (2019). This allows for highly smooth interpolation of the temperature from the layers to the levels (typically this is required in GCMs) and strong gradients in temperature are accurately interpolated through. Linear extrapolation is used for the end points to avoid numerical issues and overshoot.

We also include dry convective adjustment schemes, currently 'Ray_adj', based on Raymond Pierrehumbert's python code and mixing length theory (MLT).

Results can be plotted using the plot_TP.py and plot_flux.py python codes in the main directory.

# Auxiliary codes

In the hybrid_sigma directory you can use gen_hybrid_sigma_grid.py to generate a custom hybrid sigma grid for use here or in your GCMs.

# Namelist options

In the file 'FMS_RC.nml' you can select different options that control the simulation, below are the brief descriptions and options for the parameters.

sw_scheme: \
'sw_direct' - direct beam only \ 
'sw_adding' - adding method \
'sw_SDA' - SDA method \
'sw_Toon' - Toon89 method \
'sw_SH2' - SH2 method (experimental) \
'sw_SH4' - SH4 method (experimental) \ 
'sw_disort_ts' - two stream disort method \

lw_scheme: \
'lw_AA_E' - AA_E method \
'lw_AA_L' - AA_L method \
'lw_sc_linear' - sc linear method \
'lw_VIM' - VIM method \
'lw_Toon' - Toon89 method \ 
'lw_disort_ts' - two steam disort method \

opac_scheme: \
'Constant' - constant k_V and k_IR values \
'Heng' - Heng et al. method with fl parameter \
'TK19' - Tan & Komacek (2019) pressure dependent k_V and k_IR UHJ scheme \
'Freedman' - Uses the Rosseland mean fitting function from Freedman et al. (2014) \
'Valencia' - Uses the Rosseland mean fitting function from Valencia et al. (2013)

adj_scheme: \
'Ray_dry' - Ray Pierrehumbert's dry convective adjustment scheme \
'MLT' - Mixing length theory

--- The option 'None' for each of these scheme will turn it off (e.g. To run without conv adjustment set adj_scheme = 'None') ---

nlay, a_sh , b_sh - the number of layers, and filenames that contain the a and b constants for the hybrid sigma grid

pref - reference surface pressure (pa)

t_step - time step in seconds \
nstep - number of integer timesteps \
Rd_air - specific gas constant of the air (J kg-1 K-1)\
cp_air - specific heat capacity (constant pressure) of the air (J kg-1 K-1) \
grav - gravitational acceleration constant (m s-2) \
mu_z - cosine angle of the solar zenith angle \
Tirr - Irradiation temperature \
Tint - Internal temperature

k_V - visible band opacity (m2 kg-1) \
k_IR - IR band opacity (m2 kg-1) \
AB - Bond albedo \
fl - The Heng et al. (2011) parameter used for pressure dependent IR optical depths \
met - metallicity in dex solar (M/H)

ns - Shortwave power index (For Heng opacities) \
nl - Longwave power index (For Heng opacities)

sw_ac - shortwave single scattering albedo (constant at all layers) \
sw_gc - shortwave asymmetry factor (constant at all layers) \
lw_ac - longwave single scattering albedo (constant at all layers) \
lw_gc - longwave asymmetry factor (constant at all layers)

zcorr - include zenith angle correction (.True.) \
zcorr_meth - zenith angle correction method (1,2)  \
radius - radius of the planet at surface (m)

iIC - Initial condition selection integer (3 = Guillot (2010) profile using k_V and k_IR, 6 = adiabatic from Ts_init to Tstrat_init) \
corr - Flag to perform the adiabatic gradient correction in the initial conditions

surf - solid surface is present (.True.) \
Ts_init - initial temperature of surface \
Tstrat_init - initial temperature of stratosphere \
cp_surf - surface areal heat capacity \
sw_a_surf - shortwave surface albedo \
lw_a_surf - longwave surface albedo

# Plotting results

A python plotting routines, 'plot_TP.py'and 'plot_flux.py' are provided to plot the results of the simulation.

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in intensity (to calculate the flux), various values and weights can be changed for testing at the top of the two-stream modules.
You will need to clean and recompile the code if these are changed. Beware, some codes use the classical Gauss-Legendre method (flux = 2pi * intensity), while others include the weighting in the intensity (flux = pi * intensity).

# Personal recommendations

For shortwave scattering problems we recommend the adding method as included, or if more accuracy is needed the SDA or Toon89 methods.
For longwave scattering problems we recommend the linear absorption approximation method or if more accuracy is required the VIM or Toon89 methods.

# Future developments

Ability to include solid surface temperatures and temperature evolution, this involves some extra switches and boundary conditions. \
Additional opacity schemes from the literature can be added. \
Window fractions and functions.
