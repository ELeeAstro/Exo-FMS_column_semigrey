# Exo-FMS_column_semigrey

Major Update History
 - May 2021 - initial models
 - Dec 2021 - major overhaul
 - Jun 2022 - Improvements + Bezier short char.
 - Aug 2023 - Added AA and VIM methods, cleanup - add 4SDA method to AA_E, AA_L and VIM

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

To compile enter 'make' in the main directory. To remove compiled code enter 'make clean'. \
Some compiler options for gfortran, nvfortran and ifort are provided in the makefile.

This code performs various 'two-stream' approaches from the literature in a semi-grey context:
1. Isothermal layer approximation
2. Toon89 method (Scattering and non-scattering versions)
3. Short Characteristics method (w. linear and Bezier interpolants)
4. Heng et al. method
5. Neil Lewis's scattering code, following Pierrehumbert (2010)
6. Mendonca et al. method
7. Two-stream DISORT version (w. Tint modifications by Xianyu Tan)
8. Absorption Approximation (AA), with exponential function (_E) or linear (_L)
9. Variational Iteration Method (VIM), following Zhang et al. (2017)

You can also see the header comments in the source code for some additional information.

For the shortwave fluxes: \
Toon89 uses the Toon89 shortwave scheme (Toon et al. 1989) \
AA_E, AA_L and VIM use the alpha-4SDA (4-stream Spherical Harmonic with Doubling-Adding) method (Zhang & Li 2013) \
The rest use an approximate adding method (Mendonca et al. 2015 + references) 

We detect if any albedo is present in the column, and perform the adding method to calculate the scattered flux, otherwise if there is no albedo only the direct beam is used.

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques
as they would perform inside a GCM setting. This is also useful to see differences in each method and their various approximations.

We also include dry convective adjustment schemes, currently only 'Ray_adj', based on Raymond Pierrehumbert's python code.

# Namelist options

In the file 'FMS_RC.nml' you can select different options that control the simulation, below are the brief descriptions and options for the parameters.

ts_scheme: \
'Isothermal' - Isothermal ts method \
'Isothermal_2' - Isothermal ts method - high optical depth version \
'Toon' - Toon89 ts method without scattering \
'Toon_scatter' - Toon89 ts method with scattering \
'Shortchar_linear' - Short characteristics method  with linear interpolants \
'Shortchar_Bezier' - Short characteristics method with Bezier interpolants \
'Heng' - Heng et al. method \
'Lewis_scatter' - Neil Lewis's scattering code, following Pierrehumbert (2010) \
'Mendonca' - Mendonca et al. method \
'Disort_scatter' - two-stream DISORT version with scattering \
'AA_E' - Absorption Approximation with exponential Planck function \
'AA_L' - Absorption Approximation with linear Planck function \
'VIM' - Variational Iteration Method 

opac_scheme: \
'Constant' - constant k_V and k_IR values \
'Heng' - Heng et al. method with fl parameter \
'TK19' - Tan & Komacek (2019) pressure dependent k_V and k_IR UHJ scheme \
'Freedman' - Uses the Rosseland mean fitting function from Freedman et al. (2014) \
'Valencia' - Uses the Rosseland mean fitting function from Valencia et al. (2013)

adj_scheme: \
'Ray_dry' - Ray Pierrehumbert's dry convective adjustment scheme

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

Bezier - use Bezier interpolation for temperature levels (.True.)

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

A python plotting routine, 'plot_TP.py' is provided to plot the results of the simulation.

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in intensity (to calculate the flux), various values and weights can be changed for testing at the top of the two-stream modules.
You will need to clean and recompile the code if these are changed.

# Personal recommendations

For non-scattering problems, we generally recommend that the AA_E method be used, as it is fast, efficient, very stable and also very accurate. 
For shortwave scattering (and non-scattering) problems we recommend the alpha-4SDA method as included (or using the two-stream DISORT or Toon89 with scattering), alpha-4SDA method is generally fast and accurate (enough).
For longwave scattering problems we recommend the VIM method for general use, as a fast and efficient method.
If that fails we recommend the two-stream Toon or DISORT version, which are slower but may be more accurate in some circumstances.

For general use I suggest looking at the VIM module, since it uses alpha-4SDA (shortwave) and alpha-4VIM (longwave), both decent and fast approximation methods that include scattering.
If speed with approximate scattering is super super important then use AA_E and swap the SDA method for the adding method.

# Future developments

We will include more codes that include IR multiple-scattering in the future. \
Add more approximate sw scattering schemes \
Additional opacity schemes from the literature can be added. \
Window fractions and functions.
