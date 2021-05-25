# Exo-FMS_column_semigrey

Elspeth KH Lee - May 2021

This is one part of a series of codes that build upon different two-stream approaches and schemes, primarily useful for the GCM modelling community.
This is the semi-grey version, where one band is in the visible, representing incident radiation from a host star, and one band in the IR, representing the internal radiation propigating inside the planetary atmosphere.

Some useful references for semi-grey RT modelling are: \
Guillot (2010) \
Heng et al. (2011) \
Rauscher & Menou (2012) \
Parmentier et al. (2014, 2015) \
Lee et al. (2021)

The advantage of semi-grey scheme is that they are fast and produce relatively realistic T-p profiles if parameters are chosen carefully.
Larger errors typically occur at low pressure, where the cooling efficiencies are generally too low in the semi-grey framework, leading to highly isothermal upper atmospheres.

To compile enter 'make' in the main directory. To remove compiled code enter 'make clean'.

This code performs various two-stream approaches (non-scattering) from the literature in a semi-grey context:
1. Isothermal layer approximation
2. Toon et al. method
3. Short Characteristics method
4. Heng et al. method
5. Neil Lewis' scattering code, following Pierrehumbert (2010)
6. Mendonca et al. method (IN DEVELOPMENT)

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques
as they would perform inside a GCM setting. This is also useful to see differences in each method and their various approximations.

We also include dry convective adjustment schemes, currently only 'Ray_adj', based on Raymond Pierrehumbert's python code.

# Namelist options

In the file 'FMS_RC.nml' you can select different options that control the simulation, below are the brief descriptions and options for the paramaters.

ts_scheme: \
'Isothermal' - Isothermal ts method \
'Isothermal_2' - Isothermal ts method - high optical depth version \
'Toon' - Toon et al. ts method \
'Shortchar' -  Short characteristics method \
'Heng' - Heng et al. method \
'Lewis' - Neil Lewis method following Pierrehumbert (2010) \
'Lewis_sw' - 'Lewis' but only shortwave scattering (Shortchar IR) \
'Mendonca' - Mendonca et al. method (IN DEVELOPMENT)

opac_scheme: \
'Constant' - constant k_V and k_IR values \
'Heng' - Heng et al. method with fl parameter \
'TK19' - Tan & Komacek (2019) pressure dependent k_V and k_IR UHJ scheme \
'Freedman' - Uses the Rosseland mean fitting function from Freedman et al. (2014) \
'Valencia' - Uses the Rosseland mean fitting function from Valencia et al. (2013)

adj_scheme: \
'Ray_dry' - Ray Pierrehumbert's dry convective adjustment scheme

--- The option 'None' for each of these scheme will it off (e.g. To run without conv adjustment set adj_scheme = 'None') ---

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
fl - The Heng et al. (2011) parameter used for pressure dependent IR optical depths

sw_ac - shortwave single scattering albedo (constant at all layers) \
sw_gc - shortwave asymmetry factor (constant at all layers) \
lw_ac - longwave single scattering albedo (constant at all layers) \
lw_gc - longwave asymmetry factor (constant at all layers)

iIC - Initial condition selection integer (3 = Guillot (2010) profile using k_V and k_IR) \
corr - Flag to perform the adiabatic gradient correction in the initial conditions \
met - metallicty in dex solar (M/H)

# Plotting results

A python plotting routine, 'plot_TP.py' is provided to plot the results of the simulation.

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in intensity (to calculate the flux), various values and weights can be changed for testing at the top of the two-stream modules.
You will need to clean and recompile the code if these are changed.

# Personal recommendations

We generally recommend that the short characteristics method be used, as it is fast, efficient, very stable and also very accurate. This is currently what is used inside Exo-FMS for the Hot Jupiter simulations, and is even fast enough for high-resolution cases.

# Future developments

We will include versions that include multiple-scattering in the future. \
Additional opacity schemes from the literature can be added. \
Ability to include solid surface temperatures and temperature evolution, this involves some extra switches and boundary conditions. \
Window fractions and functions.
