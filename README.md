# Exo-FMS_column_semigrey

Major Update History
 - May 2021 - initial models
 - Dec 2021 - major overhaul
 - Jun 2022 - Improvements + Bezier short char.
 - Aug 2023 - Added AA and VIM methods, cleanup - add 4SDA method to AA_E, AA_L and VIM
 - Jul 2024 - More methods + updates to methods, split sw and lw method
 - Nov 2026 - Major update: adding-doubling (AD), Riccati short characteristics with scattering (UNDEROGOING TESTING), SH2/SH4, VIM_L, DDA, Fu97 and TAS longwave methods, layer-by-layer adding, fixed SDA, SH2/SH4 and AD shortwave methods, n-stream DISORT benchmarks, Exo-FMS_2 Toon, AA and convection schemes, gas giant lower boundary everywhere

Elspeth K.H. Lee - Nov 2026

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

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques as they would perform inside a GCM setting.
This is also useful to see differences in each method and their various approximations.

To compile enter 'make' in the src directory. To remove compiled code enter 'make clean'. \
Some compiler options for gfortran, nvfortran and intel are provided in the makefile. \
No external libraries are needed, the DISORT toolchain (including its own LAPACK routines) is included in src.

# Personal recommendations

For the longwave: \
'lw_Toon89' - the default choice. Fast and robust for any single scattering albedo. \
'lw_DDA' - for strong scattering at Toon89 cost. Setting use_sft = .True. (source function technique along the Hogan (2024) angles) is more accurate than Toon89 but around double cost. \
'lw_sc_scat_ric' - for accurate heating rates. Becomes inaccurate for optically thick, near-conservative scattering (ssa > ~0.95) unless the vertical resolution is increased (NOT PUBLISHED YET). \
'lw_AD' - for reference accuracy, slow method.  \
'lw_AA_L' - if there is no longwave scattering. The fastest method and exact for pure absorption, but inaccurate once scattering is present. \
'lw_TAS' - a fast method for weak scattering a < 0.5 \

For the shortwave: \
'sw_adding' - the default choice. The fastest multiple scattering method, works well for clouds and hazes as well as uniform atmospheres. \
'sw_Toon89' - The well tested Toon89 two-stream method.\
'sw_Liou88' - for four-stream accuracy. The most accurate method apart from sw_AD, at ~2.4x the cost of Toon89 - the fastest four-stream method here. \
'sw_AD' - for reference accuracy, but slow.  \
'sw_direct' - if there is no shortwave scattering.

# Methods

For the shortwave method we have:
1. Direct beam only method
2. Adding method with delta-Eddington layers (Meador & Weaver 1980)
3. Spherical harmonic doubling adding method (SDA, Zhang & Li 2013)
4. Toon89 multiple scattering method
5. Two- and four-stream spherical harmonics (SH2 and SH4), converted from PICASO
6. Adding-doubling (AD) 8-stream multiple scattering, converted from AddingDoublingRT
7. Delta-four-stream discrete ordinates with adding (Liou88, Liou et al. 1988)
8. Two-stream DISORT version
9. n-stream DISORT version (benchmark)

Feautrier shortwave is in development.

For the longwave method we have:
1. Absorption approximation with exponential in tau function (AA_E)
2. Absorption approximation with linear in tau function (AA_L)
3. Short characteristics (sc) with linear interpolants
4. Short characteristics 8-stream discrete ordinates with anisotropic multiple scattering, solved by Riccati invariant imbedding, with a source function technique flux integration along the Hogan (2024) angles (sc_scat_ric)
5. Toon89 multiple scattering method
6. Two- and four-stream spherical harmonics (SH2 and SH4), converted from PICASO
7. Variational iteration method (VIM), exponential and linear in tau (VIM_E and VIM_L)
8. Adding-doubling (AD) 8-stream multiple scattering, converted from AddingDoublingRT
9. Four-stream discrete ordinates adding method (DDA, Zhang et al. 2013), with an optional source function technique flux integration along the Hogan (2024) angles (use_sft, off by default)
10. Delta-two- and four-stream combination method (Fu97, Fu et al. 1997)
11. Chou et al. (1999) scaling with the Tang et al. (2018) adjustment (TAS)
12. Two-stream DISORT version
13. n-stream DISORT version (benchmark)
14. Absorption-only discontinuous finite-element method (DFE)

The direct and iterative short characteristics scattering solvers and the Feautrier method are in development and are disabled. DFE is experimental and does not include scattering.

You can also see the header comments in the source code for some additional information.

For the longwave radiation, we use the weighted essentially non oscillatory order 4 interpolation method (WENO4) from Janett et al. (2019) to interpolate the temperature from the layers to the levels (typically this is required in GCMs).
The interpolation is done in log pressure, which gives a very smooth interpolation and strong gradients in temperature are accurately interpolated through. The end points are linearly extrapolated in log(p)-log(T) space to avoid numerical issues and overshoot.

The model uses a gas giant lower boundary. Shortwave radiation is not reflected at the bottom of the column, and in the longwave the internal temperature sets the net flux at the lower boundary, F_up - F_down = sigma*Tint^4. There is no solid surface.

We also include the dry convection schemes from Exo-FMS_2: the Manabe et al. (1965) adjustment, Raymond Pierrehumbert's adjustment and mixing length theory (MLT).
These act on the column after the radiative heating of the timestep, so the column is convectively stable at the end of each step.

Results can be plotted using the plot_TP.py and plot_flux.py python codes in the main directory.

# Auxiliary codes

In the hybrid_sigma directory you can use gen_hybrid_sigma_grid.py to generate a custom hybrid sigma grid for use here or in your GCMs.

# Namelist options

In the file 'FMS_RC.nml' you can select different options that control the simulation, below are the brief descriptions and options for the parameters.

sw_scheme: \
'sw_direct' - direct beam only \
'sw_adding' - adding method \
'sw_SDA' - four-stream spherical harmonic doubling adding method \
'sw_Toon89' - Toon89 method (Exo-FMS_2 version) \
'sw_SH2' - two-stream spherical harmonic method (PICASO) \
'sw_SH4' - four-stream spherical harmonic method (PICASO) \
'sw_AD' - adding-doubling method \
'sw_Liou88' - delta-four-stream discrete ordinates method (Liou et al. 1988) \
'sw_disort_ts' - two stream disort method \
'sw_disort' - n-stream disort method (benchmark, slow)

lw_scheme: \
'lw_AA_E' - AA_E method \
'lw_AA_L' - AA_L method \
'lw_sc_linear' - sc linear method (same result as AA_L, written as short characteristics) \
'lw_sc_scat_ric' - sc 8-stream scattering method, Riccati solution \
'lw_Toon89' - Toon89 method (Exo-FMS_2 version) \
'lw_SH2' - two-stream spherical harmonic method (PICASO) \
'lw_SH4' - four-stream spherical harmonic method (PICASO) \
'lw_VIM_E' - VIM method, exponential in tau \
'lw_VIM_L' - VIM method, linear in tau \
'lw_AD' - adding-doubling method \
'lw_DDA' - four-stream discrete ordinates adding method \
'lw_Fu97' - delta-two- and four-stream combination method \
'lw_TAS' - Chou et al. (1999) scaling with the Tang et al. (2018) adjustment \
'lw_disort_ts' - two stream disort method \
'lw_disort' - n-stream disort method (benchmark, slow) \
'lw_DFE' - absorption only DFE method (requires lw_ac = 0)

opac_scheme: \
'constant' - constant k_V and k_IR values \
'Heng' - Heng et al. method with fl parameter \
'TK19' - Tan & Komacek (2019) pressure dependent k_V and k_IR UHJ scheme \
'Freedman' - Uses the Rosseland mean fitting function from Freedman et al. (2014) \
'Valencia' - Uses the Rosseland mean fitting function from Valencia et al. (2013)

adj_scheme: \
'Manabe' - Manabe et al. (1965) dry convective adjustment (recommended) \
'Pierrehumbert' - Ray Pierrehumbert's dry convective adjustment ('Ray_dry' also works) \
'MLT' - Mixing length theory, also gives Kzz

--- The option 'none' for each of these scheme will turn it off (e.g. To run without conv adjustment set adj_scheme = 'none') ---

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
radius - planetary radius at the model's lower boundary (m)

iIC - Initial condition selection integer (3 = Guillot (2010) profile using k_V and k_IR, 6 = adiabatic from Tdeep_init to Tstrat_init) \
corr - Flag to perform the adiabatic gradient correction in the initial conditions

Tdeep_init - initial deep atmosphere temperature for the adiabatic initial condition \
Tstrat_init - initial temperature of stratosphere

# n-stream DISORT benchmark

sw_disort and lw_disort use the full DISORT 4.0.99 solver and are meant as accuracy benchmarks rather than for use in a GCM, they are a lot slower than the other methods.
The number of streams is set by the nstr parameter at the top of sw_disort_mod.f90 and lw_disort_mod.f90 (default 32), you will need to clean and recompile if you change this.

DISORT assumes a surface, which a gas giant does not have. Rather than changing DISORT, the internal flux lower boundary is found from two DISORT calls with different bottom temperatures.
Since the fluxes are linear in the bottom emission, the right combination of the two gives exactly F_up - F_down = sigma*Tint^4 at the bottom.

The DISORT toolchain is compiled in double precision (see the makefile), as single precision is not accurate enough for the small internal flux compared to the large fluxes deep in a gas giant atmosphere.

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in intensity (to calculate the flux), various values and weights can be changed for testing at the top of the modules.
Most longwave methods use the 4 node Gauss-Jacobi quadrature from Hogan (2024).
The schemes where the quadrature points are also the scattering streams (sc_scat_ric, AD) use Gauss-Legendre points. sc_scat_ric uses its Gauss-Legendre streams for the scattering only, the fluxes are integrated along the Hogan (2024) angles (use_sft, on by default).
You will need to clean and recompile the code if these are changed. Beware, some codes use the classical Gauss-Legendre method (flux = 2pi * intensity), while others include the weighting in the intensity (flux = pi * intensity).

# Future developments

Additional opacity schemes from the literature can be added. \
Window fractions and functions.
