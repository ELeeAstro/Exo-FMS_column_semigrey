# Exo-FMS_column_semigrey

This is one part of a series of codes that build upon different two-stream approaches and scheme, primarily useful for the GCM modelling community.
This is the semi-grey version, where one band is in the visible, representing incident radiation from a host star, and one band in the IR, representing the internal radiation propigating inside the planetary atmosphere.

The advantage of semi-grey scheme is that they are fast and produce relativly realistic T-p profiles if paramaters are chosen carefully.
The big differences typically occur at low pressure, where the cooling efficencies are generally too low in the semi-grey framework, leading to highly isothermal upper atmospheres.

This code performs various two-stream approaches (non-scattering) from the literature in a semi-grey context:
1. Isothermal layer approximation
2. Toon et al. method
3. Short Characteristics method
4. Heng et al. method
5. Mendonca et al. method (IN DEVELOPMENT)

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques
as they would perform inside a GCM setting. This is also useful to see differences in each method and their various approximations.

We also include a dry convective adjustment schemes, currently only 'Ray_adj', based on Raymond Pierrehumbert's python code.

# Namelist options

In the file 'FMS_RC.nml' you can select different options that control the simulation

ts_scheme: \
'Isothermal' - Isothermal ts method \
'Toon' - Toon et al. ts method \
'Shortchar' -  Short characteristics method \
'Heng' - Heng et al. method \
'Mendonca' - Mendonca et al. method (IN DEVELOPMENT) 

opac_scheme:
'Constant' - constant k_V and k_IR values \
'Heng' - Heng et al. method with fl paramater \
'TK19' - Tan & Komacek (2019) pressure dependent k_V and k_IR UHJ scheme

adj_scheme
'Ray_dry' - Ray Pierrehumbert's dry convective adjustment scheme

The option 'None' for each of these scheme will it off (e.g. To run without conv adjustment set adj_scheme = 'None')

nlay, a_sh , b_sh - the number of layers, and filenames that contain the a and b constants for the hybrid sigma grid

t_step - time step in seconds \
nstep - number of integer timesteps \
Rd_air - specific gas constant of the air \
cp_air - heat capacity at constant pressure of the air \
grav - gravitional acceleration constant (m s-2) \
mu_z - cosine angle of the solar zenith angle \
Tirr - Irradiation temperature \
Tint - Internal temperature

k_V - visible band opacity (cm2 kg-1) \
k_IR - IR band opacity (cm2 kg-1) \
AB - Bond albedo \
fl - The Heng et al. (2011) paramater used for pressure depdendent IR optical depeths

iIC - Initial condition slection integer (3 = Guillot (2010) profile using k_V and k_IR) \
corr - Flag to perform the adibiatic gradient correction in the initial conditions

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in flux, various values can be changed for testing at the top of the two-stream modules.
You will need to clean and recompile the code if these are changed.

# Personal reccomendations

We generally reccomend that the short characteristics method be used, as it is fast, efficenct, very stable and also very accurate. This is currently what is used inside Exo-FMS for the Hot Jupiter simulations, and is even fast enough for high-resolution cases.
