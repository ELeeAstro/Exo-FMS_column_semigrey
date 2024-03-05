import numpy as np
from scipy.linalg import solve_banded

def legP(mu): # Legendre polynomials
    return np.array([1, mu, (3*mu**2 - 1)/2, (5*mu**3 - 3*mu)/2,
        (35*mu**4 - 30*mu**2 + 3)/8, 
        (63*mu**5 - 70*mu**3 + 15*mu)/8, 
        (231*mu**6 - 315*mu**4 + 105*mu**2 - 5)/16 ])

def solve_4_stream_banded(M, B, stream):
    """
    Solve the Spherical Harmonics Problem

    Returns
    -------
    intgrl_new : numpy.ndarray
       Integrated source function for source function technique
    flux : numpy.ndarray
        Upwards lux at bottom of atmosphere
    X : numpy.ndarray
        Coefficients of flux/intensity matrix problem
    """
    #   find constants
    diag = int(3*stream/2 - 1)
    #with objmode(X='float64[:]'):
    X = solve_banded((diag,diag), M, B)

    return X


nlayer = 2
Pubar1 = np.ones(4)
Pubar1[:] = 2.0
w_multi = np.ones(((4, nlayer)))
A = np.ones((4,4,nlayer))
Aint=np.zeros((4,nlayer))
for j in range(4):
    Aint = Aint  + w_multi[j,:]*Pubar1[j]*A[j]

print(Aint[0,:],Aint[1,:],Aint[2,:],Aint[3,:])

quit()

stream = 2
nlayer = 54
nlevel = 55
u0 = 1.0
ubar0 = u0
Pu0 = legP(-u0) # legendre polynomials for -u0
f_deltaM = 0.0
Finc = 56703.744189999998

cosb_og = 0.1
w0 = 0.1

surf_reflect = 0.0

tau = np.loadtxt('tau_e.txt')
dtau = tau[1:nlevel] - tau[0:nlayer]

a = np.zeros((stream, nlayer))
b = np.zeros((stream, nlayer))

w_multi = np.ones(((stream, nlayer)))

for l in range(1,stream):
    w = (2*l+1) * cosb_og**l
    w_multi[l,:] = (w - (2*l+1)*f_deltaM) / (1 - f_deltaM)
for l in range(stream):
    a[l,:] = (2*l + 1) -  w0 * w_multi[l,:]
    print(l,(2*l + 1), w_multi[l,0])
    b[l,:] = ( Finc * (w0 * w_multi[l,:])) * Pu0[l] / (4*np.pi)

b_surface = (0. + surf_reflect*u0*Finc*np.exp(-tau[-1]/u0))
b_top = 0.0

eta = np.zeros((2, nlayer)) # will remain zero for thermal
Del = ((1 / ubar0)**2 - a[0,:]*a[1,:])
eta[0,:] = (b[1,:] /ubar0 - a[1,:]*b[0,:]) / Del
eta[1,:] = (b[0,:] /ubar0 - a[0,:]*b[1,:]) / Del

lam = np.sqrt(a[0,:]*a[1,:])
expo = np.minimum(lam*dtau,35)
exptrm = np.exp(-expo)

#   parameters in matrices
q = lam/a[1]
Q1 = (0.5 + q)*2*np.pi
Q2 = (0.5 - q)*2*np.pi

Q1mn = Q1*exptrm;  Q2mn = Q2*exptrm
Q1pl = Q1/exptrm;  Q2pl = Q2/exptrm

zmn = (0.5*eta[0] - eta[1])*2*np.pi
zpl = (0.5*eta[0] + eta[1])*2*np.pi

expon = np.exp(-tau/ubar0)
zmn_up = zmn * expon[1:] 
zpl_up = zpl * expon[1:] 
zmn_down = zmn * expon[:-1] 
zpl_down = zpl * expon[:-1] 

   #   construct matrices
Mb = np.zeros((5, 2*nlayer))
B = np.zeros((2*nlayer))
nlevel = nlayer+1

F = np.zeros((2*nlevel, 2*nlayer))
G = np.zeros((2*nlevel))

#   first row: BC 1
Mb[2,0] = Q1[0]
Mb[1,1] = Q2[0]
B[0] = b_top - zmn_down[0]
#   last row: BC 4
n = nlayer-1
Mb[3, 2*nlayer-2] = Q2mn[n] - surf_reflect*Q1mn[n]
Mb[2, 2*nlayer-1] = Q1pl[n] - surf_reflect*Q2pl[n]
B[2*nlayer-1] = b_surface - zpl_up[n] + surf_reflect * zmn_up[n]
# remaining rows
Mb[0,3::2] = -Q2[1:]
Mb[1,2::2] = -Q1[1:]
Mb[1,3::2] = -Q1[1:]
Mb[2,1:-1:2] = Q2pl[:-1]
Mb[2,2::2] = -Q2[1:]
Mb[3,:-2:2] = Q1mn[:-1]
Mb[3,1:-1:2] = Q1pl[:-1]
Mb[4,:-2:2] = Q2mn[:-1]
B[1:-1:2] = zmn_down[1:] - zmn_up[:-1]
B[2::2] = zpl_down[1:] - zpl_up[:-1]

for i in range(2*nlayer):
    print(Mb[0,i],Mb[1,i],Mb[2,i],Mb[3,i],Mb[4,i])

X = np.zeros((stream*nlayer))
X[:] = solve_4_stream_banded(Mb[:,:], B[:], stream)

#for i in range(2*nlayer):
  #print(i,X[i])


# flux at bottom of atmosphere
F_bot = np.zeros((2*nlayer))
G_bot = 0.0
F_bot[-2] = Q2mn[-1]
F_bot[-1] = Q1pl[-1]
G_bot = zpl_up[-1]

F[0,0] = Q1[0]
F[0,1] = Q2[0]
F[1,0] = Q2[0]
F[1,1] = Q1[0]

nn = np.arange(2*nlayer)
indcs = nn[::2]
k = 0
for i in indcs:
    F[i+2,i] = Q1mn[k]
    F[i+2,i+1] = Q2pl[k]
    F[i+3,i] = Q2mn[k]
    F[i+3,i+1] = Q1pl[k]
    k = k+1

G[0] = zmn_down[0]
G[1] = zpl_down[0]

G[2::2] = zmn_up
G[3::2] = zpl_up


flux_temp = np.zeros((stream*nlevel))
flux_temp[:] = F[:,:].dot(X[:]) + G[:]
flux_bot = np.sum(F_bot[:]*X[:], axis=0) + G_bot

#print(flux_temp[0::2] + u0*Finc*np.exp(-tau[:]/u0))
#print(flux_temp[1::2])
#print(flux_bot)

print(len(flux_temp),flux_temp[:])
#print(flux_bot)

'''

            a = zeros((stream, nlayer, nwno))
            b = zeros((stream, nlayer, nwno))
            w_single = ones((stream, nlayer, nwno))
            w_multi = ones(((stream, nlayer, nwno)))
            p_single = zeros(cosb_og.shape)



            if (w_single_form==1 or w_multi_form==1): # OTHG:
                for l in range(1,stream):
                    w = (2*l+1) * cosb_og**l
                    if w_single_form==1:
                        w_single[l,:,:] = (w - (2*l+1)*f_deltaM) / (1 - f_deltaM)
                    if w_multi_form==1:
                        w_multi[l,:,:] = (w - (2*l+1)*f_deltaM) / (1 - f_deltaM)

            if ((w_single_form==0) or (w_multi_form==0)): # TTHG
                g_forward = constant_forward*cosb_og
                g_back = constant_back*cosb_og
                f = frac_a + frac_b*g_back**frac_c
                f_deltaM_ = f_deltaM
                f_deltaM_ *= (f*constant_forward**stream + (1-f)*constant_back**stream)

                for l in range(1,stream):
                    w = (2*l+1) * (f*g_forward**l + (1-f)*g_back**l)
                    if w_single_form==0:
                        w_single[l,:,:] = (w - (2*l+1)*f_deltaM_) / (1 - f_deltaM_)
                    if w_multi_form==0:
                        w_multi[l,:,:] = (w - (2*l+1)*f_deltaM_) / (1 - f_deltaM_)

            if w_single_rayleigh==1:
                w_single[1:] *= ftau_cld
                if stream==4:
                    w_single[2] += 0.5*ftau_ray 
            if w_multi_rayleigh==1: 
                w_multi[1:] *= ftau_cld
                if stream==4:
                    w_multi[2] += 0.5*ftau_ray 

            #single-scattering options
            if single_form==0: # explicit single form
                if psingle_form==1: #OTHG
                    p_single=(1-cosb_og**2)/(sqrt(1+cosb_og**2+2*cosb_og*cos_theta)**3) 
                elif psingle_form==0: #'TTHG':
                    g_forward = constant_forward*cosb_og
                    g_back = constant_back*cosb_og
                    f = frac_a + frac_b*g_back**frac_c
                    p_single=(f * (1-g_forward**2) /sqrt((1+g_forward**2+2*g_forward*cos_theta)**3) 
                                #second term of TTHG: backward scattering
                                +(1-f)*(1-g_back**2) /sqrt((1+g_back**2+2*g_back*cos_theta)**3))

                if psingle_rayleigh==1: 
                    p_single = ftau_cld*p_single + ftau_ray*(0.75*(1+cos_theta**2.0))


            for l in range(stream):
                a[l,:,:] = (2*l + 1) -  w0 * w_multi[l,:,:]
                b[l,:,:] = ( F0PI * (w0 * w_single[l,:,:])) * Pu0[l] / (4*pi)

            #boundary conditions 
            b_surface = (0. + surf_reflect*u0*F0PI*exp(-tau[-1, :]/u0))#/(2*pi)
            b_surface_SH4 = -(0. + surf_reflect*u0*F0PI*exp(-tau[-1, :]/u0))/4#/(2*pi)# need to double check BCs)
            b_top_ = b_top#/(2*pi)

            if stream==2:
                M, B, F_bot, G_bot, F, G, Q1, Q2, lam, q, eta  = setup_2_stream_fluxes(nlayer, nwno, w0, b_top_, b_surface, 
                surf_reflect, u0, dtau, tau, a, b, fluxes=flx, calculation=0) 

            if stream==4:
                M, B, F_bot, G_bot, F, G, lam1, lam2, A, eta = setup_4_stream_fluxes(nlayer, nwno, w0, b_top_, b_surface, b_surface_SH4, 
                    surf_reflect, u0, dtau, tau, a, b, fluxes=flx, calculation=0) 

                # F and G will be nonzero if fluxes=1

            flux_bot = zeros(nwno)
            intgrl_new = zeros((stream*nlayer, nwno))
            flux_temp = zeros((stream*nlevel, nwno))
            intgrl_per_layer = zeros((nlayer, nwno))
            xint_temp = zeros((nlevel, nwno))
            multi_scat = zeros((nlayer, nwno))

            #========================= Start loop over wavelength =========================
            X = zeros((stream*nlayer, nwno))
            flux_temp = zeros((stream*nlevel, nwno))
            for W in range(nwno):
                X[:,W] = solve_4_stream_banded(M[:,:,W], B[:,W], stream)
                if flx==1:
                    flux_temp[:,W] = calculate_flux(F[:,:,W], G[:,W], X[:,W])
            flux_bot = np.sum(F_bot*X, axis=0) + G_bot
    
    return xint_at_top, flux#, xint_out

def setup_2_stream_fluxes(nlayer, nwno, w0, b_top, b_surface, surf_reflect, ubar0, 
        dtau, tau, a, b, B0=0., B1=0., fluxes=0, calculation=0):#'reflected'):
    """
    Setup up matrices to solve flux problem for spherical harmonics method.

    Parameters
    ----------
    nlayer : int 
        Number of layers
    nwno : int 
        Number of wavenumber points 
    w0 : numpy.ndarray
        This is a matrix of nlayer by nwave. This describes the single scattering albedo of 
        the atmosphere. Note this is free of any Raman scattering or any d-eddington correction 
        that is sometimes included in reflected light calculations.
    b_top : array 
        The diffuse radiation into the model at the top of the atmosphere
    b_surface : array
        The diffuse radiation into the model at the bottom. Includes emission, reflection 
        of the unattenuated portion of the direct beam  
    b_surface_SH4 : array
        Second bottom BC for SH4 method.
    surf_reflect : numpy.ndarray    
        Surface reflectivity as a function of wavenumber. 
    ubar0 : ndarray of float 
        matrix of cosine of the incident angle from geometric.json
    dtau : numpy.ndarray
        This is a matrix of nlayer by nwave. This describes the per layer optical depth. 
    tau : numpy.ndarray
        This is a matrix of nlevel by nwave. This describes the cumulative optical depth. 
    a: numpy.ndarray
        Coefficients of matrix capturing legendre expansion of phase function for multiple scattering
    b: numpy.ndarray
        Coefficients of source vector capturing legendre expansion of phase function for single 
        scattering
    B0 : numpy.ndarray
        Matrix of blackbodies
    B1 : numpy.ndarray
        Eqn (26) Toon 89
    fluxes : int 
        Toggle calculation of layerwise fluxes (0 = do not calculate, 1 = calculate)
    calculation : int 
        Toggle calculation method (1 = linear, 2 = exponential)

    Returns
    -------
    numpy.ndarrays
       Matrices and vectors used to calculate fluxes and intensities at each level 
    """

    eta = zeros((2, nlayer, nwno)) # will remain zero for thermal
    if calculation==0: #reflected light
        Del = ((1 / ubar0)**2 - a[0]*a[1])
        eta = [(b[1] /ubar0 - a[1]*b[0]) / Del,
            (b[0] /ubar0 - a[0]*b[1]) / Del]

    lam = sqrt(a[0]*a[1])
    expo = lam*dtau
    expo = slice_rav(expo, 35.0) 
    exptrm = exp(-expo)

    #   parameters in matrices
    q = lam/a[1]
    Q1 = (0.5 + q)*2*pi
    Q2 = (0.5 - q)*2*pi

    Q1mn = Q1*exptrm;  Q2mn = Q2*exptrm
    Q1pl = Q1/exptrm;  Q2pl = Q2/exptrm

    if calculation == 0: #reflected light
        zmn = (0.5*eta[0] - eta[1])*2*pi
        zpl = (0.5*eta[0] + eta[1])*2*pi
        expon = exp(-tau/ubar0)
        zmn_up = zmn * expon[1:,:] 
        zpl_up = zpl * expon[1:,:] 
        zmn_down = zmn * expon[:-1,:] 
        zpl_down = zpl * expon[:-1,:] 
    elif calculation == 1: # linear thermal
        zmn_down = ((1-w0)/a[0] * (B0/2 - B1/a[1])) *2*pi           #* 2*pi
        zmn_up = ((1-w0)/a[0] * (B0/2 - B1/a[1] + B1*dtau/2)) *2*pi #* 2*pi
        zpl_down = ((1-w0)/a[0] * (B0/2 + B1/a[1])) *2*pi           #* 2*pi
        zpl_up = ((1-w0)/a[0] * (B0/2 + B1/a[1] + B1*dtau/2)) *2*pi #* 2*pi


    #   construct matrices
    Mb = zeros((5, 2*nlayer, nwno))
    B = zeros((2*nlayer, nwno))
    nlevel = nlayer+1
    F = zeros((2*nlevel, 2*nlayer, nwno))
    G = zeros((2*nlevel, nwno))

    #   first row: BC 1
    Mb[2,0,:] = Q1[0,:]
    Mb[1,1,:] = Q2[0,:]
    B[0,:] = b_top - zmn_down[0,:]

    #   last row: BC 4
    n = nlayer-1
    Mb[3, 2*nlayer-2,:] = Q2mn[n,:] - surf_reflect*Q1mn[n,:]
    Mb[2, 2*nlayer-1,:] = Q1pl[n,:] - surf_reflect*Q2pl[n,:]
    B[2*nlayer-1,:] = b_surface - zpl_up[n,:] + surf_reflect * zmn_up[n,:]

    # remaining rows
    Mb[0,3::2,:] = -Q2[1:,:]
    Mb[1,2::2,:] = -Q1[1:,:]
    Mb[1,3::2,:] = -Q1[1:,:]
    Mb[2,1:-1:2,:] = Q2pl[:-1,:]
    Mb[2,2::2,:] = -Q2[1:,:]
    Mb[3,:-2:2,:] = Q1mn[:-1,:]
    Mb[3,1:-1:2,:] = Q1pl[:-1,:]
    Mb[4,:-2:2,:] = Q2mn[:-1,:]
    B[1:-1:2,:] = zmn_down[1:,:] - zmn_up[:-1,:]
    B[2::2,:] = zpl_down[1:,:] - zpl_up[:-1,:]


    # flux at bottom of atmosphere
    F_bot = zeros((2*nlayer, nwno))
    G_bot = zeros(nwno)
    F_bot[-2,:] = Q2mn[-1,:]
    F_bot[-1,:] = Q1pl[-1,:]
    G_bot = zpl_up[-1,:]

    if fluxes == 1: # fluxes per layer
        F[0,0,:] = Q1[0,:]
        F[0,1,:] = Q2[0,:]
        F[1,0,:] = Q2[0,:]
        F[1,1,:] = Q1[0,:]

        nn = np.arange(2*nlayer)
        indcs = nn[::2]
        k = 0
        for i in indcs:
            F[i+2,i,:] = Q1mn[k,:]
            F[i+2,i+1,:] = Q2pl[k,:]
            F[i+3,i,:] = Q2mn[k,:]
            F[i+3,i+1,:] = Q1pl[k,:]
            k = k+1

        G[0,:] = zmn_down[0,:]
        G[1,:] = zpl_down[0,:]

        G[2::2,:] = zmn_up
        G[3::2,:] = zpl_up

    return Mb, B, F_bot, G_bot, F, G, Q1, Q2, lam, q, eta





'''