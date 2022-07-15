# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:06:10 2019

@author: Ryan Boukrouche
"""
import os
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.ticker
from matplotlib.offsetbox import AnchoredText
from matplotlib.lines import Line2D
from itertools import groupby
from operator import itemgetter
#from moviepy.editor import VideoClip
#from moviepy.video.io.bindings import mplfig_to_npimage
import seaborn as sns

""" ----------------- Define the number of shortwave bands ---------------- """

nSWb = "3SWb"

""" -------------- Define the surface pressure to search for -------------- """

surfpres = "11"

""" ------------------ Define the input and output paths ------------------ """

path                 = "./data_output/data_output_"+nSWb+"_"+surfpres+"/"
path_general_adiabat = "./general_adiabat/output_"+surfpres+"/data/"
path_plots           = "./plots_output"
if not os.path.exists(path_plots):
    os.makedirs(path_plots)
path_plots           = "./plots_output/plots_output_"+nSWb+"_"+surfpres+"/"
if not os.path.exists(path_plots):
    os.makedirs(path_plots)

""" ------------------ Was it a time-stepped calculation? ----------------- """

timestep             = True

""" ------------------ Find the number of timesteps used ------------------ """

path_pT_Ts01         = './data_output/data_output_'+nSWb+"_"+surfpres+'/pT_profiles/pT_Ts01'
path_pT_Ts           = './data_output/data_output_'+nSWb+"_"+surfpres+'/pT_profiles/pT_Ts'
#n_files              = len([name for name in os.listdir(path_pT_Ts01) if os.path.isfile(os.path.join(path_pT_Ts01, name))])

""" --------- Was it a pure grey calculation or a band-grey one? ---------- """

Grey                 = False

""" -------------------------- Plot animations? --------------------------- """

Animations           = True     

""" Define indices for plotting single steps or single surface temperatures """

Ts_index             = -1
step_index           = -1

""" ---------------- Define the surface temperature array ----------------- """

#if (len(np.loadtxt(path + "TsL.dat",skiprows=1,ndmin=1)) == 1):  
#    TsL                  = np.loadtxt(path + "TsL.dat",skiprows=1)
#else:
#TsL                  = np.loadtxt(path + "TsL.dat",skiprows=1,ndmin=1)
#TsL                  = np.arange(200.0, 3001.0, 50.0) 
TsL = np.array([800.0])

if (len(TsL) > 1):
    for Tsi in range(len(TsL)):    
        path_plots_PT_Ts    = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/PT_Ts/'
        if not os.path.exists(path_plots_PT_Ts):
            os.makedirs(path_plots_PT_Ts)
        path_plots_PT_Ts    = path_plots_PT_Ts + 'PT_Ts_'+f'{Tsi}/'
        if not os.path.exists(path_plots_PT_Ts):
            os.makedirs(path_plots_PT_Ts)
        path_plots_flux_Ts  = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/flux_Ts/'
        if not os.path.exists(path_plots_flux_Ts):
            os.makedirs(path_plots_flux_Ts)
        path_plots_flux_Ts  = path_plots_flux_Ts + 'flux_Ts_'+f'{Tsi}/'
        if not os.path.exists(path_plots_flux_Ts):
            os.makedirs(path_plots_flux_Ts)
        path_plots_tau_Ts   = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/tau_Ts/'
        if not os.path.exists(path_plots_tau_Ts):
            os.makedirs(path_plots_tau_Ts)
        path_plots_tau_Ts   = path_plots_tau_Ts + 'tau_Ts_'+f'{Tsi}/'
        if not os.path.exists(path_plots_tau_Ts):
            os.makedirs(path_plots_tau_Ts)
        path_plots_OPR_Ts   = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/OPR_Ts/'
        if not os.path.exists(path_plots_OPR_Ts):
            os.makedirs(path_plots_OPR_Ts)
        path_plots_OPR_Ts   = path_plots_OPR_Ts + 'OPR_Ts_'+f'{Tsi}/'
        if not os.path.exists(path_plots_OPR_Ts):
            os.makedirs(path_plots_OPR_Ts)
        path_plots_heat_Ts    = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/heat_Ts/'
        if not os.path.exists(path_plots_heat_Ts):
            os.makedirs(path_plots_heat_Ts)
        path_plots_heat_Ts    = path_plots_PT_Ts + 'heat_Ts_'+f'{Tsi}/'
        if not os.path.exists(path_plots_heat_Ts):
            os.makedirs(path_plots_heat_Ts)
else:
    path_plots_PT_Ts    = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/PT_Ts/'
    if not os.path.exists(path_plots_PT_Ts):
        os.makedirs(path_plots_PT_Ts)
    path_plots_PT_Ts    = path_plots_PT_Ts + 'PT_Ts_00/'
    if not os.path.exists(path_plots_PT_Ts):
        os.makedirs(path_plots_PT_Ts)
    path_plots_flux_Ts  = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/flux_Ts/'
    if not os.path.exists(path_plots_flux_Ts):
        os.makedirs(path_plots_flux_Ts)
    path_plots_flux_Ts  = path_plots_flux_Ts + 'flux_Ts_00/'
    if not os.path.exists(path_plots_flux_Ts):
        os.makedirs(path_plots_flux_Ts)
    path_plots_tau_Ts   = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/tau_Ts/'
    if not os.path.exists(path_plots_tau_Ts):
        os.makedirs(path_plots_tau_Ts)
    path_plots_tau_Ts   = path_plots_tau_Ts + 'tau_Ts_00/'
    if not os.path.exists(path_plots_tau_Ts):
        os.makedirs(path_plots_tau_Ts)
    path_plots_OPR_Ts   = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/OPR_Ts/'
    if not os.path.exists(path_plots_OPR_Ts):
        os.makedirs(path_plots_OPR_Ts)
    path_plots_OPR_Ts   = path_plots_OPR_Ts + 'OPR_Ts_00/'
    if not os.path.exists(path_plots_OPR_Ts):
        os.makedirs(path_plots_OPR_Ts)
    path_plots_heat_Ts    = "./plots_output/plots_output_"+nSWb+"_"+surfpres+'/heat_Ts/'
    if not os.path.exists(path_plots_heat_Ts):
        os.makedirs(path_plots_heat_Ts)
    path_plots_heat_Ts    = path_plots_heat_Ts + 'heat_Ts_00/'
    if not os.path.exists(path_plots_heat_Ts):
        os.makedirs(path_plots_heat_Ts)
            
        
    
n_files              = np.array([len([name for name in os.listdir(path_pT_Ts+f'{j:02}') if os.path.isfile(os.path.join(path_pT_Ts+f'{j:02}', name))]) for j in range(1,len(TsL)+1)])
#n_files is the array of 1800, length 57
fps                  = 1
Rcp                  = 0.25009059689312091
psurf                = float(surfpres+'e5')

def planck_wl(wl, T):
    h = 6.626075540e-34    #Planck's constant
    c = 2.99792458e8       #Speed of light
    k = 1.38065812e-23     #Boltzman thermodynamic constant
    a = (2.*h*c**2)/wl**5
    b = (h*c)/(wl*k*T)
    intensity = a/ (  (np.exp(b) - 1.0) )
    return intensity

""" ----------------------- Define the other arrays ----------------------- 

if timestep and n_step > 1:
    if isinstance(TsL, float) == True:
        TsL = np.array([TsL])
    dt_array=np.array([np.loadtxt(path + "pT_profiles/pT_Ts" + f'{j:02}' + "/steps.txt",skiprows=1,usecols=0) for j in range(1,len(TsL)+1)])
    n_step = np.shape([dt_array])[-1]
    Ts_array=np.array([np.loadtxt(path + "pT_profiles/pT_Ts" + f'{j:02}' + "/steps.txt",skiprows=1,usecols=1) for j in range(1,len(TsL)+1)])
    #dt_array=dt_array[np.nonzero(dt_array)] # if an element is 0 it's just not counted in steps
    #Ts_array=Ts_array[np.nonzero(Ts_array)]
    #enthalpy  = np.array([np.loadtxt(path + "/pT_profiles/pT_Ts" + f'{j:02}' + "/steps.txt",skiprows=1,usecols=2) for j in range(1,len(TsL)+1)])
    steps=np.zeros([len(TsL),n_step])
    #print("Ts_array = ",Ts_array)
    for j in range(len(TsL)):
        #print("steps = ",steps)
        #print("np.shape([dt_array])",np.shape([dt_array]))
        if len(TsL) == 1 and n_step == 1: # if there is only one Ts value, and only one step
            dt_array = np.array([dt_array]) # add a dimension
            Ts_array = np.array([Ts_array]) 
            #print("Badger")
        #print(dt_array[0,:])
        #print(dt_array[:,0])
        #print("Ts_array = ",Ts_array)
        #print(np.shape(Ts_array))
        steps[j,0]=dt_array[j,0]
        for i in range(n_step-1):
            steps[j,i+1]=steps[j,i]+dt_array[j,i]
"""

first_step_c  = "0000010800" 
first_step    = int(first_step_c)    # first step registered
step_interval = 10800    
n_step        = n_files*step_interval

pmid          = np.array([np.array([np.loadtxt(path + "pT_profiles/pT_Ts" + f'{j:02}' + "/p-T" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
pint          = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])                      
Tmid          = np.array([np.array([np.loadtxt(path + "pT_profiles/pT_Ts" + f'{j:02}' + "/p-T" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])       
Dry_adiabat   = np.array([TsL[j]*(pmid[j,:,:]/psurf)**Rcp for j in range(len(TsL))])
Moist_adiabat = np.array([np.loadtxt(path + "pT_profiles/pT_Ts" + f'{j:02}' + "/p-T"+ first_step_c + ".dat",skiprows=1,usecols=2) for j in range(1,len(TsL)+1)])

cffIR         = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
cffW1         = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
cffW2         = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "1SWb"):
    cffSW         = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    cff = cffIR + cffW1 + cffW2 + cffSW
elif (nSWb == "2SWb"):
    cffUV         = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    cffVIS        = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    cff = cffIR + cffW1 + cffW2 + cffUV + cffVIS
elif (nSWb == "3SWb"):
    cffUV         = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    cffVIS1       = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    cffVIS2       = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/cff" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    cff = cffIR + cffW1 + cffW2 + cffUV + cffVIS1 + cffVIS2
    
scffIR        = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
scffW1        = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
scffW2        = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "1SWb"):
    scffSW        = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    scff = scffIR + scffW1 + scffW2 + scffSW
if (nSWb == "2SWb"):
    scffUV        = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    scffVIS       = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    scff = scffIR + scffW1 + scffW2 + scffUV + scffVIS
if (nSWb == "3SWb"):
    scffUV        = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    scffVIS1      = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    scffVIS2      = np.array([np.array([np.loadtxt(path + "cff/cff_Ts" + f'{j:02}' + "/scff" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    scff = scffIR + scffW1 + scffW2 + scffUV + scffVIS1 + scffVIS2

fluxIR        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
fluxW1        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
fluxW2        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "1SWb"):
    fluxSW        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "2SWb"):
    fluxUV        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    fluxVIS       = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "3SWb"):
    fluxUV        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    fluxVIS1      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    fluxVIS2      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Total_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=6) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])

diffusefluxIR        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
diffusefluxW1        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
diffusefluxW2        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "1SWb"):
    diffusefluxSW        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    diffuseflux = diffusefluxIR + diffusefluxW1 + diffusefluxW2 + diffusefluxSW
if (nSWb == "2SWb"):
    diffusefluxUV        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    diffusefluxVIS       = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    diffuseflux = diffusefluxIR + diffusefluxW1 + diffusefluxW2 + diffusefluxUV + diffusefluxVIS
if (nSWb == "3SWb"):
    diffusefluxUV        = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    diffusefluxVIS1      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    diffusefluxVIS2      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/Diffuse_net_flux" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    diffuseflux = diffusefluxIR + diffusefluxW1 + diffusefluxW2 + diffusefluxUV + diffusefluxVIS1 + diffusefluxVIS2

IplusIR       = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
IplusW1       = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])                                
IplusW2       = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])                                
if (nSWb == "1SWb"):
    IplusSW       = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Iplus = IplusIR + IplusW1 + IplusW2 + IplusSW
    OLR = IplusIR[:,:,0] + IplusW1[:,:,0] + IplusW2[:,:,0]
    OSR = IplusSW[:,:,0]
    OPR = OLR + OSR
if (nSWb == "2SWb"):
    IplusUV       = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    IplusVIS      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Iplus = IplusIR + IplusW1 + IplusW2 + IplusUV + IplusVIS
    OLR = IplusIR[:,:,0] + IplusW1[:,:,0] + IplusW2[:,:,0]
    OSR = IplusUV[:,:,0] + IplusVIS[:,:,0]
    OPR = OLR + OSR
if (nSWb == "3SWb"):
    IplusUV       = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    IplusVIS1     = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    IplusVIS2     = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_up" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Iplus = IplusIR + IplusW1 + IplusW2 + IplusUV + IplusVIS1 + IplusVIS2
    OLR = IplusIR[:,:,0] + IplusW1[:,:,0] + IplusW2[:,:,0]
    OSR = IplusUV[:,:,0] + IplusVIS1[:,:,0] + IplusVIS2[:,:,0]
    OPR = OLR + OSR

IminusIR      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
IminusW1      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
IminusW2      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "1SWb"):
    IminusSW      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Iminus = IminusIR + IminusW1 + IminusW2 + IminusSW
if (nSWb == "2SWb"):
    IminusUV      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    IminusVIS     = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Iminus = IminusIR + IminusW1 + IminusW2 + IminusUV + IminusVIS
if (nSWb == "3SWb"):
    IminusUV      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    IminusVIS1     = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    IminusVIS2     = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/diffuse_down" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Iminus = IminusIR + IminusW1 + IminusW2 + IminusUV + IminusVIS1 + IminusVIS2

Direct_beamIR      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
Direct_beamW1      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
Direct_beamW2      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "1SWb"):
    Direct_beamSW      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Direct_beam = Direct_beamIR + Direct_beamW1 + Direct_beamW2 + Direct_beamSW 
if (nSWb == "2SWb"):
    Direct_beamUV      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Direct_beamVIS     = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Direct_beam = Direct_beamIR + Direct_beamW1 + Direct_beamW2 + Direct_beamUV + Direct_beamVIS
if (nSWb == "3SWb"):
    Direct_beamUV      = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Direct_beamVIS1     = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Direct_beamVIS2     = np.array([np.array([np.loadtxt(path + "fluxes/flux_Ts" + f'{j:02}' + "/direct_beam" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    Direct_beam = Direct_beamIR + Direct_beamW1 + Direct_beamW2 + Direct_beamUV + Direct_beamVIS1 + Direct_beamVIS2

"""
    Planck_total  = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=1) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
    Planck_IR      = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=2) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
    Planck_W1      = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=3) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
    Planck_W2      = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=4) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
    if (nSWb == "1SWb"):
        Planck_SW      = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=5) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
    if (nSWb == "2SWb"):
        Planck_UV      = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=5) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
        Planck_VIS     = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=6) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
    if (nSWb == "3SWb"):
        Planck_UV      = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=5) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
        Planck_VIS1    = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=6) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
        Planck_VIS2    = np.array([np.array([np.loadtxt(path + "fluxes/B_Ts" + f'{j:02}' + "/B" + f'{i:04}' + ".dat",skiprows=1,usecols=7) for i in range(1,n_step+1)]) for j in range(1,len(TsL)+1)])
"""

tau_total     = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=0) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
tau_IR         = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
tau_W1         = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
tau_W2         = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "1SWb"):
    tau_SW         = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "2SWb"):
    tau_UV         = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    tau_VIS        = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
if (nSWb == "3SWb"):
    tau_UV         = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    tau_VIS1       = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
    tau_VIS2       = np.array([np.array([np.loadtxt(path + "tau/tau_Ts" + f'{j:02}' + "/tau" + f'{i:010}' + ".dat",skiprows=1,usecols=6) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])

dT_rad        = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=1) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_conv       = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=2) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_conv_moist = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=3) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_total      = dT_rad + dT_conv + dT_conv_moist
dT_rad_planetary      = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=4) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_rad_stellar      = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=5) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_rad_IR      = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=6) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_rad_W1      = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=7) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_rad_W2      = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=8) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_rad_UV      = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=9) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_rad_VIS1      = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=10) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])
dT_rad_VIS2      = np.array([np.array([np.loadtxt(path + "heating_rates/heat_Ts" + f'{j:02}' + "/heat" + f'{i:010}' + ".dat",skiprows=1,usecols=11) for i in range(first_step,n_step[j-1]+1,step_interval)]) for j in range(1,len(TsL)+1)])

def literature_comparison(ax):
   
    Goldblatt13_Ts  = []
    Goldblatt13_OLR = []
    with open("./comparison_data/Goldblatt13_data.txt", 'r') as data_file:
        for line in data_file:
            if not line.startswith('#'):
                line = line.rstrip('\n')
                line = line.split(",")
                Goldblatt13_Ts.append(float(line[0]))
                Goldblatt13_OLR.append(float(line[1]))
    Kopparapu13_Ts  = []
    Kopparapu13_OLR = []
    with open("./comparison_data/Kopparapu13_data.txt", 'r') as data_file:
        for line in data_file:
            if not line.startswith('#'):
                line = line.rstrip('\n')
                line = line.split(",")
                Kopparapu13_Ts.append(float(line[0]))
                Kopparapu13_OLR.append(float(line[1]))
    Hamano15_Ts  = []
    Hamano15_OLR = []
    with open("./comparison_data/Hamano15_data.txt", 'r') as data_file:
        for line in data_file:
            if not line.startswith('#'):
                line = line.rstrip('\n')
                line = line.split(",")
                Hamano15_Ts.append(float(line[0]))
                Hamano15_OLR.append(float(line[1]))

    ### Plot and annotate literature comparison
    ax.plot(Goldblatt13_Ts, Goldblatt13_OLR, color='#646464', ls=":", lw=1.0, zorder=0.1)
    ax.text(1900, 320, "Goldblatt+ 13", va="bottom", ha="right", fontsize=7, color='#646464', bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    # ax.plot(Kopparapu13_Ts, Kopparapu13_OLR, color='#646464', ls="-.", lw=1.0, zorder=0.1)
    ax.plot(Hamano15_Ts, Hamano15_OLR, color='#646464', ls="-.", lw=1.0, zorder=0.1)
    # ax.text(2180, 330, "Hamano+ 15", va="top", ha="left", fontsize=7, color='#646464', bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    ax.text(2180, 350, "Kopparapu+ 13 / Hamano+ 15", va="top", ha="left", fontsize=7, color='#646464', bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    # ax.text(1500, 330, "Kopparapu+ 13", va="top", ha="left", fontsize=7, color=#646464, bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))


""" *********************************************************************** """
""" ------- Plotting temperature profiles and contribution functions ------ """
""" *********************************************************************** """
""" Only for an array of 57 Ts values """

if (len(TsL) == 57):
    TsL_plot_indices = [8,26,42,56]
    fig, (ax1, ax2) = pl.subplots(1, 2, figsize=(13,6))#(22,10))
    fig.subplots_adjust(wspace=0)
    colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
    for Tsi in TsL_plot_indices:
        Ts = TsL[Tsi]
        color=next(colors)
        ax1.semilogy(Tmid[Tsi,step_index,:],pmid[Tsi,step_index,:], color=color, label=f'{TsL[Tsi]}')
        ax1.semilogy(Moist_adiabat[Tsi,:],pmid[Tsi,step_index,:],'r--')
        #ax1.text(2500., 1e5, 'A', fontsize=14, family='serif')
        ax1.set_xlabel(r"Temperature, $T$ [K]")
        ax1.set_ylabel(r"Pressure, $P$ [Pa]")
        ax1.invert_yaxis()
        ax1.set_xlim(0.,3000.)
        ax1.set_ylim(pmid[-1,step_index,-1],1.)
        sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
        ax1.legend(frameon=True,loc='best',title=r'$T_\mathrm{surf}$ (K)')#,fontsize=14)

    colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
    for Tsi in TsL_plot_indices:
        color=next(colors)
        ax2.loglog(cff[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color=color)       
        if (scff[Tsi,step_index] > 0.):
            ax2.plot(scff[Tsi,step_index], pmid[Tsi,step_index,-1], marker='^',color=color, label=f'{Ts}')
            
    ax2.set_xlabel(r"Flux contribution, $\mathcal{CF}_\mathrm{F} \: [\mathrm{W} \: \mathrm{m}^{-2} \mathrm{dex}^{-1}]$")
    #ax2.text(1e0, 1e6, 'B', fontsize=14, family='serif')
    #ax2.set_ylabel(r"Pressure, $P$ [Pa]")
    
    # Minor ticks
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12)
    ax2.yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
    ax2.yaxis.set_minor_locator(locmin)
    ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    
    ax2.invert_yaxis()
    ax2.set_ylim(pmid[-1,step_index,-1],1.)
    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
    ax2.legend(frameon=True,loc='best',title=r'$\mathcal{SCF}_\mathrm{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$')
    
    #custom_lines = [Line2D([0], [0], color='k', ls='-', lw=1),
    #                Line2D([0], [0], color='k', ls=':', lw=1),
    #                Line2D([0], [0], color='k', ls='-.', lw=1),
    #                Line2D([0], [0], color='k', ls='--', lw=1)]
    
    #x2labels = ax2.get_xticklabels()
    #x2labels[0] = "" # remove the first tick label
    
    # Remove the y-axis labels
    labels = [item.get_text() for item in ax2.get_yticklabels()]
    
    empty_string_labels = ['']*len(labels)
    ax2.set_yticklabels(empty_string_labels)
    
    # Remove the first x-axis label
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    string_labels = ['','',r'$\mathrm{10^{-5}}$',r'$\mathrm{10^{-4}}$',r'$\mathrm{10^{-3}}$',r'$\mathrm{10^{-2}}$',r'$\mathrm{10^{-1}}$',r'$\mathrm{10^{0}}$',r'$\mathrm{10^{1}}$']
    ax2.set_xticklabels(string_labels)
    
    #ax2.legend(custom_lines, ['SW', 'IR', 'W1', 'W2'],title='Spectral region')
    
    """ ======================= General figure settings ======================= """
    
    ax1.grid(axis='y')
    ax2.grid(axis='y')
    
    FigAbox = AnchoredText('A', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
    FigAbox.patch.set(alpha=0.5)
    FigBbox = AnchoredText('B', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
    FigBbox.patch.set(alpha=0.5)
    ax1.add_artist(FigAbox)
    ax2.add_artist(FigBbox)
    
    fig.savefig(path_plots+'PT_CFF.pdf',bbox_inches='tight')
    #pl.show()
    pl.close()
    
else:
    
    fig, (ax1, ax2) = pl.subplots(1, 2, figsize=(13,6))#(22,10))
    fig.subplots_adjust(wspace=0)
    Tsi = Ts_index
    Ts = TsL[Tsi]
    ax1.semilogy(Tmid[Tsi,step_index,:],pmid[Tsi,step_index,:], color='royalblue', label='Midpoint temperature')
    ax1.semilogy(Moist_adiabat[Tsi,:],pmid[Tsi,step_index,:],'r--', label='Dew point')
    #ax1.text(2500., 1e5, 'A', fontsize=14, family='serif')
    ax1.set_xlabel(r"Temperature, $T$ [K]")
    ax1.set_ylabel(r"Pressure, $P$ [Pa]")
    ax1.invert_yaxis()
    ax1.set_xlim(0.,3000.)
    ax1.set_ylim(pmid[-1,step_index,-1],1.)
    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
    ax1.legend(frameon=True,loc='best')#,fontsize=14)
    
    ax2.loglog(cff[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='royalblue')
    if (scff[Tsi,step_index] > 0.):
        ax2.plot(scff[Tsi,step_index], pmid[Tsi,step_index,-1], marker='o',color='royalblue', label=f'{Ts}')
        
    ax2.set_xlabel(r"Flux contribution, $\mathcal{CF}_\mathrm{F} \: [\mathrm{W} \: \mathrm{m}^{-2} \mathrm{dex}^{-1}]$")
    #ax2.text(1e0, 1e6, 'B', fontsize=14, family='serif')
    #ax2.set_ylabel(r"Pressure, $P$ [Pa]")
    
    # Minor ticks
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12)
    ax2.yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
    ax2.yaxis.set_minor_locator(locmin)
    ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    
    ax2.invert_yaxis()
    ax2.set_ylim(pmid[-1,step_index,-1],1.)
    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
    #ax2.legend(frameon=True,loc='best',title=r'$\mathcal{SCF}_\mathrm{F} \: [\mathrm{W} \: \mathrm{m}^{-2} \mathrm{dex}^{-1}]$')
    
    #custom_lines = [Line2D([0], [0], color='k', ls='-', lw=1),
    #                Line2D([0], [0], color='k', ls=':', lw=1),
    #                Line2D([0], [0], color='k', ls='-.', lw=1),
    #                Line2D([0], [0], color='k', ls='--', lw=1)]
    
    #x2labels = ax2.get_xticklabels()
    #x2labels[0] = "" # remove the first tick label
    
    # Remove the y-axis labels
    labels = [item.get_text() for item in ax2.get_yticklabels()]
    
    empty_string_labels = ['']*len(labels)
    ax2.set_yticklabels(empty_string_labels)
    
    # Remove the first x-axis label
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    string_labels = ['','',r'$\mathrm{10^{-5}}$',r'$\mathrm{10^{-4}}$',r'$\mathrm{10^{-3}}$',r'$\mathrm{10^{-2}}$',r'$\mathrm{10^{-1}}$',r'$\mathrm{10^{0}}$',r'$\mathrm{10^{1}}$']
    ax2.set_xticklabels(string_labels)
    
    #ax2.legend(custom_lines, ['SW', 'IR', 'W1', 'W2'],title='Spectral region')
    
    """ ======================= General figure settings ======================= """
    
    ax1.grid(axis='y')
    ax2.grid(axis='y')
    
    FigAbox = AnchoredText('A', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
    FigAbox.patch.set(alpha=0.5)
    FigBbox = AnchoredText('B', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
    FigBbox.patch.set(alpha=0.5)
    ax1.add_artist(FigAbox)
    ax2.add_artist(FigBbox)
    
    fig.savefig(path_plots+'PT_CFF.pdf',bbox_inches='tight')
    #pl.show()
    pl.close()


""" *********************************************************************** """
""" --------------- Plotting flux and heating rate profiles --------------- """
""" *********************************************************************** """

fig, ((ax1, ax3), (ax2, ax4)) = pl.subplots(2, 2, figsize=(13,6))#(22,10))
fig.subplots_adjust(wspace=0.)
Tsi   = Ts_index # for a single surface temperature
color = 'b'
""" ========================== A1 ========================== """

#colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
#for Tsi in TsL_plot_indices:
#    color=next(colors)   
ax1.semilogy(diffuseflux[Tsi,step_index,:],pint[Tsi,step_index,:], color='b', lw=2, label='Net flux')
ax1.semilogy(Iplus[Tsi,step_index,:],pint[Tsi,step_index,:], color='c', ls='-', label='Upward diffuse flux')
ax1.semilogy(Iminus[Tsi,step_index,:],pint[Tsi,step_index,:], color='orange', ls='-', label='Downward diffuse flux')
ax1.semilogy(Direct_beam[Tsi,step_index,:],pint[Tsi,step_index,:], color='r', ls='-', label='Direct beam')
#ax1.text(2500., 1e6, 'A1', fontsize=14, family='serif')
ax1.set_xlabel(r"Diffuse net flux, $\mathcal{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$")
ax1.set_ylabel(r"Pressure, $P$ [Pa]")
ax1.invert_yaxis()
pl.xscale('symlog')
#ax1.set_xlim(-1e-2,1e7)
ax1.set_ylim(pint[-1,step_index,-1],1.)
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
ax1.legend(frameon=True,loc='best')#,fontsize=14)

""" ========================== B1 ========================== """

#colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
#for Tsi in TsL_plot_indices:
#    color=next(colors)
ax2.semilogy(IplusIR[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-.',color='lightblue', label='W1')
ax2.semilogy(IplusW1[Tsi,step_index,:],pint[Tsi,step_index,:],ls='--',color='royalblue', label='W2')
ax2.semilogy(IplusW2[Tsi,step_index,:],pint[Tsi,step_index,:],ls=':',color='navy', label='IR')
if (nSWb == "1SWb"):
    ax2.semilogy(IplusSW[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='SW')
if (nSWb == "2SWb"):
    ax2.semilogy(IplusUV[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='UV')
    ax2.semilogy(IplusVIS[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS')
if (nSWb == "3SWb"):
    ax2.semilogy(IplusUV[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='UV')
    ax2.semilogy(IplusVIS1[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS1')
    ax2.semilogy(IplusVIS2[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS2')

ax2.set_xlabel(r"Upward flux, $\mathcal{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$")
ax2.set_ylabel(r"Pressure, $P$ [Pa]")
ax2.invert_yaxis()

#ax2.set_xlim(1e-3,1e7)
ax2.set_ylim(pint[-1,step_index,-1],1.)
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
ax2.legend(frameon=True,loc='best',title='Spectral band')

""" ========================== A2 ========================== """

#colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
#for Tsi in TsL_plot_indices:
#    color=next(colors)
pl.xscale('symlog')
ax3.semilogy(Iminus[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color=color)
ax3.semilogy(IminusIR[Tsi,step_index,:],pint[Tsi,step_index,:],ls=':',color='navy', label='IR')
ax3.semilogy(IminusW1[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-.',color='lightblue', label='W1')
ax3.semilogy(IminusW2[Tsi,step_index,:],pint[Tsi,step_index,:],ls='--',color='royalblue', label='W2')
if (nSWb == "1SWb"):
    ax3.semilogy(IminusSW[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='SW')
if (nSWb == "2SWb"):
    ax3.semilogy(IminusUV[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='UV')
    ax3.semilogy(IminusVIS[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS')
if (nSWb == "3SWb"):
    ax3.semilogy(IminusUV[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='UV')
    ax3.semilogy(IminusVIS1[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS1')
    ax3.semilogy(IminusVIS2[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS2')

#ax3.set_xlabel(r"Heating rate, $\mathcal{H} \: [\mathrm{K} \: \mathrm{day}^{-1}]$")
ax3.set_xlabel(r"Downward flux, $\mathcal{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$")

ax3.invert_yaxis()
pl.xscale('symlog')
#ax3.set_xlim(1e-3,1e7)
ax3.set_ylim(pint[-1,step_index,-1],1.)
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
ax3.legend(frameon=True,loc='best',title='Spectral band')

# Remove the y-axis labels
labels = [item.get_text() for item in ax3.get_yticklabels()]

empty_string_labels = ['']*len(labels)
ax3.set_yticklabels(empty_string_labels)

# Remove the x-axis labels
#labels = [item.get_text() for item in ax3.get_xticklabels()]

#empty_string_labels = ['']*len(labels)
#ax3.set_xticklabels(empty_string_labels)

""" ========================== B2 ========================== """

#colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
#for Tsi in TsL_plot_indices:
#    color=next(colors)
    
ax4.semilogy(Direct_beam[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color=color)
ax4.semilogy(Direct_beamIR[Tsi,step_index,:],pint[Tsi,step_index,:],ls=':',color='navy', label='IR')
ax4.semilogy(Direct_beamW1[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-.',color='lightblue', label='W1')
ax4.semilogy(Direct_beamW2[Tsi,step_index,:],pint[Tsi,step_index,:],ls='--',color='royalblue', label='W2')
if (nSWb == "1SWb"):
    ax4.semilogy(Direct_beamSW[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='SW')
if (nSWb == "2SWb"):
    ax4.semilogy(Direct_beamUV[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='UV')
    ax4.semilogy(Direct_beamVIS[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS')
if (nSWb == "3SWb"):
    ax4.semilogy(Direct_beamUV[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='UV')
    ax4.semilogy(Direct_beamVIS1[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS1')
    ax4.semilogy(Direct_beamVIS2[Tsi,step_index,:],pint[Tsi,step_index,:],ls='-',color='g', label='VIS2')

ax4.set_xlabel(r"Direct beam, $\mathcal{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$")
ax4.invert_yaxis()
ax4.set_ylim(pint[-1,step_index,-1],1.)
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)

h, l = ax4.get_legend_handles_labels()
ph = [pl.plot([],marker="", ls="")[0]]*2

handles = ph[:1] + h[0:4] + ph[1:] + h[4:8]
labels = ["Upward fluxes"] + l[0:4] + ["Downward fluxes"] + l[4:8]
leg = pl.legend(handles, labels, ncol=2, loc='best', borderpad=.4)

for vpack in leg._legend_handle_box.get_children():
    for hpack in vpack.get_children()[:1]:
        hpack.get_children()[0].set_width(0)

# Remove the y-axis labels
labels = [item.get_text() for item in ax4.get_yticklabels()]

empty_string_labels = ['']*len(labels)
ax4.set_yticklabels(empty_string_labels)

# Remove the first x-axis label
#labels = [item.get_text() for item in ax4.get_xticklabels()]
#string_labels = ['','',r'$\mathrm{10^{-5}}$',r'$\mathrm{10^{-4}}$',r'$\mathrm{10^{-3}}$',r'$\mathrm{10^{-2}}$',r'$\mathrm{10^{-1}}$',r'$\mathrm{10^{0}}$',r'$\mathrm{10^{1}}$']
#ax4.set_xticklabels(string_labels)

""" ======================= General figure settings ======================= """

#ax2.legend(custom_lines, ['SW', 'IR', 'W1', 'W2'],title='Spectral region')

#FigAbox = AnchoredText(r'\textbf{A}', loc=4, prop={'size':24},frameon=True)
#FigAbox.patch.set(alpha=0.5)
#FigBbox = AnchoredText(r'\textbf{B}', loc=4, prop={'size':24},frameon=True)
#FigBbox.patch.set(alpha=0.5)
#ax1.add_artist(FigAbox)
#ax2.add_artist(FigBbox)

FigA1box = AnchoredText('A1', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
FigA1box.patch.set(alpha=0.5)
FigB1box = AnchoredText('B1', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
FigB1box.patch.set(alpha=0.5)
FigA2box = AnchoredText('A2', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
FigA2box.patch.set(alpha=0.5)
FigB2box = AnchoredText('B2', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
FigB2box.patch.set(alpha=0.5)
ax1.add_artist(FigA1box)
ax2.add_artist(FigB1box)
ax3.add_artist(FigA2box)
ax4.add_artist(FigB2box)

fig.savefig(path_plots+'flux.pdf',bbox_inches='tight')
#pl.show()
pl.close()

""" *********************************************************************** """
""" --------------- Plotting detailed heating rate profiles --------------- """
""" *********************************************************************** """

fig, ((ax1, ax2, ax3)) = pl.subplots(3, sharex=True, sharey=True) #, figsize=(13,6))#(22,10))
fig.subplots_adjust(wspace=0.)
color = 'b'
""" ========================== A1 ========================== """

ax1.semilogy(dT_rad[Tsi,step_index,:],pmid[Tsi,step_index,:], color='r', lw=2, label='dT_rad')
ax1.semilogy(dT_conv[Tsi,step_index,:],pmid[Tsi,step_index,:], color='g', ls='-', label='dT_conv')
ax1.semilogy(dT_conv_moist[Tsi,step_index,:],pmid[Tsi,step_index,:], color='b', ls='-', label='dT_conv_moist')
#ax1.text(2500., 1e6, 'A1', fontsize=14, family='serif')
#ax1.set_xlabel(r"Temperature tendency, $\mathrm{dT} \: [\mathrm{K} \: \mathrm{s}^{-1}]$")
ax1.set_ylabel(r"Pressure, $P$ [Pa]")
ax1.invert_yaxis()
pl.xscale('symlog')
#ax1.set_xlim(-1e-2,1e7)
ax1.set_ylim(pmid[-1,step_index,-1],1.)
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
ax1.legend(frameon=True,loc='best')#,fontsize=14)

""" ========================== A2 ========================== """

ax2.semilogy(dT_rad_planetary[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='b', label='dT_rad_planetary')
ax2.semilogy(dT_rad_stellar[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='r', label='dT_rad_stellar')
#ax2.set_xlabel(r"Temperature tendency, $\mathrm{dT} \: [\mathrm{K} \: \mathrm{s}^{-1}]$")
ax2.set_ylabel(r"Pressure, $P$ [Pa]")
ax2.invert_yaxis()

#ax2.set_xlim(1e-3,1e7)
ax2.set_ylim(pmid[-1,step_index,-1],1.)
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
ax2.legend(frameon=True,loc='best')

""" ========================== A3 ========================== """

pl.xscale('symlog')
ax3.semilogy(dT_rad_IR[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='r', label='IR')
ax3.semilogy(dT_rad_W1[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='lightblue', label='W1')
ax3.semilogy(dT_rad_W2[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='royalblue', label='W2')
ax3.semilogy(dT_rad_UV[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='purple', label='UV')
ax3.semilogy(dT_rad_VIS1[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='orange', label='VIS1')
ax3.semilogy(dT_rad_VIS2[Tsi,step_index,:],pmid[Tsi,step_index,:],ls='-',color='g', label='VIS2')
ax3.set_xlabel(r"Temperature tendency, $\mathrm{dT} \: [\mathrm{K} \: \mathrm{s}^{-1}]$")
ax3.set_ylabel(r"Pressure, $P$ [Pa]")
ax3.invert_yaxis()
pl.xscale('symlog')
#ax3.set_xlim(1e-3,1e7)
ax3.set_ylim(pmid[-1,step_index,-1],1.)
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
ax3.legend(frameon=True,loc='best',bbox_to_anchor=(0.5, 1.2),fontsize=7,title='Spectral band')

""" ======================= General figure settings ======================= """

FigA1box = AnchoredText('A1', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
FigA1box.patch.set(alpha=0.5)
FigA2box = AnchoredText('A2', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
FigA2box.patch.set(alpha=0.5)
FigA3box = AnchoredText('A3', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
FigA3box.patch.set(alpha=0.5)
ax1.add_artist(FigA1box)
ax2.add_artist(FigA2box)
ax3.add_artist(FigA3box)

fig.savefig(path_plots+'heat.pdf',bbox_inches='tight')
#pl.show()
pl.close()

""" *********************************************************************** """
""" ---------------------- Plotting Planck functions ---------------------- """
""" *********************************************************************** """
"""
if not Grey:
    planck=pl.figure(dpi=150)
    pl.semilogy(np.pi*Planck_total[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='Planck_total')
    pl.semilogy(np.pi*Planck_IR[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k:', label='Planck_IR')
    pl.semilogy(np.pi*Planck_W1[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k--', label='Planck_W1')
    pl.semilogy(np.pi*Planck_W2[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k-.', label='Planck_W2')
    if (nSWb == "1SWb"):
        pl.semilogy(np.pi*Planck_SW[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='Planck_SW')
    if (nSWb == "2SWb"):
        pl.semilogy(np.pi*Planck_UV[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='Planck_UV')
        pl.semilogy(np.pi*Planck_VIS[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='Planck_VIS')
    if (nSWb == "3SWb"):
        pl.semilogy(np.pi*Planck_UV[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='Planck_UV')
        pl.semilogy(np.pi*Planck_VIS1[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='Planck_VIS1')
        pl.semilogy(np.pi*Planck_VIS2[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='Planck_VIS2')

    #pl.xlim(-1e-2,1e-4 )
    pl.xlabel(r'Band radiance, $\pi B [W m^{-2}]$')
    pl.ylabel(r'Pressure, $p$ [Pa]')
    pl.tight_layout()
    pl.gca().invert_yaxis()
    pl.legend()
    pl.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
    planck.savefig(path_plots+f'planck{round(TsL[Ts_index])}.pdf',bbox_inches='tight')
    #pl.show()
    pl.close()
"""

tau=pl.figure(dpi=150)
pl.semilogy(tau_total[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='tau_total')
pl.semilogy(tau_IR[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k:', label='tau_IR')
pl.semilogy(tau_W1[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k--', label='tau_W1')
pl.semilogy(tau_W2[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k-.', label='tau_W2')
if (nSWb == "1SWb"):
    pl.semilogy(tau_SW[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='tau_SW')
if (nSWb == "2SWb"):
    pl.semilogy(tau_UV[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='tau_UV')
    pl.semilogy(tau_VIS[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='tau_VIS')
if (nSWb == "3SWb"):
    pl.semilogy(tau_UV[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='tau_UV')
    pl.semilogy(tau_VIS1[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='tau_VIS1')
    pl.semilogy(tau_VIS2[Ts_index,step_index,:],pint[Ts_index,step_index,:],'k', label='tau_VIS2')

#pl.xlim(-1e-2,1e-4 )
pl.xlabel(r'Optical depth, $\tau$')
pl.ylabel(r'Pressure, $p$ [Pa]')
pl.tight_layout()
pl.gca().invert_yaxis()
pl.legend()
pl.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
tau.savefig(path_plots+f'tau{round(TsL[Ts_index])}.pdf',bbox_inches='tight')
#pl.show()
pl.close()

""" *********************************************************************** """
""" ------------- Plotting TOA fluxes vs surface temperatures ------------- """
""" *********************************************************************** """

if len(TsL) > 1:
    OPRTsplot, ax = pl.subplots(dpi=150)
    pl.semilogy(TsL,OPR[:,step_index],'b', label='Planetary Outgoing Radiation (OPR)')
    pl.semilogy(TsL,Direct_beam[:,step_index,0],'r', label='Incoming Shortwave Radiation (ISR)')
    #pl.semilogy(TsL,np.ones(len(flux_down[:,step_index,0]))*187.,'r--',lw=.5,zorder=3)
    #pl.semilogy(TsL,np.ones(len(flux_down[:,step_index,0]))*267.,'r--',lw=.5,zorder=3)
    #pl.semilogy(TsL,np.ones(len(flux_down[:,step_index,0]))*311.,'r--',lw=.5,zorder=3,label='Incoming Stellar Radiation (ISR)')
    pl.semilogy(TsL,OLR[:,step_index],'c', label='Outgoing Longwave Radiation (OLR)')
    pl.semilogy(TsL,OSR[:,step_index],'g', label='Outgoing Shortwave Radiation (OSR)')
    #pl.semilogy(TsL,IplusIR[:,step_index,0],'c:', label='Upward LW IR')
    #pl.semilogy(TsL,IplusW1[:,step_index,0],'c--', label='Upward LW W1')
    #pl.semilogy(TsL,IplusW2[:,step_index,0],'c-.', label='Upward LW W2')
    #literature_comparison(ax) # plot literature comparison (Kopparapu13, Hamano15, Goldblatt13)
    intersection = np.where(abs(OSR[:,step_index]-OLR[:,step_index])<=50.)
    pl.vlines(x=TsL[intersection],ymin=1e1,ymax=OLR[intersection],ls='-',lw=.5,color="k")
    pl.text(800., 40., 'Longwave dominated', fontsize=9, family='serif', color='grey', alpha=0.8)
    pl.text(1800., 40., 'Shortwave dominated', fontsize=9, family='serif', color='grey', alpha=0.8)
    #pl.text(2700., 327., '+6067 Myr', fontsize=7, family='serif', color='r', alpha=0.8)
    #pl.text(2700., 195., '+4567 Myr', fontsize=7, family='serif', color='r', alpha=0.8)
    #pl.text(2700., 130., '+100 Myr', fontsize=7, family='serif', color='r', alpha=0.8)

    pl.text(520., 15., 'Runaway greenhouse', fontsize=9, family='serif', color='grey', alpha=0.8)
    pl.annotate('', fontsize=9, family='serif', xy=(1250., 13.),\
                xycoords='data',\
                    xytext=(510., 13.),\
                        textcoords='data',\
                            arrowprops=dict(arrowstyle= '-|>',\
                                            color='grey',\
                                                lw=2,\
                                                    ls='-')\
                                )
             
             
    pl.text(2220., 15., 'Primordial magma ocean', fontsize=9, family='serif', color='grey', alpha=0.8)
    pl.annotate('', fontsize=9, family='serif', xy=(2150., 13.),\
                xycoords='data',\
                    xytext=(3000., 13.),\
                        textcoords='data',\
                            arrowprops=dict(arrowstyle= '-|>',\
                                            color='grey',\
                                                lw=2,\
                                                    ls='-')\
                                )
             
    literature_comparison(ax)
    pl.xlabel('Surface temperature, $T_\mathrm{s}$ [K]')
    pl.ylabel(r'Top of atmosphere flux, F [W m$^{-2}$]')
    pl.tight_layout()
    #pl.grid()
    pl.legend()
    pl.xlim(500.,TsL[-1])
    pl.ylim(1e1,1e6)
    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
    #pl.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    OPRTsplot.savefig(path_plots+'OPR_Ts.pdf',bbox_inches='tight')
    #pl.show()
    pl.close()


"""" ====================================================================== """
""" ----------------------------- ANIMATIONS ------------------------------ """
"""" ====================================================================== """

if Animations:
    
    """ ******************************************************************* """
    """ ---- Animating temperature profiles and contribution functions ---- """
    """ ******************************************************************* """
    
    if n_files[0] > 1:
        
        if (len(TsL) == 57):   
            
            """ *********************************************************************** """
            """ ------- Plotting temperature profiles and contribution functions ------ """
            """ *********************************************************************** """
            """ Only for an array of 57 Ts values """
                       
            for i in range(n_files[0]):
                
                TsL_plot_indices = [8,26,42,56] # to plot all 4 Ts values on one step plot, all 4 must have the same number of steps                 
                        
                fig, (ax1, ax2) = pl.subplots(1, 2, figsize=(13,6))#(22,10))
                fig.subplots_adjust(wspace=0)
                colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
                for Tsi in TsL_plot_indices:
                    color=next(colors)
                    ax1.semilogy(Tmid[Tsi,i,:],pmid[Tsi,i,:], color=color, label=f'{TsL[Tsi]}')
                    ax1.semilogy(Moist_adiabat[Tsi,:],pmid[Tsi,i,:],'r--')
                    #ax1.text(2500., 1e5, 'A', fontsize=14, family='serif')
                    ax1.set_xlabel(r"Temperature, $T$ [K]")
                    ax1.set_ylabel(r"Pressure, $P$ [Pa]")
                    ax1.invert_yaxis()
                    ax1.set_xlim(0.,3000.)
                    ax1.set_ylim(pmid[-1,i,-1],1.)
                    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                    ax1.legend(frameon=True,loc='best',title=r'$T_\mathrm{surf}$ (K)')#,fontsize=14)
                
                    ax2.loglog(cff[Tsi,i,:],pmid[Tsi,i,:],ls='-',color=color)                       
                    if (scff[Tsi,i] > 0.):
                        ax2.plot(scff[Tsi,i], pmid[Tsi,i,-1], marker='o',color=color, label=f'{Tsi}')               
                
                    ax2.set_xlabel(r"Flux contribution, $\mathcal{CF}_\mathrm{F} \: [\mathrm{W} \: \mathrm{m}^{-2} \mathrm{dex}^{-1}]$")
                    #ax2.text(1e0, 1e6, 'B', fontsize=14, family='serif')
                    #ax2.set_ylabel(r"Pressure, $P$ [Pa]")
                
                    # Minor ticks
                    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12)
                    ax2.yaxis.set_major_locator(locmaj)
                    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
                    ax2.yaxis.set_minor_locator(locmin)
                    ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                
                    ax2.invert_yaxis()
                   
                    ax2.set_ylim(pmid[-1,step_index,-1],1.)
                    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                    ax2.legend(frameon=True,loc='best',title=r'$\mathcal{SCF}_\mathrm{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$')
                
                    #custom_lines = [Line2D([0], [0], color='k', ls='-', lw=1),
                    #                Line2D([0], [0], color='k', ls=':', lw=1),
                    #                Line2D([0], [0], color='k', ls='-.', lw=1),
                    #                Line2D([0], [0], color='k', ls='--', lw=1)]
                
                    #x2labels = ax2.get_xticklabels()
                    #x2labels[0] = "" # remove the first tick label
                
                    # Remove the y-axis labels
                    labels = [item.get_text() for item in ax2.get_yticklabels()]
                
                    empty_string_labels = ['']*len(labels)
                    ax2.set_yticklabels(empty_string_labels)
                
                    # Remove the first x-axis label
                    labels = [item.get_text() for item in ax2.get_xticklabels()]
                    string_labels = ['','',r'$\mathrm{10^{-5}}$',r'$\mathrm{10^{-4}}$',r'$\mathrm{10^{-3}}$',r'$\mathrm{10^{-2}}$',r'$\mathrm{10^{-1}}$',r'$\mathrm{10^{0}}$',r'$\mathrm{10^{1}}$']
                    ax2.set_xticklabels(string_labels)
                
                    #ax2.legend(custom_lines, ['SW', 'IR', 'W1', 'W2'],title='Spectral region')
                
                    """ ======================= General figure settings ======================= """
                
                    ax1.grid(axis='y')
                    ax2.grid(axis='y')
                
                    FigAbox = AnchoredText('A', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                    FigAbox.patch.set(alpha=0.5)
                    FigBbox = AnchoredText('B', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                    FigBbox.patch.set(alpha=0.5)
                    ax1.add_artist(FigAbox)
                    ax2.add_artist(FigBbox)
                
                    fig.savefig(path_plots+'PT_Ts/PT_Ts'+f'_{Tsi}/'+'PT_CFF'+f'_{i}.jpeg',bbox_inches='tight')
                    #pl.show()
                    pl.close()
                
        else:
            
            for Tsi in range(len(TsL)):
                for i in range(n_files[Tsi]):
                
                    fig, (ax1, ax2) = pl.subplots(1, 2, figsize=(13,6))#(22,10))
                    fig.subplots_adjust(wspace=0)
                    ax1.semilogy(Tmid[Tsi,i,:],pmid[Tsi,i,:], color='royalblue', label='Midpoint temperature')
                    ax1.semilogy(Moist_adiabat[Tsi,:],pmid[Tsi,i,:],'r--', label='Dew point')
                    #ax1.text(2500., 1e5, 'A', fontsize=14, family='serif')
                    ax1.set_xlabel(r"Temperature, $T$ [K]")
                    ax1.set_ylabel(r"Pressure, $P$ [Pa]")
                    ax1.invert_yaxis()
                    ax1.set_xlim(0.,3000.)
                    ax1.set_ylim(pmid[-1,step_index,-1],1.)
                    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                    ax1.legend(frameon=True,loc='best')#,fontsize=14)
                                  
                    ax2.loglog(cff[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='royalblue')
                    if (scff[Tsi,i] > 0.):
                        ax2.plot(scff[Tsi,i], pmid[Tsi,i,-1], marker='o',color=color, label=f'{Ts}')
                    
                    ax2.set_xlabel(r"Flux contribution, $\mathcal{CF}_\mathrm{F} \: [\mathrm{W} \: \mathrm{m}^{-2} \mathrm{dex}^{-1}]$")
                    #ax2.text(1e0, 1e6, 'B', fontsize=14, family='serif')
                    #ax2.set_ylabel(r"Pressure, $P$ [Pa]")
                
                    # Minor ticks
                    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12)
                    ax2.yaxis.set_major_locator(locmaj)
                    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
                    ax2.yaxis.set_minor_locator(locmin)
                    ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                
                    ax2.invert_yaxis()
                   
                    ax2.set_ylim(pmid[-1,step_index,-1],1.)
                    sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                    #ax2.legend(frameon=True,loc='best',title=r'$\mathcal{SCF}_\mathrm{F} \: [\mathrm{W} \: \mathrm{m}^{-2} \mathrm{dex}^{-1}]$')
                
                    #custom_lines = [Line2D([0], [0], color='k', ls='-', lw=1),
                    #                Line2D([0], [0], color='k', ls=':', lw=1),
                    #                Line2D([0], [0], color='k', ls='-.', lw=1),
                    #                Line2D([0], [0], color='k', ls='--', lw=1)]
                
                    #x2labels = ax2.get_xticklabels()
                    #x2labels[0] = "" # remove the first tick label
                
                    # Remove the y-axis labels
                    labels = [item.get_text() for item in ax2.get_yticklabels()]
                
                    empty_string_labels = ['']*len(labels)
                    ax2.set_yticklabels(empty_string_labels)
                
                    # Remove the first x-axis label
                    labels = [item.get_text() for item in ax2.get_xticklabels()]
                    string_labels = ['','',r'$\mathrm{10^{-5}}$',r'$\mathrm{10^{-4}}$',r'$\mathrm{10^{-3}}$',r'$\mathrm{10^{-2}}$',r'$\mathrm{10^{-1}}$',r'$\mathrm{10^{0}}$',r'$\mathrm{10^{1}}$']
                    ax2.set_xticklabels(string_labels)
                
                    #ax2.legend(custom_lines, ['SW', 'IR', 'W1', 'W2'],title='Spectral region')
                
                    """ ======================= General figure settings ======================= """
                
                    ax1.grid(axis='y')
                    ax2.grid(axis='y')
                
                    FigAbox = AnchoredText('A', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                    FigAbox.patch.set(alpha=0.5)
                    FigBbox = AnchoredText('B', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                    FigBbox.patch.set(alpha=0.5)
                    ax1.add_artist(FigAbox)
                    ax2.add_artist(FigBbox)
                
                    fig.savefig(path_plots+'PT_Ts/PT_Ts'+f'_{Tsi:02}/'+'PT_CFF'+f'_{i}.jpeg',bbox_inches='tight')
                    #pl.show()
                    pl.close()
            
        for Tsi in range(len(TsL)):
            for i in range(n_files[Tsi]):
                    
                """ *********************************************************************** """
                """ --------------- Plotting flux and heating rate profiles --------------- """
                """ *********************************************************************** """
                    
                fig, ((ax1, ax3), (ax2, ax4)) = pl.subplots(2, 2, figsize=(13,6))#(22,10))
                fig.subplots_adjust(wspace=0.)
                color = 'b'
                """ ========================== A1 ========================== """
                    
                #colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
                #for Tsi in TsL_plot_indices:
                #    color=next(colors)
                ax1.semilogy(diffuseflux[Tsi,i,:],pint[Tsi,i,:], color='b', lw=2, label='Net flux')
                ax1.semilogy(Iplus[Tsi,i,:],pint[Tsi,i,:], color='c', ls='-', label='Upward diffuse flux')
                ax1.semilogy(Iminus[Tsi,i,:],pint[Tsi,i,:], color='orange', ls='-', label='Downward diffuse flux')
                ax1.semilogy(Direct_beam[Tsi,i,:],pint[Tsi,i,:], color='r', ls='-', label='Direct beam')
                #ax1.text(2500., 1e6, 'A1', fontsize=14, family='serif')
                ax1.set_xlabel(r"Diffuse net flux, $\mathcal{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$")
                ax1.set_ylabel(r"Pressure, $P$ [Pa]")
                ax1.invert_yaxis()
                pl.xscale('symlog')
                #ax1.set_xlim(-1e-2,1e7)
                ax1.set_ylim(pint[-1,step_index,-1],1.)
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                ax1.legend(frameon=True,loc='best')#,fontsize=14)
                    
                """ ========================== B1 ========================== """
                    
                #colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
                #for Tsi in TsL_plot_indices:
                #    color=next(colors)
                ax2.semilogy(IplusIR[Tsi,i,:],pint[Tsi,i,:],ls='-.',color='lightblue', label='W1')
                ax2.semilogy(IplusW1[Tsi,i,:],pint[Tsi,i,:],ls='--',color='royalblue', label='W2')
                ax2.semilogy(IplusW2[Tsi,i,:],pint[Tsi,i,:],ls=':',color='navy', label='IR')
                if (nSWb == "1SWb"):
                    ax2.semilogy(IplusSW[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='SW')
                if (nSWb == "2SWb"):
                    ax2.semilogy(IplusUV[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='UV')
                    ax2.semilogy(IplusVIS[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS')
                if (nSWb == "3SWb"):
                    ax2.semilogy(IplusUV[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='UV')
                    ax2.semilogy(IplusVIS1[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS1')
                    ax2.semilogy(IplusVIS2[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS2')
                    
                ax2.set_xlabel(r"Upward flux, $\mathcal{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$")
                ax2.set_ylabel(r"Pressure, $P$ [Pa]")
                ax2.invert_yaxis()
                    
                #ax2.set_xlim(1e-3,1e7)
                ax2.set_ylim(pint[-1,step_index,-1],1.)
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                ax2.legend(frameon=True,loc='best',title='Spectral band')
                    
                """ ========================== A2 ========================== """
                    
                #colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
                #for Tsi in TsL_plot_indices:
                #    color=next(colors)
                pl.xscale('symlog')
                ax3.semilogy(Iminus[Tsi,i,:],pint[Tsi,i,:],ls='-',color=color)
                ax3.semilogy(IminusIR[Tsi,i,:],pint[Tsi,i,:],ls=':',color='navy', label='IR')
                ax3.semilogy(IminusW1[Tsi,i,:],pint[Tsi,i,:],ls='-.',color='lightblue', label='W1')
                ax3.semilogy(IminusW2[Tsi,i,:],pint[Tsi,i,:],ls='--',color='royalblue', label='W2')
                if (nSWb == "1SWb"):
                    ax3.semilogy(IminusSW[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='SW')
                if (nSWb == "2SWb"):
                    ax3.semilogy(IminusUV[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='UV')
                    ax3.semilogy(IminusVIS[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS')
                if (nSWb == "3SWb"):
                    ax3.semilogy(IminusUV[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='UV')
                    ax3.semilogy(IminusVIS1[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS1')
                    ax3.semilogy(IminusVIS2[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS2')
                    
                #ax3.set_xlabel(r"Heating rate, $\mathcal{H} \: [\mathrm{K} \: \mathrm{day}^{-1}]$")
                ax3.set_xlabel(r"Downward flux, $\mathcal{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$")
                    
                ax3.invert_yaxis()
                pl.xscale('symlog')
                #ax3.set_xlim(1e-3,1e7)
                ax3.set_ylim(pint[-1,step_index,-1],1.)
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                ax3.legend(frameon=True,loc='best',title='Spectral band')
                    
                # Remove the y-axis labels
                labels = [item.get_text() for item in ax3.get_yticklabels()]
                    
                empty_string_labels = ['']*len(labels)
                ax3.set_yticklabels(empty_string_labels)
                    
                # Remove the x-axis labels
                #labels = [item.get_text() for item in ax3.get_xticklabels()]
                    
                #empty_string_labels = ['']*len(labels)
                #ax3.set_xticklabels(empty_string_labels)
                    
                """ ========================== B2 ========================== """
                    
                #colors=iter(pl.cm.Blues(np.linspace(0.3,1,len(TsL_plot_indices))))
                #for Tsi in TsL_plot_indices:
                #    color=next(colors)
                    
                ax4.semilogy(Direct_beam[Tsi,i,:],pint[Tsi,i,:],ls='-',color=color)
                ax4.semilogy(Direct_beamIR[Tsi,i,:],pint[Tsi,i,:],ls=':',color='navy', label='IR')
                ax4.semilogy(Direct_beamW1[Tsi,i,:],pint[Tsi,i,:],ls='-.',color='lightblue', label='W1')
                ax4.semilogy(Direct_beamW2[Tsi,i,:],pint[Tsi,i,:],ls='--',color='royalblue', label='W2')
                if (nSWb == "1SWb"):
                    ax4.semilogy(Direct_beamSW[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='SW')
                if (nSWb == "2SWb"):
                    ax4.semilogy(Direct_beamUV[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='UV')
                    ax4.semilogy(Direct_beamVIS[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS')
                if (nSWb == "3SWb"):
                    ax4.semilogy(Direct_beamUV[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='UV')
                    ax4.semilogy(Direct_beamVIS1[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS1')
                    ax4.semilogy(Direct_beamVIS2[Tsi,i,:],pint[Tsi,i,:],ls='-',color='g', label='VIS2')
                    
                ax4.set_xlabel(r"Direct beam, $\mathcal{F} \: [\mathrm{W} \: \mathrm{m}^{-2}]$")
                ax4.invert_yaxis()
                ax4.set_ylim(pint[-1,step_index,-1],1.)
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                    
                h, l = ax4.get_legend_handles_labels()
                ph = [pl.plot([],marker="", ls="")[0]]*2
                    
                handles = ph[:1] + h[0:4] + ph[1:] + h[4:8]
                labels = ["Upward fluxes"] + l[0:4] + ["Downward fluxes"] + l[4:8]
                leg = pl.legend(handles, labels, ncol=2, loc='best', borderpad=.4)
                    
                for vpack in leg._legend_handle_box.get_children():
                    for hpack in vpack.get_children()[:1]:
                        hpack.get_children()[0].set_width(0)
                    
                # Remove the y-axis labels
                labels = [item.get_text() for item in ax4.get_yticklabels()]
                    
                empty_string_labels = ['']*len(labels)
                ax4.set_yticklabels(empty_string_labels)
                    
                # Remove the first x-axis label
                #labels = [item.get_text() for item in ax4.get_xticklabels()]
                #string_labels = ['','',r'$\mathrm{10^{-5}}$',r'$\mathrm{10^{-4}}$',r'$\mathrm{10^{-3}}$',r'$\mathrm{10^{-2}}$',r'$\mathrm{10^{-1}}$',r'$\mathrm{10^{0}}$',r'$\mathrm{10^{1}}$']
                #ax4.set_xticklabels(string_labels)
                    
                """ ======================= General figure settings ======================= """
                    
                #ax2.legend(custom_lines, ['SW', 'IR', 'W1', 'W2'],title='Spectral region')
                    
                #FigAbox = AnchoredText(r'\textbf{A}', loc=4, prop={'size':24},frameon=True)
                #FigAbox.patch.set(alpha=0.5)
                #FigBbox = AnchoredText(r'\textbf{B}', loc=4, prop={'size':24},frameon=True)
                #FigBbox.patch.set(alpha=0.5)
                #ax1.add_artist(FigAbox)
                #ax2.add_artist(FigBbox)
                    
                FigA1box = AnchoredText('A1', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                FigA1box.patch.set(alpha=0.5)
                FigB1box = AnchoredText('B1', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                FigB1box.patch.set(alpha=0.5)
                FigA2box = AnchoredText('A2', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                FigA2box.patch.set(alpha=0.5)
                FigB2box = AnchoredText('B2', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                FigB2box.patch.set(alpha=0.5)
                ax1.add_artist(FigA1box)
                ax2.add_artist(FigB1box)
                ax3.add_artist(FigA2box)
                ax4.add_artist(FigB2box)
                    
                fig.savefig(path_plots+'flux_Ts/flux_Ts'+f'_{Tsi:02}/'+'flux'+f'_{i}.jpeg',bbox_inches='tight')
                #pl.show()
                pl.close()
            
        for Tsi in range(len(TsL)):
            for i in range(n_files[Tsi]):

                """ *********************************************************************** """
                """ --------------- Plotting detailed heating rate profiles --------------- """
                """ *********************************************************************** """

                fig, ((ax1, ax2, ax3)) = pl.subplots(3, sharex=True, sharey=True)#, figsize=(13,6))#(22,10))
                fig.subplots_adjust(wspace=0.)
                color = 'b'
                """ ========================== A1 ========================== """

                ax1.semilogy(dT_rad[Tsi,i,:],pmid[Tsi,i,:], color='r', lw=2, label='dT_rad')
                ax1.semilogy(dT_conv[Tsi,i,:],pmid[Tsi,i,:], color='g', ls='-', label='dT_conv')
                ax1.semilogy(dT_conv_moist[Tsi,i,:],pmid[Tsi,i,:], color='b', ls='-', label='dT_conv_moist')
                #ax1.text(2500., 1e6, 'A1', fontsize=14, family='serif')
                #ax1.set_xlabel(r"Temperature tendency, $\mathcal{dT} \: [\mathrm{K} \: \mathrm{s}^{-1}]$")
                ax1.set_ylabel(r"Pressure, $P$ [Pa]")
                ax1.invert_yaxis()
                pl.xscale('symlog')
                #ax1.set_xlim(-1e-2,1e7)
                ax1.set_ylim(pint[-1,step_index,-1],1.)
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                ax1.legend(frameon=True,loc='best')#,fontsize=14)

                """ ========================== A2 ========================== """

                ax2.semilogy(dT_rad_planetary[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='b', label='dT_rad_planetary')
                ax2.semilogy(dT_rad_stellar[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='r', label='dT_rad_stellar')
                #ax2.set_xlabel(r"Temperature tendency, $\mathcal{dT} \: [\mathrm{K} \: \mathrm{s}^{-1}]$")
                ax2.set_ylabel(r"Pressure, $P$ [Pa]")
                ax2.invert_yaxis()

                #ax2.set_xlim(1e-3,1e7)
                ax2.set_ylim(pint[-1,step_index,-1],1.)
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                ax2.legend(frameon=True,loc='best',title='Spectral band')

                """ ========================== A3 ========================== """

                pl.xscale('symlog')
                ax3.semilogy(dT_rad_IR[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='r', label='IR')
                ax3.semilogy(dT_rad_W1[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='lightblue', label='W1')
                ax3.semilogy(dT_rad_W2[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='royalblue', label='W2')
                ax3.semilogy(dT_rad_UV[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='purple', label='UV')
                ax3.semilogy(dT_rad_VIS1[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='orange', label='VIS1')
                ax3.semilogy(dT_rad_VIS2[Tsi,i,:],pmid[Tsi,i,:],ls='-',color='g', label='VIS2')
                ax3.set_xlabel(r"Temperature tendency, $\mathrm{dT} \: [\mathrm{K} \: \mathrm{s}^{-1}]$")
                ax3.set_ylabel(r"Pressure, $P$ [Pa]")
                ax3.invert_yaxis()
                pl.xscale('symlog')
                #ax3.set_xlim(1e-3,1e7)
                ax3.set_ylim(pint[-1,step_index,-1],1.)
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                ax3.legend(frameon=True,loc='best',bbox_to_anchor=(0.5, 1.2),fontsize=7,title='Spectral band')

                """ ======================= General figure settings ======================= """

                FigA1box = AnchoredText('A1', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                FigA1box.patch.set(alpha=0.5)
                FigA2box = AnchoredText('A2', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                FigA2box.patch.set(alpha=0.5)
                FigA3box = AnchoredText('A3', loc=4, prop={'size':16, 'family':'serif'},frameon=True)
                FigA3box.patch.set(alpha=0.5)
                ax1.add_artist(FigA1box)
                ax2.add_artist(FigA2box)
                ax3.add_artist(FigA3box)

                fig.savefig(path_plots+'heat_Ts/heat_Ts'+f'_{Tsi:02}/'+'heat'+f'_{i}.jpeg',bbox_inches='tight')
                 #pl.show()
                pl.close()

        for Tsi in range(len(TsL)):
            for i in range(n_files[Tsi]):
                    
                tau=pl.figure(dpi=150)
                pl.semilogy(tau_total[Tsi,i,:],pint[Tsi,i,:],'k', label='tau_total')
                pl.semilogy(tau_IR[Tsi,i,:],pint[Tsi,i,:],'k:', label='tau_IR')
                pl.semilogy(tau_W1[Tsi,i,:],pint[Tsi,i,:],'k--', label='tau_W1')
                pl.semilogy(tau_W2[Tsi,i,:],pint[Tsi,i,:],'k-.', label='tau_W2')
                if (nSWb == "1SWb"):
                    pl.semilogy(tau_SW[Tsi,i,:],pint[Tsi,i,:],'k', label='tau_SW')
                if (nSWb == "2SWb"):
                    pl.semilogy(tau_UV[Tsi,i,:],pint[Tsi,i,:],'k', label='tau_UV')
                    pl.semilogy(tau_VIS[Tsi,i,:],pint[Tsi,i,:],'k', label='tau_VIS')
                if (nSWb == "3SWb"):
                    pl.semilogy(tau_UV[Tsi,i,:],pint[Tsi,i,:],'k', label='tau_UV')
                    pl.semilogy(tau_VIS1[Tsi,i,:],pint[Tsi,i,:],'k', label='tau_VIS1')
                    pl.semilogy(tau_VIS2[Tsi,i,:],pint[Tsi,i,:],'k', label='tau_VIS2')
                
                
                #pl.xlim(-1e-2,1e-4 )
                pl.xlabel(r'Optical depth, $\tau$')
                pl.ylabel(r'Pressure, $p$ [Pa]')
                pl.tight_layout()
                pl.gca().invert_yaxis()
                pl.legend()
                pl.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                tau.savefig(path_plots+'tau_Ts/tau_Ts'+f'_{Tsi:02}/'+'tau'+f'_{i}.jpeg',bbox_inches='tight')
                    
                #pl.show()
                pl.close()
            
        """ *********************************************************************** """
        """ ------------- Plotting TOA fluxes vs surface temperatures ------------- """
        """ *********************************************************************** """
            
        if len(TsL) > 1:
            for i in range(n_files[Tsi]):
                                
                OPRTsplot, ax = pl.subplots(dpi=150)
                pl.semilogy(TsL,OPR[:,i],'b', label='Planetary Outgoing Radiation (OPR)')
                pl.semilogy(TsL,Direct_beam[:,i,0],'r', label='Incoming Shortwave Radiation (ISR)')
                #pl.semilogy(TsL,np.ones(len(flux_down[:,i,0]))*187.,'r--',lw=.5,zorder=3)
                #pl.semilogy(TsL,np.ones(len(flux_down[:,i,0]))*267.,'r--',lw=.5,zorder=3)
                #pl.semilogy(TsL,np.ones(len(flux_down[:,i,0]))*311.,'r--',lw=.5,zorder=3,label='Incoming Stellar Radiation (ISR)')
                pl.semilogy(TsL,OLR[:,i],'c', label='Outgoing Longwave Radiation (OLR)')
                pl.semilogy(TsL,OSR[:,i],'g', label='Outgoing Shortwave Radiation (OSR)')
                #pl.semilogy(TsL,IplusIR[:,i,0],'c:', label='Upward LW IR')
                #pl.semilogy(TsL,IplusW1[:,i,0],'c--', label='Upward LW W1')
                #pl.semilogy(TsL,IplusW2[:,i,0],'c-.', label='Upward LW W2')
                #literature_comparison(ax) # plot literature comparison (Kopparapu13, Hamano15, Goldblatt13)
                intersection = np.where(abs(OSR[:,i]-OLR[:,i])<=50.)
                pl.vlines(x=TsL[intersection],ymin=1e1,ymax=OLR[intersection],ls='-',lw=.5,color="k")
                pl.text(800., 40., 'Longwave dominated', fontsize=9, family='serif', color='grey', alpha=0.8)
                pl.text(1800., 40., 'Shortwave dominated', fontsize=9, family='serif', color='grey', alpha=0.8)
                pl.text(2700., 327., '+6067 Myr', fontsize=7, family='serif', color='r', alpha=0.8)
                pl.text(2700., 195., '+4567 Myr', fontsize=7, family='serif', color='r', alpha=0.8)
                pl.text(2700., 130., '+100 Myr', fontsize=7, family='serif', color='r', alpha=0.8)
                
                pl.text(520., 15., 'Runaway greenhouse', fontsize=9, family='serif', color='grey', alpha=0.8)
                pl.annotate('', fontsize=9, family='serif', xy=(1250., 13.),\
                            xycoords='data',\
                                xytext=(510., 13.),\
                                    textcoords='data',\
                                        arrowprops=dict(arrowstyle= '-|>',\
                                                        color='grey',\
                                                            lw=2,\
                                                                ls='-')\
                                            )
                             
                             
                pl.text(2220., 15., 'Primordial magma ocean', fontsize=9, family='serif', color='grey', alpha=0.8)
                pl.annotate('', fontsize=9, family='serif', xy=(2150., 13.),\
                            xycoords='data',\
                                xytext=(3000., 13.),\
                                    textcoords='data',\
                                        arrowprops=dict(arrowstyle= '-|>',\
                                                        color='grey',\
                                                            lw=2,\
                                                                ls='-')\
                                            )
                             
                literature_comparison(ax)
                pl.xlabel('Surface temperature, $T_\mathrm{s}$ [K]')
                pl.ylabel(r'Top of atmosphere flux, F [W m$^{-2}$]')
                pl.tight_layout()
                #pl.grid()
                pl.legend()
                pl.xlim(500.,TsL[-1])
                pl.ylim(1e1,np.max(OPR[:,step_index]))
                sns.despine(fig=None, ax=None, top=True, right=True, left=False, bottom=False, offset=None, trim=False)
                #pl.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                OPRTsplot.savefig(path_plots+'OPR_Ts_'+f'{Tsi:02}/'+'OPR_Ts'+f'_{i}.jpeg',bbox_inches='tight')
                #pl.show()
                pl.close()
