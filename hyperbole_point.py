# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:36:27 2021

@author: ftholin
"""


import sys,time
import os.path

import subprocess
from shutil import copyfile


import numpy as np

import matplotlib
import matplotlib.pyplot as plt


#################################################
#
# Field and potential due to an hyperbole point
# cf. Celestin 2008 eq. IV7-IV14
#
#################################################

def ksi2(a,x,r):
    return (a**2+x**2+r**2-((a**2+x**2+r**2)**2-4*a**2*x**2)**0.5)/(2*a**2)

def ksi(a,x,r):
    return np.sqrt(ksi2(a,x,r))

def eta(a,x,r):
    return r/(a*(1-ksi2(a,x,r)**2)**0.5)

def psi(a,x,r):
    return a*((eta(a,x,r)**2-ksi2(a,x,r)+1)/(1-ksi2(a,x,r)**2))**0.5

def Acst(V,gap,rad):
    ksi0=gap/(gap+rad)
    return 2.0*V/np.log((1.0+ksi0)/(1.0-ksi0))

def potential(V,rad,gap,x,r):
    alpha=rad+gap
    return 0.5*Acst(V,gap,rad)*np.log((1+ksi(alpha,x,r))/(1-ksi(alpha,x,r)))+V
                 
def electric_field(V,rad,gap,x,r):
    alpha=rad+gap
    return 1/psi(alpha,x,r)*Acst(V,gap,rad)/(1.0-ksi2(alpha,x,r))

#Anode potential:
Va=55000

### gap between tip and summetry plane
gap=1.6e-2

### curvature radius at the tip
rad=500.0e-6

### r
r=0.0
Nx=200
x=np.zeros(Nx)
V=np.zeros(Nx)
E=np.zeros(Nx)

for i in range(Nx):
    x[i]=i*gap/Nx    
    V[i]=potential(Va,rad,gap,x[i],r)
    E[i]=electric_field(Va,rad,gap,x[i],r)

fig, ax = plt.subplots(figsize=(6, 3.75))#plt.xlim(0, 0.03)

##===================== FONTS ========================= 

plt.rc('font',   size      = 12) # controls default text sizes
plt.rc('axes',   titlesize = 12) # fontsize of the axes title
plt.rc('axes',   labelsize = 12) # fontsize of the x and y labels
plt.rc('xtick',  labelsize = 12) # fontsize of the xtick labels
plt.rc('ytick',  labelsize = 12) # fontsize of the ytick labels
plt.rc('legend', fontsize  = 12) # legend fontsize
plt.rc('figure', titlesize = 12) # fontsize of the figure title

#=== label
ax.set_xlabel('x [ mm ]')
ax.set_ylabel('$E$ [ V/m ]', color='blue')
ax2 = ax.twinx()
ax2.set_ylabel('$V$ [ V ]', color = 'red') 

#=== plots:
line1,=ax2.plot(1.0e3*x, V, dashes=[2,2],  lw=2, c='red', alpha=1,zorder=2, label='$V_\mathrm{analytic}$')
line2,=ax.plot(1.0e3*x, E, dashes=[3,2],  lw=2, c='blue', alpha=1,zorder=2, label='$E_\mathrm{analytic}$')

plt.legend(handles=[line1, line2],loc='upper center', frameon=False,
          bbox_to_anchor=(0.5, 1.2), shadow=False, ncol=2)

line1.remove()
line2.remove()

#### En fonction du rayon de courbure

Nrad=10
radii=np.zeros(Nrad)
Erad=np.zeros((Nx,Nrad))



for j in range(Nrad):  
    radii[j]=100.0e-6+j*100.0e-6
    for i in range(Nx):    
        Erad[i,j]=electric_field(Va,radii[j],gap,x[i],r)

ax2.remove()
ax.set_ylabel('$E$ [ V/m ]', color='black')
ax2 = ax.twinx()
ax2.set_ylabel('$ratio$', color = 'grey') 

j=0        
line1,=ax.plot(1.0e3*x, Erad[:,j], dashes=[3,2],  lw=2, c='blue', alpha=1,zorder=2, label='$E_1\,(R={:4d}$'.format(int(1e6*radii[j]))+'$\,\mu\mathrm{m})$')

j=9        
line2,=ax.plot(1.0e3*x, Erad[:,j], dashes=[1,2],  lw=2, c='red', alpha=1,zorder=2, label='$E_2\,(R={:4d}$'.format(int(1e6*radii[j]))+'$\,\mu\mathrm{m})$')

#line3,=ax2.plot(1.0e3*x, 0.5*(Erad[:,0]-Erad[:,9])/(Erad[:,9]+Erad[:,0]),  lw=4, c='grey', alpha=0.3,zorder=2, label='$E_2$/$E_1$')

line3,=ax2.plot(1.0e3*x, Erad[:,9]/Erad[:,0], lw=4, c='grey', alpha=0.3,zorder=2, label='$E_2$/$E_1$')

plt.legend(handles=[line1, line2, line3],loc='upper center', frameon=False, bbox_to_anchor=(0.5, 1.2), shadow=False, ncol=3)




#line2,=ax.plot(1.0e3*x, Erad[:,1], dashes=[3,2],  lw=2, c='blue', alpha=1,zorder=2, label='$E$')
#line3,=ax.plot(1.0e3*x, Erad[:,2], dashes=[3,2],  lw=2, c='red', alpha=1,zorder=2, label='$E$')
#line4,=ax.plot(1.0e3*x, Erad[:,3], dashes=[3,2],  lw=2, c='green', alpha=1,zorder=2, label='$E$')
#line5,=ax.plot(1.0e3*x, Erad[:,4], dashes=[3,2],  lw=2, c='orange', alpha=1,zorder=2, label='$E$')


    