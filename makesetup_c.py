# -*- coding: utf-8 -*-
import numpy as np
from namelist import *
from meteo_utilities import rrmixv1
import sys

def maketopo(topo,nxb):
    """
    Topography definition

    Input:  maketopo(topo,nxb)
    Output: topo
    """
    if idbg == 1:
        print('Topography ...\n')

    x = np.arange(0,nxb,dtype='float32')

    x0 = (nxb - 1)/2. + 1
    x = (x+1 - x0)*dx
    

    
    
    toponf = topomx*np.exp(-(x/float(topowd))**2)
    

#    topo[1:-1,0] =  toponf[1:-1] + 0.25*(toponf[0:-2] - 2.*toponf[1:-1] +
#                                        toponf[2:]) 
    topo[:,0] =  - toponf[:] + 0.25*(toponf[:] - 2.*toponf[:] + 
                                        toponf[:]) + 1500
    return topo


def makeprofile(sold,snow,uold,unow,mtg,mtgnew,qvold=0,qvnow=0,
                qcold=0,qcnow=0,qrold=0,qrnow=0,ncold=0,ncnow=0,nrold=0,nrnow=0):
    """
    Make upstream profiles and initial conditions for
    isentropic density (sigma) and velocity (u)

    Input:      makeprofile(sold,snow,uold,unow,mtg,mtgnew)
    Output:     th0, exn0, prs0, z0, mtg0, s0, u0, sold, ...
                snow, uold, unow, mtg, mtgnew
    """
    #global dth
    if idbg == 1:
        print('Create initial profile ...\n')

    exn0 = np.zeros(nz1)
    z0   = np.zeros(nz1)
    mtg0 = np.zeros(nz)
    prs0 = np.zeros(nz1)
    exn0 = np.zeros(nz1)
    rh0 = np.zeros(nz)
    qv0 = np.zeros(nz)

    if imoist==1:
        qc0 = np.zeros(nz)
        qr0 = np.zeros(nz)
    if imoist==1 and imicrophys==2:
        nc0 = np.zeros(nz)
        nr0 = np.zeros(nz)

    # Upstream profile for Brunt-Vaisalla frequency (unstaggered)
    #------------------------------------------------------------
    bv0  = bv00*np.ones(nz1)

    # Upstream profile of theta (staggered)
    # -----------------------------------------------------------
    th0 = th00*np.ones(nz1) + dth*np.arange(0,nz1)

    # Upstream profile for Exner function and pressure (staggered)
    #-------------------------------------------------------------
    exn0[0] = exn00
    for k in range(1,nz1):
        exn0[k] = exn0[k-1] - 16*(g**2)*(th0[k] - th0[k-1])/((bv0[k-1] +
                  bv0[k])**2 *(th0[k-1] + th0[k])**2)

    prs0[:] = pref*(exn0[:]/cp)**cpdr

    # Upstream profile for geometric height (staggered)
    #-------------------------------------------------------------
    z0[0] = z00
    for k in range(1,nz1):
        z0[k] = z0[k-1] + 8*g*(th0[k] -th0[k-1])/((th0[k-1] +
                th0[k])*(bv0[k-1] + bv0[k])**2)

    # Upstream profile for Montgomery potential (unstaggered)
    #--------------------------------------------------------
    mtg0[0] = g*z0[0] + th00*exn0[0] + dth*exn0[0]/2.
    mtg0[1:nz] = mtg0[0:nz-1] + dth*exn0[1:nz]

    # Upstream profile for isentropic density (unstaggered)
    #------------------------------------------------------
    s0 = -1./g*(prs0[1:] - prs0[0:-1])/float(dth)

    # Upstream profile for velocity (unstaggered)
    #--------------------------------------------
    u0   = float(u00)*np.ones(nz)


    if ishear == 1:
        if idbg == 1:
            print('Using wind shear profile ...\n')
        # *** Exercise 3.3 Downslope windstorm ***
        # *** use indices k_shl, k_sht, and wind speeds u00_sh, u00
        #
        u_shear_zone = np.linspace(u00_sh,u00,k_sht - k_shl)
        # *** edit here ...
        u0[0:k_shl] = u00_sh
        u0[k_shl:k_sht] = u_shear_zone
        u0[k_sht::] = u00
    
        # *** Exercise 3.3 Downslope windstorm ***
    else:
        if idbg == 1:
            print('Using uniform wind profile ...\n')


    # Upstream profile for moisture (unstaggered)
    # -------------------------------------------


    if imoist ==1:

        # *** Exercise 4.1 Initial Moisture profile ***
        # *** define new indices and create the profile ***
        # *** for rh0; then use function rrmixv1 to compute qv0 ***
        

        # *** edit here ...
        kw = 10
        kc = 11
        rhmax = 0.98
        for k in np.arange(kc - kw + 1, kc + kw, 1):
            
            rh0[k] = rhmax * (np.cos(np.abs(k - kc) / kw * np.pi / 2))**2
            
            qv0[k] = rrmixv1(0.5*(prs0[k] + prs0[k+1])/100, 
               0.5*(th0[k]/cp*exn0[k] + th0[k+1]/cp*exn0[k+1]), rh0[k], 2)

        # *** Exercise 4.1 Initial Moisture profile ***

        # Upstream profile for number densities (unstaggered)
        # -------------------------------------------
        if imicrophys==2:
            nc0 = np.zeros(nz)
            nr0 = np.zeros(nz)

    # Initial conditions for isentropic density (sigma), velocity u, and moisture qv
    # ---------------------------------------------------------------------

    sold =   s0*np.ones_like(sold,dtype = np.float)
    snow =   s0*np.ones_like(sold,dtype = np.float)
    mtg =    mtg0*np.ones_like(sold,dtype = np.float)
    mtgnew = mtg0*np.ones_like(sold,dtype = np.float)
    uold =   u0*np.ones_like(uold,dtype = np.float)
    unow =   u0*np.ones_like(uold,dtype = np.float)

    if imoist==1:
        #if imicrophys!=0: 
        qvold = qv0*np.ones_like(qvold,dtype = np.float)
        qvnow = qv0*np.ones_like(qvold,dtype = np.float)
        qcold = qc0*np.ones_like(qcold,dtype = np.float)
        qcnow = qc0*np.ones_like(qcold,dtype = np.float)
        qrold = qr0*np.ones_like(qrold,dtype = np.float)
        qrnow = qr0*np.ones_like(qrold,dtype = np.float)

        # droplet density for 2-moment scheme
        if imicrophys==2:
            ncold = nc0*np.ones_like(ncold,dtype = np.float)
            ncnow = nc0*np.ones_like(ncold,dtype = np.float)
            nrold = nr0*np.ones_like(nrold,dtype = np.float)
            nrnow = nr0*np.ones_like(nrold,dtype = np.float)

    if imoist == 0:
        return th0,exn0,prs0,z0,mtg0,s0,u0,sold,snow,uold,unow,mtg,mtgnew
    else:
        if imicrophys == 0 or imicrophys==1:
            return th0,exn0,prs0,z0,mtg0,s0,u0,sold,snow,uold,unow,mtg, \
                   mtgnew,qv0,qc0,qr0,qvold,qvnow,qcold,qcnow,qrold,qrnow
        elif imicrophys == 2:
            return th0,exn0,prs0,z0,mtg0,s0,u0,sold,snow,uold,unow,mtg, \
                   mtgnew,qv0,qc0,qr0,qvold,qvnow,qcold,qcnow,qrold,qrnow, \
                   ncold,ncnow,nrold,nrnow

# END OF MAKESETUP.PY
