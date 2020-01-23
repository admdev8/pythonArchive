#!/usr/bin/env python 
#
# This script contains methods which were used to gnerate the fermi energy 
# given various models of the density of interface states (Dit) in graphene 
# FET devices. Several models Dit density are included for computational
# experimentation 
#
# Hysteresis modeling in graphene field effect transistors
#
# Journal of Applied Physics 117, 074501 (2015); https://doi.org/10.1063/1.4913209
#
import math as m
import numpy as np
import matplotlib.pyplot as plt

#############################
# CARRIER DENSITY VARIABLES #
#############################
kb   = 8.617e-5  #eV/K
hbar = 6.58e-16  #eV*s
vf   = 1.0e8     #cm/s
g    = 0.400     #eV
T    = 350.0     #K
gsgv = 4.        #none 
e    = 1.602e-19 #C

# Capacitance calculations
er   = 6.0
ep   = e*8.85e-14  #F/cm with J->eV  
t    = 13.0e-7     #cm 
Dit  = 0           #cm^{-2}


##################################
#        Constant D_it           #
##################################
def fermiWave(Vmax, vd = 0, n = 1): 

    eta  = ((gsgv)/(4*m.pi))/((hbar*vf)**2)
    Cox  = (er*ep)/t
  
    # Calculations
    C = (e**2)/Cox #eV^{-1}
    a = (eta*C/2.0)
    b = C*((eta*g/2)+Dit)

    ef = lambda dv: (1/(2*a))*(m.sqrt(4*a*dv + b**2) - b)
    _v = np.asarray([Vmax*m.cos(i)-vd for i in np.linspace(0,2*n*m.pi,50)])
    v  = np.asarray([Vmax*m.cos(i) for i in np.linspace(0,2*n*m.pi,50)])

    f = []
    for i in _v:
        if i >= 0: 
            f.append( ef(i) ) 
        else: 
            f.append( -1*ef(abs(i)))
    return v,f,v

##################################
#        Gaussian D_it           #
##################################
def erf(z):
    if z == 0:
        return 0
    else:
        return ((2*z*np.exp(-z**2))/m.sqrt(m.pi)) * (1./_erf(z))
def _erf(z, n = 0):
    if n == 50:
        return 2.*z
    else:
        if n%2 == 0:
            return (1 + 2.*n) - (2.*(n+1)*z**2 )/_erf(z, n+1) 
        else:
            return (1 + 2.*n) + (2.*(n+1)*z**2 )/_erf(z, n+1) 
def solve(v, ef, rhs, dhs):
    for i in range(50):
        if abs(v - rhs(ef)) < 1e-4:
            return ef
        else:
            ef +=  (v - rhs(ef))/dhs(ef) 
    return 0.0


###################################
#         GAUSSIAN WAVE           #
###################################
def gaussERF(ef, et, st, D):
    return D*(erf( (ef - et)/(m.sqrt(2)*st)) + erf( et / (m.sqrt(2)*st)) )

def gaussDIS(ef,et, st,D):
    return D*(1./m.sqrt(2)/st)*np.exp( -((ef - et)/(m.sqrt(2)*st))**2)

def fermiNormal(Vmax, traps, Vd = 0, npoints = 100):
    # Constants
    eta  = ((gsgv)/(4*m.pi))/((hbar*vf)**2)
    Cox  = (er*ep)/t
    C    = (e**2)/Cox #eV^{-1}

    vg  = np.asarray([Vmax*m.cos(i)-Vd for i in np.linspace(0,2*m.pi,npoints)])
    _vg   = np.asarray([Vmax*m.cos(i) for i in np.linspace(0,2*m.pi,npoints)])

    # Gaussian things 
    et, st, D = traps[0], traps[1], traps[2] 
    rhsp = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*gaussERF(ef,et,st,D)
    dhsp = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*gaussDIS(ef,et,st,D)
    rhsm = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*gaussERF(ef,-et,st, D)
    dhsm = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*gaussDIS(ef,-et,st,D)

    Ewave = [solve(vg[0], 0.0 , rhsp, dhsp)]
    for v in vg[1:]: 
        if v>0:
            Ewave.append( solve(v, Ewave[-1], rhsp, dhsp) )
        else: 
            Ewave.append( -1.0*solve(-v, abs(Ewave[-1]),  rhsm, dhsm) )
    check = []
    for ef in Ewave: 
        if ef>0:
            check.append(rhsp(ef)) 
        else:
            check.append(-1.0*rhsm(-1*ef)) 
    return _vg[1:], Ewave[1:], check[1:]



###################################
#      DOUBLE GAUSSIAN  WAVE      #
###################################
def gaussERF(ef, et, st, D):
    return D*(erf( (ef - et)/(m.sqrt(2)*st)) + erf( et / (m.sqrt(2)*st)) )

def gaussDIS(ef,et, st,D):
    return D*(1./m.sqrt(2)/st)*np.exp( -((ef - et)/(m.sqrt(2)*st))**2)

def fermiBivariate(Vmax, traps1, traps2, Vd = 0, npoints = 100):
    # Constants
    eta  = ((gsgv)/(4*m.pi))/((hbar*vf)**2)
    Cox  = (er*ep)/t
    C    = (e**2)/Cox #eV^{-1}

    vg  = np.asarray([Vmax*m.cos(i)-Vd for i in np.linspace(0,2*m.pi,npoints)])
    _vg   = np.asarray([Vmax*m.cos(i) for i in np.linspace(0,2*m.pi,npoints)])

    # Gaussian things 
    et1, st1, D1 = traps1[0], traps1[1], traps1[2] 
    et2, st2, D2 = traps2[0], traps2[1], traps2[2] 
    rhsp = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*gaussERF(ef,et1,st1,D1) + C*gaussERF(ef,et2,st2,D2)
    dhsp = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*gaussDIS(ef,et1,st1,D1) + C*gaussDIS(ef,et2,st2,D2)
    rhsm = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*gaussERF(ef,-et1,st1,D1) + C*gaussERF(ef,-et2,st2,D2)
    dhsm = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*gaussDIS(ef,-et1,st1,D1) + C*gaussDIS(ef,-et2,st2,D2)


    Ewave = [solve(vg[0], 0.2, rhsp, dhsp)]
    for v in vg[1:]: 
        if v>0:
            Ewave.append( solve(v, Ewave[-1], rhsp, dhsp) )
        else: 
            Ewave.append( -1.0*solve(-v, abs(Ewave[-1]),  rhsm, dhsm) )
    check = []
    for ef in Ewave: 
        if ef>0:
            check.append(rhsp(ef)) 
        else:
            check.append(-1.0*rhsm(-1*ef)) 
    return _vg[1:], Ewave[1:], check[1:]



################################
#        CAUCHY WAVE           #
################################
def cauchyERF(ef, et, g, D):
    return D* ( (1./m.pi)* np.arctan((ef - et)/g) - (1./m.pi)* np.arctan(-et/g))

def cauchyDIS(ef, et, g, D):
    return D* (1./ ( m.pi*g * ( 1  +  ((ef - et)/g)**2 )) )

def fermiCauchy(Vmax, traps, Vd = 0, npoints = 100):
    
    # Constants
    eta  = ((gsgv)/(4*m.pi))/((hbar*vf)**2)
    Cox  = (er*ep)/t
    C    = (e**2)/Cox #eV^{-1}

    vg  = np.asarray([Vmax*m.cos(i)-Vd for i in np.linspace(0,2*m.pi,npoints)])
    _vg   = np.asarray([Vmax*m.cos(i) for i in np.linspace(0,2*m.pi,npoints)])
    # Cauchy things 
    
    et, g, D = traps[0], traps[1], traps[2] 
    rhsp = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*cauchyERF(ef, et, g, D)
    dhsp = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*cauchyDIS(ef, et, g, D)
    rhsm = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*cauchyERF(ef, -et, g, D)
    dhsm = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*cauchyDIS(ef, -et, g, D)
   
    Ewave = [solve(vg[0], 0.1 , rhsp, dhsp)]
    for v in vg[1:]: 
        if v>0:
            Ewave.append( solve(v, Ewave[-1], rhsp, dhsp) )
        else: 
            Ewave.append( -1.0*solve(-v, abs(Ewave[-1]),  rhsm, dhsm) )
    check = []
    for ef in Ewave: 
        if ef>0:
            check.append(rhsp(ef)) 
        else:
            check.append(-1.0*rhsm(-1*ef)) 

    return _vg[1:], Ewave[1:], check[1:]


def fermiNormalCauchy(Vmax, Vd, traps1, traps2, points):
    # Constants
    eta  = ((gsgv)/(4*m.pi))/((hbar*vf)**2)
    Cox  = (er*ep)/t
    C    = (e**2)/Cox #eV^{-1}

    # Generate the voltage wave if needed
    if isinstance(points, int):
        vg   = np.asarray([Vmax*m.cos(i)-Vd for i in np.linspace(0,2*m.pi,npoints)])
        _vg  = np.asarray([Vmax*m.cos(i) for i in np.linspace(0,2*m.pi,npoints)])
    else:
        #vm_list = [i for i, j in enumerate(points) if abs(j - max(points)) < 0.01]
        #_vg = points[vm_list[0]:vm_list[1]]
        vg   = [v - Vd for v in points]
  
    # Gaussian things 
    et1, st1, D1 = traps1[0], traps1[1], traps1[2] 
    et2, g2, D2 = traps2[0], traps2[1], traps2[2] 
    rhsp = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*gaussERF(ef,et1,st1,D1) + C*cauchyERF(ef,et2,g2,D2)
    dhsp = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*gaussDIS(ef,et1,st1,D1) + C*cauchyDIS(ef,et2,g2,D2)
    rhsm = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*gaussERF(ef,-et1,st1,D1) + C*cauchyERF(ef,-et2,g2,D2)
    dhsm = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*gaussDIS(ef,-et1,st1,D1) + C*cauchyDIS(ef,-et2,g2,D2)

    Ewave = [solve(vg[0], 0.2, rhsp, dhsp)]
    for v in vg[1:]: 
        if v>0:
            Ewave.append( solve(v, Ewave[-1], rhsp, dhsp) )
        else: 
            Ewave.append( -1.0*solve(-v, abs(Ewave[-1]),  rhsm, dhsm) )
    check = []
    for ef in Ewave: 
        if ef>0:
            check.append(rhsp(ef)) 
        else:
            check.append(-1.0*rhsm(-1*ef)) 
    return points, Ewave, check


def fermiNormalCauchyNDF(vg, traps1, traps2, Vd = 0):
    # Constants
    eta  = ((gsgv)/(4*m.pi))/((hbar*vf)**2)
    Cox  = (er*ep)/t
    C    = (e**2)/Cox #eV^{-1}
    
    # Apply fermi shift
    vg   = vg - Vd
  
    # Gaussian things 
    et1, st1, D1 = traps1[0], traps1[1], traps1[2] 
    et2,  g2, D2 = traps2[0], traps2[1], traps2[2] 
    rhsp = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*gaussERF(ef,et1,st1,D1) + C*cauchyERF(ef,et2,g2,D2)
    dhsp = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*gaussDIS(ef,et1,st1,D1) + C*cauchyDIS(ef,et2,g2,D2)
    rhsm = lambda ef: ef + (eta*C/2.0)*(g*ef + ef**2) + C*gaussERF(ef,-et1,st1,D1) + C*cauchyERF(ef,-et2,g2,D2)
    dhsm = lambda ef: 1  + (eta*C/2.0)*(g + 2*ef) + C*gaussDIS(ef,-et1,st1,D1) + C*cauchyDIS(ef,-et2,g2,D2)

    if vg>0:
        ef = solve(vg, 0.1, rhsp, dhsp) 
    else: 
        ef = -1.0*solve(-vg, 0.1,  rhsm, dhsm)
    
    if ef>0:
        check = rhsp(ef) 
    else:
        check = -1.0*rhsm(-1*ef) 
    return ef, check

if __name__ == "__main__":
    traps = [-0.200, .01, 6.0e14]
    vg, ef, ck = fermiCauchy(7.0,traps, 1.75, 101 )

    traps = [-0.290, 0.140 , 4.0e13]
    vg2, ef2, ck2 = fermiNormal(7.0,traps, 1.75, 101)

    trapsA = [-0.270, 0.140 , 2.0e13]
    trapsB = [-0.300, 0.110 , 5.0e13]
    vg3, ef3, ck3 = fermiBivariate(7.0,trapsA, trapsB, 1.75, 101)
    
    trapsA = [-0.270, 0.140 , 2.0e13]
    trapsB = [-0.210, 0.001 , 1.0e15]
    vg4, ef4, ck4 = fermiNormalCauchy(7.0,trapsA, trapsB, 1.75, 101)

    fig = plt.figure(1) 
    ax1= fig.add_subplot(111)
    ax2= ax1.twinx()
    ax1.plot(vg, 'k') 
    ax1.plot(ck, 'kx') 
    ax1.plot(ck2, 'kx') 
    ax1.plot(ck3, 'kx') 
    ax1.plot(ck4, 'kx') 

    h1,=ax2.plot(ef,'r') 
    h2,=ax2.plot(ef2,'g')
    h3,=ax2.plot(ef3,'b')
    h3,=ax2.plot(ef4,'c')
    ax1.set_ylabel("Waveform $(V)$")
    ax2.set_ylabel("Fermi Energy $(eV)$")
    #plt.legend([h1,h2,h3],["Gaussain","Cauchy","Bimodal"])
    plt.show()

    #X = np.linspace(-3,3) 
    #Y = [erf(i) for i in X]
    #plt.plot(X,Y) 
    #plt.show()
