#        Python Settings
from __future__ import division
from scipy.optimize import newton
import numpy as np
#np.set_printoptions(threshold=np.nan)
import scipy.special as sp
import re
import matplotlib.pyplot as plt
import math
import os

###############################################################################
#                             Problem 2
###############################################################################
SigmaT = 1.0
SigmaS0 = 0.5
SigmaS1 = [-0.1, 0.0, 0.1]
SigmaA = 0.5
phileft = 1
length = 0
x = np.linspace(0,10,1000)
delta = x[1]-x[0]
def MCDiffusion(SigmaS1):
    D = 1/(3*(SigmaT - SigmaS1))
    L = np.zeros((len(x), len(x)))
    for i in range(0,len(x)):
        for j in range(0,len(x)):
            if i == j:
                L[i,j] = 2*(D/delta) + SigmaA*delta
            elif j == i+1:
                L[i,j] = -D/delta
            elif j == i-1:
                L[i,j] = -D/delta
    L[0,0] = (2*D)/((4*D+delta)*delta)
    L[len(x)-1,len(x)-1] = (2*D)/((4*D+delta)*delta)
    S = np.zeros((len(x), 1))
    S[0,0] = 1/delta
    flux = np.linalg.solve(L,S)
    print flux[len(x)-1]
    return flux

flux0 = MCDiffusion(SigmaS1[0])
flux1 = MCDiffusion(SigmaS1[1])
flux2 = MCDiffusion(SigmaS1[2])

#plt.figure
#params = {'mathtext.default': 'regular' }          
#plt.figure(figsize = (15,10.5))
#plt.rcParams.update(params)
#
#plt.xlabel("x (cm)",  fontname="Arial", fontsize=30)
#plt.tick_params(which='major', length=15, labelsize=25)
#plt.tick_params(which='minor', length=7)
##grid(b=True, which='major', color='light grey', linestyle='-')
#plt.grid(True, which='minor', color='lightgrey', linestyle='-')
#plt.grid(True, which='major', color='dimgrey', linestyle='-')
#
#plt.ylabel('Flux', fontname="Arial", fontsize=30)
#plt.title ("Mesh Centered Diffusion",fontsize=30)
#plt.rc('font',family='Arial')
#
#p1 = plt.plot(x, flux0, 'k-', label = 'flux -0.1', linewidth = 4)
#p2 = plt.plot(x, flux1, 'g-', label = 'flux 0.0', linewidth = 4)
#p3 = plt.plot(x, flux2, 'r-', label = 'flux 0.1', linewidth = 4)
#
#plt.legend(loc=1,prop={'size':20})

###############################################################################
###############################################################################


###############################################################################
#                           Problem 3
###############################################################################

def SolveLegendre(l, x):
    ans = sp.legendre(l)(x)
    return ans
def SolveLegendreDerivative(l, x):
    legendre = sp.legendre(l)
    ans = np.polyder(legendre)(x)
    return ans
def Weight(l, x):
    w=2*(1-x**2)/((l + 1)**2*(SolveLegendre(l+1, x))**2)
    return w
def FindZeros(l):
    N = l
    legendre = sp.legendre(l)
    Zeros =[]
    for i in range(1,N+1):
        x = np.cos((np.pi*(i-0.25))/(N+0.5))*(1-(1/8)*(N**(-2) - N**(-3)))
        ans = newton(legendre, x)
        Zeros.append(ans)
    return Zeros
def Verify(l,x, count):
    check1 = (l + 1)*SolveLegendre(l+1,x) - (2*l+1)*x*SolveLegendre(l,x) + l*SolveLegendre(l-1,x)
    check2 = (1-x**2)*SolveLegendreDerivative(l,x) + l * x*SolveLegendre(l,x) - l * SolveLegendre(l-1,x)
    if round(check1,10) == 0 and round(check2,10) == 0:
        count = count +1
    return count


l = [32]
g = SolveLegendre(18,0)
#print g
q = SolveLegendreDerivative(3,0)
#print q
for k in l:
    t = FindZeros(k)
z =[]
summ=0
for i in range(0,l[0]):
    z.append(Weight(l[0], t[i]))
    # The summation acts as a check on the weighting
    summ = z[i] + summ


###############################################################################
#                             Problem 4 from HW4
###############################################################################

eig=1

for blah in range(10):
    SigmaT = 1.0
    SigmaS0 = 0.5
    SigmaS1 = [-0.1, 0.0, 0.1]
    SigmaA = 0.5
    Sigmaf = 0.4
    Sigmac = 0.1
    nu=2.5
    x = np.linspace(0,10, 1000)
    delta=x[1]-x[0]
    k=range(0,32)
    #while abs(phi - 4.54e-5)/(4.54e-5) > 0.001:
    flux=[]
    totNFM = 1000
    psi = np.zeros((totNFM+1,1))
    print psi
    phi = np.zeros((totNFM, 1))
    Q = np.zeros((1, totNFM))
    A = np.zeros((1,16))
    B = np.zeros((1,16))
    alpha = []
    for i in range(0,16):
        alpha.append(abs(t[i])/t[i])
    for k in range(0,16):
        A[0,k] = (2*t[k] - SigmaT * delta)/(2*t[k]+SigmaT*delta)*alpha[k]
        B[0,k] = (2*delta)/(2*t[k]+SigmaT*delta)*alpha[k]
    eps_phi = 0.001
    maxit = 200
    err = 1
    it = 0
    while (err>eps_phi and it < maxit):
        phi0 = phi
        for i in range(0,totNFM):
            Q[0,i] = (SigmaS1[1]*phi[i] + i/eig * nu * Sigmaf *phi[i]+ 1)
        for j in range(0, totNFM):
            calc=0
            for k in range(0,16):
                calc = A[0,k] *psi[j,0] + B[0,k]*Q[0,j]
            psi[j+1,0] = calc
        for i in range(1,totNFM):
            calc=0
            for k in range(0,16):
                calc = z[k] * 0.5 *psi[i,0] + calc
            phi[i,0] = calc
            
        err_phi = abs((phi-phi0)/phi0)
        it = it+1
    if (it <= maxit):
        print 'Converged in ' + str(it) + ' iterations'
    else:
        print 'Failed to converge'
        
    x = x[1:]
    eig=phi[0]
    phi = phi[1:]
    
#    for i in range(0,15):
#        print phi[i]
        
        
    
        
    
plt.figure
params = {'mathtext.default': 'regular' }          
plt.figure(figsize = (15,10.5))
plt.rcParams.update(params)

plt.xlabel("x (cm)",  fontname="Arial", fontsize=30)
plt.tick_params(which='major', length=15, labelsize=25)
plt.tick_params(which='minor', length=7)
#grid(b=True, which='major', color='light grey', linestyle='-')
plt.grid(True, which='minor', color='lightgrey', linestyle='-')
plt.grid(True, which='major', color='dimgrey', linestyle='-')

plt.ylabel('Flux', fontname="Arial", fontsize=30)
plt.title ("Mesh Centered Diffusion",fontsize=30)
plt.rc('font',family='Arial')
#plt.errorbar(phits_eve,phits_tave,yerr=phits_FError, linestyle="None",capsize=8)
#plt.errorbar(HZETRN_eave, HZETRN_Results,xerr=HZETRN_err, linestyle="None", ecolor='r', capsize=0)

p1 = plt.plot(x, phi, 'k-', label = 'flux', linewidth = 4)

plt.legend(loc=1,prop={'size':20})
#    

print eig
    
###############################################################################
###############################################################################