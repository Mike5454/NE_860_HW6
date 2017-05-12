from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time


# In this problem we are using monte-carlo to solve for flux in a homogenous slab
# the slab is 10cm, and in the slab is an isotropic souce of 1/sec. The slab has
# vacuum conditions on both boundaries. I will use the collision estimator
# and the track length estimator using 10 cells and I will also compute the 
# current out both sides. Use enough histories for less than 1% ucertainty and 
# reporting the time 

# *************** Constants ******************
SigmaT = 1.0 #cm^-1
SigmaS = 0.5 #cm^-1
SigmaA = 0.5 #cm^-1

L = 10 #cm
ncells=10
Delta = L/ncells #cm
n = 10000000

# ****** RNG and time ******
#np.random.seed(1234)
rand = np.random.rand
t = time.time()

# ****** Initial Matrix ******
init = np.zeros(10)
N_collision=init
N_tracklength=init
S_collision=init
S_tracklength=init

# ****** Interaction terms ******
leakR=0
leakL=0
leak=0
Absorbed=0
Scatter=0
path=0
dmat = np.zeros((1,10))
dmat2 = np.zeros((1,10))

# ****** Monte Carlo ******
# For all particles
for i in range(n):
    # Sample direction, position, and velocity
    mu = 2*rand()-1
    weight = 1
    x = L*rand()
    # While particle is alive
    while True :
        # Define boundary
        if mu > 0 :
            boundary = L-x
        else :
            boundary = x
        # Sample distance traveled
        distance = -np.log(rand())/SigmaT * abs(mu)
        # Move particle to new location
        if abs(boundary) > abs(distance) :
            # Find where particle collides
            x2 = x + distance * mu/abs(mu)
#            for i in range(0,int(L/Delta)):
                # Add path length traveled to matrix
            mu = abs(mu)
            # For variance, I assume variance will be highest at the edge since it has the lowest flux value
            s1 = dmat[0,int(str(x)[0])]
            s3 = dmat2[0,int(str(x)[0])]
            # These sum the pathlengths in each cell for all cases
            if str(x)[0] < str(x2)[0]:
                dmat[0,int(str(x)[0])] = dmat[0,int(str(x)[0])] + (int(str(x)[0])+1 - x)/mu
                dmat[0,int(str(x2)[0])] = dmat[0,int(str(x2)[0])] + (x2 - int(str(x2)[0]))/mu
            elif str(x2)[0] < str(x)[0]:
                dmat[0,int(str(x2)[0])] = dmat[0,int(str(x2)[0])] + (int(str(x2)[0])+1 - x2)/mu
                dmat[0,int(str(x)[0])] = dmat[0,int(str(x)[0])] + (x - int(str(x)[0]))/mu
            elif str(x2)[0] == str(x)[0]:
                dmat[0,int(str(x)[0])] = dmat[0,int(str(x)[0])] + (abs(x2 - x))/mu
            # in the case a neutron skips an entire cell or more
            if abs(int(str(x)[0]) - int(str(x2)[0])) > 1:
                if int(str(x)[0]) < int(str(x2)[0]):
                    for i in range(int(str(x)[0])+1, int(str(x2)[0])):
                        dmat[0,i] = dmat[0,i] + 1/mu
                else:
                    for i in range(int(str(x2)[0])+1, int(str(x)[0])):
                        dmat[0,i] = dmat[0,i] + 1/mu
            # Collision Estimator and uncertainty
            dmat2[0,int(str(x2)[0])] = dmat2[0,int(str(x2)[0])] + 1
            s2 = dmat[0,int(str(x)[0])]
            s4 = dmat2[0,int(str(x)[0])]
            # Does it scatter?
            if rand() < SigmaS/SigmaT :
                # Give another direction
                mu = 2.0*rand()-1.0
                #weight2 = mu*weight + (1-mu**2)**(0.5)*(1-weight**2)**(0.5)*mu
            else :
                # If absorbed, particle is dead
                Absorbed += 1
                break
        else :
            # If leak, particle is dead
            ###################################
            #!!!!!! ASSUMPTION STATED!!!!!!!#
            ##################################
            # For the pathlength contributions from the boundary leakage particles I assumed
            # most of the contribution would come from the end cells 
            # i.e. I did not include the contributions of neutrons leaked put the boundary from cells 2-9
            # This simplifies the part of the code by a fair amound and with little difference
            # to the results since those collisions are not as probable.
            if mu > 0:
                leakR = leakR + weight
                dmat[0,L-1] = dmat[0,L-1] + boundary/abs(mu)
            else:
                leakL = leakL + weight
                dmat[0,0] = dmat[0,0] + boundary/abs(mu)
            break
        
# normalize 
flux = []
flux2=[]
for i in range(0,ncells):
    flux.append(dmat[0,i]/n)
    flux2.append(dmat2[0,i]/n)
    
uncertainty = abs(s1-s2)/s2
uncertainty2 = abs(s3-s4)/s4
    
print('\n')
#print('Cell Number ' + str(1))
print(" Leakage right: %.6e" % (float(leakR)/float(n)))
print(" Leakage left: %.6e" % (float(leakL)/float(n)))
print(" Fraction absorbed: %.6e" % (float(Absorbed)/float(n)))
print(" Particle balance: %.6e" % (float(Absorbed+leakR+leakL)/float(n)))
print(" Relative Uncertainty Pathlength: " + str(uncertainty*100))
print(" Relative Uncertainty Collisional: " + str(uncertainty2*100))

t = time.time()-t    
print("     Elapsed time: %.6e" % t)
print '\n'
print 'Pathlength Flux Values'
print flux
print 'Collisional Flux Values'
print flux2

plt.figure
params = {'mathtext.default': 'regular' }          
plt.figure(figsize = (15,10.5))
plt.rcParams.update(params)

plt.xlabel("Cell",  fontname="Arial", fontsize=30)
#plt.xscale('log')
plt.tick_params(which='major', length=15, labelsize=25)
plt.tick_params(which='minor', length=7)
#grid(b=True, which='major', color='light grey', linestyle='-')
plt.grid(True, which='minor', color='lightgrey', linestyle='-')
plt.grid(True, which='major', color='dimgrey', linestyle='-')

plt.ylabel(r'$\phi(x_i)$', fontname="Arial", fontsize=30)
#plt.yscale('log')

#plt.title ("CDF",fontsize=30)
plt.rc('font',family='Arial')
#plt.errorbar(phits_eve,phits_tave,yerr=phits_FError, linestyle="None",capsize=8)
#plt.errorbar(HZETRN_eave, HZETRN_Results,xerr=HZETRN_err, linestyle="None", ecolor='r', capsize=0)

p1 = plt.plot(range(1,11), flux, 'k-', label = 'Track Length', linewidth = 7)
p2 = plt.plot(range(1,11), flux2, 'r-', label = 'Collision', linewidth = 7)

plt.legend(loc=8,prop={'size':25})
plt.show
