# let's define a number of harmonic potentials
#
#     u_i(x) = 1/2 * k_i (x - x_i)^2 + c
#
# each with a different equilibrium position x_i and force constant k_i

import numpy as np

def u_i(x, k, x0, c):
    """Returns the (reduced) energy of the harmonic potential in units kT."""
    return 0.5*k*(x-x0)**2 + c

def delta_f_harmonic(k_i, k_j, dims=1):
    """Returns the difference in free energy f_ij = f_j - f_i (in units kT) of two harmonic potentials
    with force constants k_i and k_j.  For  $n$-dimensional harmonic potentials,

    -\ln (Z_j/Z_i) = - (n/2) * \ln [ (k_i/k_j) ]"""
    
    return -1.0*(float(dims)/2.0) * np.log( k_i/k_j )
    

#x_i = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.])   # length (unitless)
xshift = 0.5
x_i = np.array([0., 1., 2., 3., 4.,  5.+xshift, 6.+xshift, 7.+xshift, 8.+xshift, 9.+xshift])   # length (unitless)

ka, kb = 10.0, 500.0 
#k_i = np.array([ka, ka, ka, ka, ka, kb, ka, ka, ka, ka, ka])  # kT/(length)^2
k_i = np.array([ka, ka, ka, ka, ka,      ka, ka, ka, ka, ka])  # kT/(length)^2


#c_i = np.array([0.0, 0.0, 0.0, 0.0, 0.0, delta_f_harmonic(kb, ka), 0.0, 0.0, 0.0, 0.0, 0.0])  # kT/(length)^2
c_i = np.array([0.0, 0.0, 0.0, 0.0, 0.0,              0.0, 0.0, 0.0, 0.0, 0.0])  # kT/(length)^2
print(c_i)
n_ensembles = len(x_i)


# Let's make a plot of these

from matplotlib import pyplot as plt

plt.figure()
xvalues = np.arange(-2., 11., 0.05)
for i in range(n_ensembles):
    plt.plot(xvalues, u_i(xvalues, k_i[i], x_i[i], c_i[i]), label='i=%d'%i)
plt.xlim(0,10)
plt.ylim(-1,1)
plt.yticks(np.arange(-1, 1.5, 0.5))
plt.xlabel("x",fontsize=14)
plt.ylabel("$u_i(x)$ (in $k_BT$)",fontsize=14)
#plt.legend(loc='best')
proj_path = '/home/tun50867/work/git/one_over_t_project'
plt.savefig(f'{proj_path}/plots/harmonic_potential.pdf',format='pdf', dpi=600, bbox_inches='tight')

