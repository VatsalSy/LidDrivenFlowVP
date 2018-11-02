# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
matplotlib.rcParams['text.latex.unicode'] = True


def gettingfield(filename, LEVEL):
    print('Getting field values')
    exe = ["./getData", filename, str(LEVEL), str(mu0), str(tauy), str(mumax)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    Xtemp = []
    Ytemp = []
    psitemp = []
    ftemp = []
    for n1 in range(1,len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Xtemp.append(float(temp3[0]))
            Ytemp.append(float(temp3[1]))
            psitemp.append(float(temp3[2]))
            ftemp.append(float(temp3[3]))
    X = np.asarray(Xtemp)
    Y = np.asarray(Ytemp)
    psi = np.asarray(psitemp)
    f = np.asarray(ftemp)

    N = int(np.sqrt(len(X)))
    X.resize((N, N))
    Y.resize((N, N))
    psi.resize((N, N))
    f.resize((N, N))

    print('Got field values')
    return X, Y, psi, f

# ----------------------------------------------------------------------------------------------------------------------


place = "tau0"
name = "tau0.png"
LEVEL = 6
mu0 = 1.0
tauy = 0.0
mumax = 10**4

X, Y, psi, f = gettingfield(place, LEVEL)
f = np.log(f)/np.log(10)
## Part to plot
plt.close()
fig, ax = plt.subplots()
fig.set_size_inches(19.20, 10.80)
rc('axes', linewidth=2)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
# cMap = matplotlib.colors.ListedColormap(['b','r'])
cntrl = ax.pcolormesh(X, Y, f, cmap="RdBu_r",edgecolor='face', vmax = f.max(), vmin = f.min())
cb1 = fig.add_axes([0.8, 0.1, 0.03, 0.8])
c1 = plt.colorbar(cntrl,cax=cb1)
c1.set_label(r'$log\left(\|D_{ij}\|\right)$', fontsize=30,labelpad=30)
c1.ax.tick_params(labelsize=20)
c1.ax.tick_params(labelright=True)
ax.contour(X, Y, psi, 20, colors='black', linewidths=2)
ax.set_xlabel(r'$X/D$', fontsize=30)
ax.set_ylabel(r'$Y/D$', fontsize=30)
ax.set_aspect('equal')
ax.set_ylim(X.min(), X.max())
ax.set_xlim(Y.min(), Y.max())

# plt.show()
plt.savefig(name,dpi=300)
plt.close()
