# LidDrivenFlowVP
This repository contains implementation of Viscoplastic viscosity in Basilsik C.
Vola et al. (2003) used Augumented Lagrangian Method to simulate lid driven cavity flow. The method has been applied in Basilisk C (http://basilisk.fr/sandbox/popinet/lid-bingham.c)
Here, I implement the regularization method of implementing Viscoplastic fluid.

![VP equation](https://latex.codecogs.com/gif.latex?%5Cmu_%7Bequi%7D%20%3D%20%5Cleft%28%5Cfrac%7B%5Ctau_y%7D%7B2%5C%7CD_%7Bij%7D%5C%7C%7D%5Cright%29%20&plus;%20%5Cmu_0)
