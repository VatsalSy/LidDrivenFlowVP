/* Title: Getting data from Basilisk file
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"

char filename[80];
int LEVEL;
double mu_0, tauy, mumax;
scalar unyielded[], psi[], omega[], D2p[];

event init(t = 0)
{
  restore (file = filename);
  boundary(all);
  N = pow(2,LEVEL);
  // solve for the streamfunction
  psi[top] = dirichlet(0);
  psi[bottom] = dirichlet(0);
  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);
  foreach() {
    omega[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
    psi[] = 0.;
  }
  boundary ({omega,psi});
  poisson (psi, omega);
  // double muTemp = mu_0;
  foreach() {
    double D11 = (u.x[1,0] - u.x[-1,0]);
    double D22 = (u.y[0,1] - u.y[0,1]);
    double D12 = 0.5*(u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,1]);
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(2*Delta);
    D2p[] = D2/sqrt(2.0);
    // if (D2 > 0){
    //   double temp = tauy/(sqrt(2.0)*D2) + mu_0;
    //   muTemp = min(temp, mumax);
    // } else {
    //   if (tauy > 0.0){
    //     muTemp = mumax;
    //   } else {
    //     muTemp = mu_0;
    //   }
    // }
    // D2p[] = muTemp;
    // if (D2 == 0) {
    //   unyielded[] = 1.0;
    // } else {
    //   double temp = tauy/(sqrt(2.0)*D2) + mu_0;
    //   if (temp >= mumax){
    //     unyielded[] = 1.0;
    //   }
    // }
  }
  FILE * fp = ferr;
  output_field ({psi, D2p}, fp);
  fclose (fp);
}

int main(int a, char const *arguments[])
{
  // Top moving wall
  u.t[top] = dirichlet(1);
  /**
  For the other no-slip boundaries this gives */
  u.t[bottom] = dirichlet(0);
  u.t[left]   = dirichlet(0);
  u.t[right]  = dirichlet(0);
  uf.n[left]   = 0;
  uf.n[right]  = 0;
  uf.n[top]    = 0;
  uf.n[bottom] = 0;
  sprintf (filename, "%s", arguments[1]);
  LEVEL = atoi(arguments[2]);
  mu_0 = atof(arguments[3]);
  tauy = atof(arguments[4]);
  mumax = atof(arguments[5]);
  run();
}
