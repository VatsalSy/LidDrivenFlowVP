/* Title: Plannar Coeutte flow
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"

/* Values of yeld stress, viscosity and coefficient.
   Newtonian: $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 1
   Bingham $\mu_0 = 1.0; $\tau_y > 0.0$ and n = 1
*/
char filename[80];
double tauy;
int counter;
#define mu_0 (1.0)
// The regularisation value of viscosity
#define mumax (1e4)

int imax = 1e5;
#define LEVEL 6
#define Maxdt (1e-4)
#define ERROR (1e-10)

scalar un[];
face vector muv[];

// initialization event
event init (t = 0) {
  // preparing viscosity to be used as Non-Newtonian fluid
  mu = muv;
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
  }
  foreach(){
    un[] = u.x[];
  }
  dump (file = "start");
}


// int main(int arg, char const *arguments[])
int main()
{
  // sprintf (filename, "%s", arguments[1]);
  // tauy = 1000/sqrt(2.0);
  init_grid (1<<LEVEL);
  L0 = 1.0;
  origin (-0.5, -0.5);
  DT = Maxdt;
  TOLERANCE = 1e-5;
  // CFL number
  CFL = 0.5;
  for (counter = 0; counter < 5; counter++){
    if (counter == 0){
      tauy = 0.0;
    } else {
      tauy = pow(10, counter-1)/sqrt(2);
    }
    fprintf(ferr, "tauy = %g\n", tauy);
    sprintf (filename, "tau%d", counter);
    run();
  }
  // run();
}

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

/**
We look for a stationary solution. */
event logfile (i++; i <= imax) {
  double du = change (u.x, un);
  fprintf(ferr, "i = %d: dt = %g err = %g\n", i, dt, du);
  if (i > 0 && du < ERROR){
    dump (file = filename);
    return 1; /* stop */
  }
  if (i>imax-10){
    dump (file = filename);
  }
}
//
/**
## Implementation of generalized Newtonian viscosity
*/

event properties(i++) {
  /*
   $$D_{11} = \frac{\partial u}{\partial x}$$
   $$D_{12} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{21} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{22} = \frac{\partial v}{\partial y}$$
   The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
   $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
   the equivalent viscosity is
   $$\mu_{eq}= \mu_0 + \frac{\tau_y}{\sqrt{2} D_2 }$$
   **Note:** $\|D\| = D_2/\sqrt{2}$
   Finally, mu is the min of of $\mu_{eq}$ and a large $\mu_{max}$, then the fluid flows always, it is not a solid, but a very viscous fluid.
   $$ \mu = min\left(\mu_{eq}, \mu_{max}\right) $$
  */
  double muTemp = mu_0;
  foreach_face() {
    double D11 = (u.x[] - u.x[-1,0]);
    double D22 = ((u.y[0,1]-u.y[0,-1])+(u.y[-1,1]-u.y[-1,-1]))/4.0;
    double D12 = 0.5*(((u.x[0,1]-u.x[0,-1])+(u.x[-1,1]-u.x[-1,-1]))/4.0 + (u.y[] - u.y[-1,0]));
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);
    if (D2 > 1e-6) {
      double temp = tauy/(sqrt(2.0)*D2) + mu_0;
      muTemp = min(temp, mumax);
    } else {
      if (tauy > 0.0){
        muTemp = mumax;
      } else {
        muTemp = mu_0;
      }
    }
    muv.x[] = fm.x[]*(muTemp);
  }
  boundary ((scalar *){muv});
}
