/* Title: Getting y data from Basilisk file
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
char filename[80];
event init(t = 0)
{
  restore (file = filename);
  boundary(all);
  FILE * fp = ferr;
  int LEVEL = 7;
  for (double y = -0.5; y < 0.5; y += 1./pow(2.,LEVEL)){
    fprintf(ferr, "%g %g\n", y, interpolate(u.x, 0.0, y));
  }
  fflush (fp);
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
  run();
}
