#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"

namespace tsd {

// Plotting of a covariance ellipse around a 2d state
//
// Calling Sequence
// plot_cov_ellipse(m1,m2,s1,s2,s12)
//
// Parameters
//  m1,m2: Means (scalars)
//  s1,s2: Variances of each variable (scalars)
//  s12: Cross-variance (scalar)
//
// See also
//  plot_obs
//  plot_states
// Code d'après "Signal processing with SCILAB - INRIA"
void plot_cov_ellipse(tsd::vue::Figure &f, float m1, float m2, float s1, float s2, float s12)
{
  soit n = 500;
  soit t = linspace(0, 2*π, n);//0:.1:.1+%pi*2;
  soit c = 2*cos(t);
  soit s = 2*sin(t);
  soit ρ = s12 / sqrt(s1*s2);
  soit cr = sqrt(s1)*c + m1 * Vecf::ones(n);
  soit sr = sqrt(s2)*(ρ*c+sqrt(1-ρ*ρ)*s)+m2*Vecf::ones(n);
  f.plot(cr,sr, "-b");
}
}
