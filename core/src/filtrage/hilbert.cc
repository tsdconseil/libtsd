#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"

namespace tsd::filtrage {

//  A. Reilly, G. Frazer, and B. Boashash, "Analytic signal generation Tips and traps", IEEE Trans. Signal Processing, vol. 42, pp.3241-3245, Nov. 1994.
//  h(n)=(2/(n*pi))*(sin(n*pi/2))^2 for n not equal 0 and h(0)=0.
ArrayXf design_rif_hilbert(int n, const std::string &fen)
{
  ArrayXf h = ArrayXf::Zero(n);

  for(auto i = 0; i < n; i++)
  {
    float x = ((float) i) - (n/2);
    if(x == 0)
      h(i) = 0;
    else
      h(i) = (2/(π*x)) * carré(sin(π * x / 2));
  }

  return h * fenetre(fen, n);
}





struct HilbertTransformeur: Filtre<float, cfloat, HilbertTransformeurConfig>
{
  sptr<FiltreGen<float>> retard, hilbert;

  HilbertTransformeur(int ntaps, const std::string &fenetre)
  {
    HilbertTransformeurConfig config{ntaps, fenetre};
    configure(config);
  }
  int configure_impl(const HilbertTransformeurConfig &config)
  {
    ArrayXf h = design_rif_hilbert(config.ntaps, config.fenetre);
    hilbert   = filtre_rif<float>(h);
    retard    = ligne_a_retard<float>(config.ntaps/2);
    return 0;
  }
  void step(IArrayXf x, ArrayXcf &y)
  {
    auto x2 = hilbert->step(x);
    // Ligne à retard sur y, pour être synchrone avec le signal filtré
    auto y2 = retard->step(x);
    y = y2 + cfloat(0,1) * x2;
  }
};


sptr<Filtre<float, cfloat, HilbertTransformeurConfig>> hilbert_transformeur(int ntaps, const std::string &fenetre)
{
  return std::make_shared<HilbertTransformeur>(ntaps, fenetre);
}

ArrayXcf hilbert(IArrayXf x, int ncoefs)
{
  return hilbert_transformeur(ncoefs)->step(x);
}

}





