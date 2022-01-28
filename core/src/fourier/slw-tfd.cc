#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"


namespace tsd::fourier
{

  struct SLWTFDConfig
  {
    int n;
  };

struct SLWTFD
{
  ArrayXcf  twid;
  ArrayXcf  dft;
  ArrayXXcf table;
  int idx = 0;

  int configure(const SLWTFDConfig &config)
  {
    dft   = ArrayXcf::Zero(config.n);
    table = ArrayXXcf::Zero(config.n,config.n);
    twid  = polar(linspace(0,(2*π*(config.n-1))/config.n,config.n));
    return 0;
  }

  void step(IArrayXcf &x, ArrayXcf &y)
  {
    dft = twid * (dft + x - table.col(idx));
    y = dft;
    table.col(idx) = x;
    idx++;
  }

};




}
