///////////////////////////////////////////////////////////////////////////
// Simple vérification que l'API en langue anglaise compile bien.        //
// (les vrais tests sont faits sur l'API racines - en langue française)  //
///////////////////////////////////////////////////////////////////////////


#include "dsp/dsp.hpp"
#include "dsp/fourier.hpp"
#include "dsp/filter.hpp"
#include "dsp/time.hpp"
#include "dsp/figure.hpp"

using namespace dsp;
using namespace dsp::fourier;
using namespace dsp::filter;
using namespace dsp::time;
using namespace dsp::view;

struct Fr
{
  int champs_fr;
};

struct En: Fr
{
  int &en_field = champs_fr;
};





int test_dsp()
{
  En en;
  en.en_field = 5;

  msg("en field = {}, fr = {}", en.en_field, en.champs_fr);

  {
    FFTFilterConfig c;
    c.freq_domain_processing = [](ArrayXcf &)
    {

    };
    auto f = filter_fft(c);
  }

  {
    auto t = DateTime::now();
    msg("Now = {}", t);
    t += Duration::microsecondes(100);
  }

  {
    Figure f;
    ArrayXf t = linspace(0, 19, 20);
    f.plot(t, t.sqrt(), "b-o");
    f.show();
  }

  return 0;
}

