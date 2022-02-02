///////////////////////////////////////////////////////////////////////////
// Simple vérification que l'API en langue anglaise compile bien.        //
// (les vrais tests sont faits sur l'API racines - en langue française)  //
///////////////////////////////////////////////////////////////////////////


#include "dsp/dsp.hpp"
#include "dsp/fourier.hpp"
#include "dsp/filter.hpp"
#include "dsp/time.hpp"
#include "dsp/figure.hpp"
#include "dsp/telecom.hpp"

using namespace dsp;
using namespace dsp::fourier;
using namespace dsp::filter;
using namespace dsp::time;
using namespace dsp::view;
using namespace dsp::telecom;



/*struct Fr
{
  int champs_fr;
};

struct En: Fr
{
  int &en_field = champs_fr;
};*/





int test_dsp()
{
  /*En en;
  en.en_field = 5;
  msg("en field = {}, fr = {}", en.en_field, en.champs_fr);*/


  // dsp::filter
  {
    // TODO
  }

  // dsp::fourier
  {
    FFTFilterConfig c;
    c.freq_domain_processing = [](ArrayXcf &)
    {

    };
    auto f = filter_fft(c);
  }

  // dsp::time
  {
    auto t = DateTime::now();
    msg("Now = {}", t);
    t += Duration::microseconds(100);
  }

  // dsp::view
  {
    Figure f;
    ArrayXf t = linspace(0, 19, 20);
    f.plot(t, t.sqrt(), "b-o");
    f.show();
  }

  // dsp::telecom
  {
    auto spec = ShapingFilterSpec::srrc(0.3);
    auto filt = spec.matched_filter(63, 4);
    Figure f;
    f.plot(spec.get_coefs(63, 4), "", "Matched filter");
    f.show();
  }

  return 0;
}

