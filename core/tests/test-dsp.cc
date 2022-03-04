///////////////////////////////////////////////////////////////////////////
// Simple vérification que l'API en langue anglaise compile bien.        //
// (les vrais tests sont faits sur l'API racines - en langue française)  //
///////////////////////////////////////////////////////////////////////////


#include "dsp/dsp-all.hpp"
#include "dsp/time.hpp"

using namespace dsp::time;



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

  // dsp::stats
  {
    ArrayXf x = randn(10);
    ArrayXf y = dsp::stats::levinson_real(x);
  }

  return 0;
}

