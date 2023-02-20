///////////////////////////////////////////////////////////////////////////
// Simple vérification que l'API en langue anglaise compile bien.        //
// (les vrais tests sont faits sur l'API racines - en langue française)  //
///////////////////////////////////////////////////////////////////////////


#include "dsp/dsp-all.hpp"
#include "dsp/time.hpp"

using namespace dsp::time;



entier test_dsp()
{
  /*En en;
  en.en_field = 5;
  msg("en field = {}, fr = {}", en.en_field, en.champs_fr);*/


  {
    msg("Vérif affectation colonne d'un tableau...");
    let M  = Tabf::zeros(5,10);
    let v1 = Vecf::ones(5);
    let v2 = 2 * v1;

    msg("Affectation colonne 0 : get col(0)...");
    //let col0 = M.col(0);
    Vecf col0(M.col(0));
    col0.fr.dump_infos();
    msg("Affectation colonne 0 : copy (?)...");
    col0 = v1;
    msg("Affectation colonne 1...");
    M.col(1) = v2;

    msg("M = \n{}", M);
    msg("v1 = {}", v1);
    msg("v2 = {}", v2);

    dsp_assert(M(0,0) == 1);
    dsp_assert(M(0,1) == 2);
  }



  // dsp::filter
  {
    // TODO
  }

  // dsp::fourier
  {
    FFTFilterConfig c;
    c.freq_domain_processing = [](tsd::Veccf &)
    {

    };
    let f = filter_fft(c);
  }

  // dsp::time
  {
    let t = DateTime::now();
    msg("Now = {}", t);
    t += Duration::microseconds(100);
  }

  // dsp::view
  {
    Figure f;
    let t = linspace(0, 19, 20);
    f.plot(t, sqrt(t), "b-o");
    f.show();
  }

  // dsp::telecom
  {
    let spec = ShapingFilterSpec::srrc(0.3);
    let filt = spec.matched_filter(63, 4);
    Figure f;
    f.plot(spec.get_coefs(63, 4), "", "Matched filter");
    f.show();
  }

  // dsp::stats
  {
    let x = randn(10);
    let y = dsp::stats::levinson_real(x);
  }

  retourne 0;
}

