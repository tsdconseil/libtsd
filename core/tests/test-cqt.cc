#include "tsd/tsd.hpp"
#include "tsd/figure.hpp"
#include "tsd/fourier.hpp"
#include "tsd/fourier/cqt.hpp"
#include "tsd/wav.hpp"
#include "tsd/tests.hpp"


using namespace tsd;
using namespace tsd::tf;
using namespace tsd::tf::cqt;
using namespace tsd::fourier;
using namespace tsd::vue;


ArrayXf cqt_chirp(float fs, float fmin, float fmax, float duration)
{
// Build a linear chirp test signal.
//
// Calling Sequence
// x = cqt_chirp(fs,fmin,fmax,duration)
//
// Parameters
// fs: Sampling frequency (in Hz)
// fmin: Start frequency (in Hz)
// fmax: End frequency (in Hz)
// duration: Output signal duration, in second
// x: Output chirp signal
//
// Description
// This function can be used for test purposes.
// It builds a simple linear chirp signal, as defined by :
//
// <latex>$$x = \sin\left(2\pi \int f(t)dt\right)$$</latex>,
// with <latex>$$f(t) = f_{min} + \frac{t}{\textrm{duration}}\cdot\left(f_{max}-f_{min}\right)$$</latex>
//
// <refsection><title>Example</title></refsection>
// <programlisting>
//// Sampling frequency = 8 KHz
//fs   = 8e3;
//// Frequency range = from 1 Hz to 10 Hz
//fmin = 1; fmax = 10;
//// Generate a linear chirp (1 second long)
//x = cqt_chirp(fs, fmin, fmax, duration = 1);
//clf(); plot(x);
// </programlisting>
// <mediaobject><imageobject><imagedata fileref="ex_cqt_chirp.png" format="PNG"/></imageobject><caption><para>Example of linear chirp</para></caption></mediaobject>
//
// Authors
//  J.A., full documentation available on <ulink url="http://www.tsdconseil.fr/log/scriptscilab/cqt">http://www.tsdconseil.fr/log/scriptscilab/cqt</ulink>

    int n = ceil(duration * fs);
    ArrayXf f = linspace(fmin, fmax, n) / fs; // Frequency (in Hz)
    ArrayXf phase = 2 * π * cumsum(f); // Phase (in radians)
    return phase.sin();
}





int test_cqt()
{
  msg_majeur("Test CQT...");

  ArrayXf chirp = cqt_chirp(48e3, 1e3, 20e3, 10);
  //cqt.step(chirp);
  //auto [x, fe] = tsd::audio::wav_charge("c:/dbi/samples/piano scale.wav");

  ArrayXXf A = tsd::tf::periodogramme_cqt(chirp, 48e3, 27.5, 20e3, std::pow(2.0f, 1.0f/12), 20);

  if(tests_debug_actif)
  {
    Figure fig("CQT");
    fig.plot_img(A, "mono");
    fig.afficher();
  }
  return 0;
}
