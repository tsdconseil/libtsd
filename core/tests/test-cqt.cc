#include "tsd/tsd-all.hpp"
#include "tsd/fourier/cqt.hpp"
#include "tsd/tests.hpp"

using namespace tsd::tf;
using namespace tsd::tf::cqt;

void test_cqt()
{
  msg_majeur("Test CQT...");

  soit fe = 48.0e3f;
  soit chirp = sigchirp(1e3/fe, 20e3/fe, 10 * fe, 'l');

  //soit [x, fe] = tsd::audio::wav_charge("c:/dbi/samples/piano scale.wav");
  soit A = tsd::tf::periodogramme_cqt(chirp, fe, 27.5, 20e3, pow(2.0f, 1.0f/12), 20);

  si(tests_debug_actif)
  {
    Figure fig("CQT");
    fig.plot_img(A, "mono");
    fig.afficher();
  }
}
