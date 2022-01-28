#include "tsd/tsd.hpp"
#include "tsd/wav.hpp"
#include <assert.h>

int test_wav()
{
  //tsd::audio::WavLecteur lecteur;
  //tsd_assert(lecteur.charge("./data/audio/man.wav") == 0);
  //lecteur.

  auto [x, fech] = tsd::audio::wav_charge("./data/man.wav");
  tsd_assert(x.rows() > 0);
  tsd_assert(fech > 0);

  {
    tsd::audio::WavEcrivain ecrivain;
    ecrivain.init("./build/tmp1.wav", {1, fech, tsd::audio::WavConfig::PCM_16});
    ecrivain.ecris(x);
  }

  {
    tsd::audio::WavEcrivain ecrivain;
    ecrivain.init("./build/tmp2.wav", {1, fech, tsd::audio::WavConfig::PCM_FLOAT});
    x /= 32768.0f;
    ecrivain.ecris(x);
  }

  {
    tsd::audio::WavEcrivain ecrivain;
    ecrivain.init("./build/tmp3.wav", {2, fech, tsd::audio::WavConfig::PCM_FLOAT});

    tsd::ArrayXcf z = x;
    z += tsd::cfloat(0,1) * x;

    ecrivain.ecris(z);
  }

  msg_avert("test_wav() : TODO : automatiser : le test !!!!");

  return 0;
}
