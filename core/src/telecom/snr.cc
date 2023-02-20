#include "tsd/telecom.hpp"

namespace tsd::telecom {

template <typename T> entier sgn(T val) {
    retourne (T(0) < val) - (val < T(0));
}

struct Matzner: EstimateurSNR
{
  float e2, e4;
  float γ;

  Matzner(float γ)
  {
    reset(γ);
  }
  void reset(float γ)
  {
    this->γ = γ;
    e2 = e4 = 0;
  }
  void step(const Veccf &x, Vecf &S, Vecf &N)
  {
    soit n = x.rows();
    S.resize(n);
    N.resize(n);
    pour(auto i = 0; i < n; i++)
    {
      soit e2i = std::norm(x(i));
      soit e4i = e2i * e2i;

      e2 = γ * e2i + (1 - γ) * e2;
      e4 = γ * e4i + (1 - γ) * e4;

      // n-PSK : ke = 1 (constant energy), kn = 2
      S(i) = sqrt(2 * e2*e2 - e4);
      N(i) = e2 - S(i);
      // EbN0 = S / (3 * N);     /// ????
    }
  }
};

/** @brief SNR estimation using absolute value of signal / pour BPSK signal */
struct BPSKAbsSNREstimator
{
  float eabs, en;
  float γ;

  BPSKAbsSNREstimator(float γ = 0.5f)
  {
    reset(γ);
  }

  void reset(float γ)
  {
    this->γ = γ;
    en = eabs = 0;
    msg("BPSK abs SNR estimation: γ={}.", γ);
  }

  /** @param z Signal d'entrée
   *  @param[out] a   Absolute value of the signal gain (energy of signal is the square of this)
   *  @param[out] en  Energy of the noise */
  void step(const Veccf &x, Vecf &S, Vecf &N)
  {
    entier n = x.rows();
    S.resize(n);
    N.resize(n);
    pour(auto i = 0; i < n; i++)
    {
      soit eabsi = abs(real(x(i)));
      soit eni   = norm(x(i) - sgn(real(x(i))) * eabs);
      eabs = γ * eabsi + (1-γ) * eabs;
      en = γ * eni   + (1-γ) * en;
      S(i) = eabs;
      N(i) = en;
    }
    //a  = abs(z);
    //en = std::norm(z - sgn(z) * a);
  }



};


sptr<EstimateurSNR> snr_Matzner(float γ)
{
  retourne std::make_shared<Matzner>(γ);
}


#if 0
/** @brief Detector based on the assumption that the squared signal should have
 *  a constant instantaneous frequency. */
template<typename T = float>
class BPSKCFreqSLoop
{
public:

  BPSKCFreqSLoop(qu0x16_t γ = 0.5f)
  {
    reset(γ);
  }

  void reset(qu0x16_t γ)
  {
    this->γ = γ;
    ugamma = qu0x16_t::one - γ;
    si(γ != 0)
      ugamma.v++;
    en = eabs = 0;

    lz = cfloat(1,0);
    dsp_trace("BPSK abs SNR estimation: γ = %f, uγ = %f.", (float) this->γ, (float) ugamma);
  }

  /** @param[out] a absolute value of the signal gain (energy of signal is the square of this)
   *  @param[out] en Energy of the noise */
  void step(T &a, T &en, const std::complex<T> &z)
  {
    cfloat sq, nmsq = z * z * conj(lz) * conj(lz);

    nmsq = γ * sq + ugamma * msq;

    en =

    msq = nmsq;

    float freq = arg(msq);

    a  = std::pow(std::norm(sq), 0.25);


    lz = z;
  }

private:

  std::complex<T> lz, msq;

  T eabs, en;
  qu0x16_t γ, ugamma;
};
#endif


}

