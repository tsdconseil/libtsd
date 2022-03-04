#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/telecom.hpp"
#include "tsd/vue.hpp"

using namespace std;
using namespace tsd::vue;

namespace tsd::telecom {





ArrayXf doppler_distri(ArrayXd f, float fd, double fc)
{
  int n = f.rows();
  ArrayXf S(n);
  msg("fd = {}, fc = {}", fd, fc);
  for(auto i = 0; i < n; i++)
  {
    if(abs(f(i)-fc) >= fd)
      S(i) = 0;
    else
      S(i) = 1 / (π*fd*sqrt(1.0-pow((f(i)-fc)/fd, 2.0)));
  }
  return S;
}


// ?
ArrayXf doppler_filtre(float fd, float fs)
{
  //auto fs2 = 4 * fd; // OSF = 4
  int ntaps = 512;
  ArrayXd f = linspace(0, 2*fd, ntaps/2).cast<double>();
  ArrayXf S = doppler_distri(f, fd, 0);
  ArrayXf h = tsd::filtrage::design_rif_freq(ntaps, S);
  h /= h.sum();
  return h;
}



struct CanalDispersif: Filtre<cfloat, cfloat, CanalDispersifConfig>
{
  ArrayXf hd;

  sptr<FiltreGen<cfloat>> rif, reechan;

  ArrayXcf reste;
  bool premier_appel = true;

  int configure_impl(const CanalDispersifConfig &config)
  {
    return 0;
  }

  CanalDispersif(const CanalDispersifConfig &config)
  {
    Configurable<CanalDispersifConfig>::configure(config);

    hd = doppler_filtre(config.fd, config.fe);

    float fs2 = 4 * config.fd;

    tsd::filtrage::analyse_filtre(hd, config.fd).afficher();

    rif     = tsd::filtrage::filtre_rif<float,cfloat>(hd);
    reechan = tsd::filtrage::filtre_reechan<cfloat>(config.fe / fs2);


  }

  ArrayXcf gen_bruit(int n)
  {
    auto &config = Configurable<CanalDispersifConfig>::config;
    ArrayXcf b(n);
    b.real() = randn(n) / sqrt(2);
    b.imag() = randn(n) / sqrt(2);

    if(config.type == TypeCanal::RICE)
    {
      // x = a + n
      // K = a²/σ²
      // => a = σ sqrt(K)
      b += sqrt(config.K);
      // normalisation :
      //noise = noise ./ sqrt(chn.K ^ 2 + 1);
      b /= sqrt(b.square().mean());
    }
    return b;
  }

  void step(const Eigen::Ref<const Vecteur<cfloat>> x, Vecteur<cfloat> &y)
  {
    auto &config = Configurable<CanalDispersifConfig>::config;
    int n   = x.rows();
    int nr  = reste.rows();
    // (1) Generate m samples
    float fs2 = 4 * config.fd;
    int m = ceil((n-nr+1) * fs2 / config.fe);

    //int nh = hd.rows();

    if(premier_appel)
    {
      premier_appel = false;
      ArrayXcf x1 = rif->step(gen_bruit(hd.rows()));
      reechan->step(x1);
    }

    ArrayXcf x0 = gen_bruit(m);
    ArrayXcf x1 = rif->step(x0);
    ArrayXcf x2 = reechan->step(x1);



    //tsd_assert(x2.rows() >= x.rows());

    /////////////////////////

    y.resize(n);

    if(nr > 0)
      y.head(nr) = x.head(nr) * reste;

    tsd_assert(x2.rows() >= n - nr);

    y.tail(n - nr) = x.tail(n - nr) * x2.head(n - nr);

    reste = x2.tail(x2.rows() - (n - nr));

    {
      Figures f;
      f.subplot().plot(x, "", "x (entrée)");
      f.subplot().plot(x0, "", "x0 (bruit)");
      f.subplot().plot(x1, "", "x1 (filtre Doppler)");
      f.subplot().plot(x2, "", "x2 (ré-échan)");
      f.subplot().plot(y, "", "y (produit)");
      f.afficher();
    }

    // ArrayXcf x1 = tsd::filtrage::filtrer(hd, noise);
    // x1 = x1.segment(3*nh/2, m).eval();
    // (4) upsample to fs
    //ArrayXcf x2 = resample(x1, fs / fs2);
  }
};


sptr<Filtre<cfloat, cfloat, CanalDispersifConfig>> canal_dispersif(const CanalDispersifConfig &config)
{
  return std::make_shared<CanalDispersif>(config);
}





float bruit_thermique(float bp, float T)
{
  T += 273.15; // Conversion en Kelvin
  float kb = 1.380650e-23; // en Joules
  return kb * T * bp; // en W
}


}
