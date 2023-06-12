#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/telecom.hpp"
#include "tsd/vue.hpp"

using namespace std;
using namespace tsd::vue;

namespace tsd::telecom {



// TODO: A passer dans tsd
Veccf randnc(entier n)
{
  Veccf y(n);
  y.set_real(randn(n) / sqrt(2));
  y.set_imag(randn(n) / sqrt(2));
  retourne y;
}

Vecf doppler_distri(const Vecd &f, float fd, double fc)
{
  retourne Vecf::int_expr(f.rows(), IMAP(
      (abs(f(i)-fc) >= fd) ? 0 : 1 / (π*fd*sqrt(1.0-pow((f(i)-fc)/fd, 2.0)));
      ));
}


// ?
Vecf doppler_filtre(float fd, float fs)
{
  //soit fs2 = 4 * fd; // OSF = 4
  soit ntaps = 512;
  soit f = linspace(0, 2*fd, ntaps/2).as<double>();
  soit S = doppler_distri(f, fd, 0);
  soit h = tsd::filtrage::design_rif_freq(ntaps, S);
  h /= h.somme();
  retourne h;
}



struct CanalDispersif: Filtre<cfloat, cfloat, CanalDispersifConfig>
{
  Vecf hd;

  sptr<FiltreGen<cfloat>> rif, reechan;

  Veccf reste;
  bouléen premier_appel = oui;

  void configure_impl(const CanalDispersifConfig &config)
  {
  }

  CanalDispersif(const CanalDispersifConfig &config)
  {
    Configurable<CanalDispersifConfig>::configure(config);

    hd = doppler_filtre(config.fd, config.fe);

    soit fs2 = 4 * config.fd;

    tsd::filtrage::plot_filtre(hd, non, config.fd).afficher();

    rif     = tsd::filtrage::filtre_rif<float,cfloat>(hd);
    reechan = tsd::filtrage::filtre_reechan<cfloat>(config.fe / fs2);
  }

  Veccf gen_bruit(entier n)
  {
    soit &config = Configurable<CanalDispersifConfig>::config;

    soit b = randnc(n);

    si(config.type == TypeCanal::RICE)
    {
      // x = a + n
      // K = a²/σ²
      // => a = σ sqrt(K)
      b += sqrt(config.K);
      // normalisation :
      //noise = noise ./ sqrt(chn.K ^ 2 + 1);
      b /= sqrt(square(b).moyenne());
    }
    retourne b;
  }

  void step(const Vecteur<cfloat> &x, Vecteur<cfloat> &y)
  {
    soit &config = Configurable<CanalDispersifConfig>::config;
    soit n   = x.rows(), nr  = reste.rows();

    // (1) Génère m échantillons
    soit fs2 = 4 * config.fd;
    soit m = (entier) ceil((n-nr+1) * fs2 / config.fe);

    si(premier_appel)
    {
      premier_appel = non;
      reechan->step(rif->step(gen_bruit(hd.rows())));
    }

    soit x0 = gen_bruit(m);
    soit x1 = rif->step(x0);
    soit x2 = reechan->step(x1);

    y.resize(n);

    si(nr > 0)
      y.head(nr) = x.head(nr) * reste;

    assertion(x2.rows() >= n - nr);

    y.tail(n - nr) = x.tail(n - nr) * x2.head(n - nr);

    reste = x2.tail(x2.rows() - (n - nr));

    si(config.debug_actif)
    {
      Figures f;
      f.subplot().plot(x,  "", "x (entrée)");
      f.subplot().plot(x0, "", "x0 (bruit)");
      f.subplot().plot(x1, "", "x1 (filtre Doppler)");
      f.subplot().plot(x2, "", "x2 (ré-échan)");
      f.subplot().plot(y,  "", "y (produit)");
      f.afficher("Canal dispersif");
    }
  }
};


sptr<Filtre<cfloat, cfloat, CanalDispersifConfig>> canal_dispersif(const CanalDispersifConfig &config)
{
  retourne std::make_shared<CanalDispersif>(config);
}



float bruit_thermique(float bp, float T)
{
  T += 273.15; // Conversion en Kelvin
  soit kb = 1.380650e-23f; // en Joules
  retourne kb * T * bp; // en W
}


}
