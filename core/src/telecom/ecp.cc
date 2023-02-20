#include "tsd/telecom.hpp"
#include "tsd/vue.hpp"
#include <cmath>

using namespace tsd::vue;

namespace tsd::telecom {


struct ECP: Filtre<cfloat, cfloat, ECPConfig>
{
  sptr<SourceGen<cfloat>> ol;

  ECP(const ECPConfig &config)
  {
    configure(config);
  }
  entier configure_impl(const ECPConfig &config)
  {
    msg("ECP : fe = {} Hz, doppler = {} Hz (soit {:.3f} degrés / échantillon), Eb/N0 = {:.1f} dB, fsymb = {}, fbit = {}.",
        config.fe, config.décalage_fréquence, 360.0f * config.décalage_fréquence / config.fe,
        config.Eb_N0,
        config.fsymb, config.fbit);

    ol = source_ohc(config.décalage_fréquence / config.fe);

    retourne 0;
  }
  void step(const Veccf &x, Veccf &y)
  {
    //soit &config = Configurable<ECPConfig>::config;

    msg("ECP ({} échantillons)...", x.rows());

    // (1) Niveau du signal
    soit niveau = sqrt(abs2(x).moyenne());

    msg("  Décalage phase & fréquence...");
    // (2) Décalage phase & fréquence
    soit n = x.rows();

    y = x * ol->step(n);
    y *= std::polar(1.0f, config.décalage_phase);


    // (3) Bruit blanc

    // (3a) Compute standard deviation sigma from Eb/N0
    soit ebn0_lin = pow(10.0f, config.Eb_N0 / 10.0f);

    // eb/n0 = 1 / (2*σ)
    // => σ = level * sqrt(0.5 / eb/n0)

    // OVS effect:
    // signal amplitude *= OVS
    // noise  amplitude *= sqrt(OVS)
    // ===> Eb/N0 in amplitude *= sqrt(OVS)
    // ===> Eb/N0 in energy *= OVS
    // ====> Must compensate

    // Attention : ceci fonctionne pour Es/N0
    soit σ = niveau * sqrt(0.5f * (config.fe / config.fbit) / ebn0_lin);

    msg("ECP : niveau détecté = {}, Eb/N0 = {:.1f} dB, σ = {}", niveau, config.Eb_N0, σ);

    msg("  bruit...");
    // (3b) Add the noise
    y = tsd::telecom::bruit_awgn(y, σ);

    // (4) Apply some delay

    soit delais = (entier) config.délais_horloge;

    si(delais < 0)
    {
      Veci r = randi(config.fe / config.fsymb, 1);
      delais = r(0);
    }

    si(delais > 0)
      y = y.tail(n - delais).eval();


    // TODO : head / tail
    //pour(i = 0; i + config.clock_delay < n; i++)
      //y(i) = y(i + (entier) config.clock_delay);


    /*dsp_check(n > config.clock_delay,
              "n = %d, clock_delay = %f.",
              n, (float) config.clock_delay);*/

    //tsd_assert(n > config.clock_delay);

    //y.resize(n - config.clock_delay);

    si(config.debug_actif)
    {
      Figures f;
      f.subplot().plot(x, "", "Entrée");
      f.subplot().plot_psd(x, config.fe);
      f.subplot().plot(y, "",
          format("Sortie Eb/N0={:.1g} dB, niveau = {:.2g}, σ = {:.3g}", config.Eb_N0, niveau, σ));
      f.subplot().plot_psd(y, config.fe);
      f.afficher("ECP");

    }
  }
};

sptr<Filtre<cfloat, cfloat, ECPConfig>> ecp_création(const ECPConfig &config)
{
  retourne std::make_shared<ECP>(config);
}


}
