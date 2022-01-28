// Démodulation générique ((a priori) modulation numérique)

#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"

using namespace std;
using namespace tsd::vue;

#define VERBOSE(AA)


namespace tsd::telecom {

using namespace tsd::filtrage;


struct ModGen : Modulateur
{
  // Fréquence d'échantillonnage interne
  float fe1;

  // Facteur de sur-échantillonnage
  unsigned int osf;

  // Oscillateur local si Fi != 0
  sptr<SourceGen<cfloat>> ol;

  // (Eventuel) adaptateur de rythme
  sptr<FiltreGen<cfloat>> ra;

  float latence = 0;

  bool valide = false;

  ModConfig config;

  sptr<FiltreGen<cfloat>> filtre_mise_en_forme;

  sptr<FormeOnde> wf;

  void def_fo(sptr<FormeOnde> wf)
  {
    this->wf = wf;
  }

  ModGen(const ModConfig &config)
  {
    valide = (configure(config) == 0);
  }

  float delais()
  {
    return latence;
  }

  int configure(const ModConfig &config)
  {
    latence = 0;
    this->config = config;

    wf = config.wf;

    msg("<h3>Configuration modulateur</h3>");
    msg("Configuration : fe={} Hz, fi={} Hz, fsymb={} Hz.",
          config.fe, config.fi, config.fsymb);

    fe1 = config.fe;

    bool auto_adapation_rythme = false;

    if(auto_adapation_rythme)
    {
      if((fe1 > (16 * config.fsymb)) || (fmod(config.fe,config.fsymb) != 0))
      {
        fe1 = 16 * config.fsymb;
        msg("Fréquence d'échantillonnage interne : {} Hz, adaptation de rythme : {}",
            fe1, config.fe / fe1);
        ra = filtre_reechan<cfloat>(config.fe / fe1);
        msg_avert("TODO : calcul latence après ra");
      }
    }


    osf = fe1 / config.fsymb;

    msg("Facteur de sur-échantillonnage : {}", osf);

    //stdo.printf("Configuration modulateur : fe=%g Hz, fi=%g Hz, fsymb=%g Hz (osf = %d)",
    //      config.fe, config.fi, config.fsymb, osf);



    //if(fmod(config.fe,config.fsymb) != 0)
    //{
      //erreur("mod_init: sample frequency must be a multiple of symbol frequency.");
      //return -1;
    //}

    int ncoefs = config.ncoefs_filtre_mise_en_forme;

    filtre_mise_en_forme = config.wf->filtre.filtre_mise_en_forme(ncoefs, osf);
    ol = source_ohc(config.fi / config.fe);

    {
      ArrayXf h = config.wf->filtre.get_coefs(ncoefs, osf);

      //tsd_assert(h.rows() == ncoefs);

      // Temps vers le milieu du premier bit transmis
      // Filtre = filtre_rif_ups, R = osf
      //latence = (h.rows() - 1) / 2.0f + osf - 1;
      latence = filtre_rif_ups_délais(h.rows(), osf);
      msg("Modulateur : calcule délais = {} (ncoefs1={}, ncoefs={}, osf={})", latence, ncoefs, h.rows(), osf);
    }

    if(config.debug_actif)
    {
      ArrayXf h = config.wf->filtre.get_coefs(ncoefs, osf);
      auto f = analyse_filtre(h, config.fe);
      f.afficher("Filtre de mise en forme");
    }

    this->config.wf->cnt = 0;

    return 0;
  }

  ArrayXcf flush(int nech)
  {
    if(nech <= 0)
      return ArrayXcf();
    // Nombre d'échantillons d'entrée = 1/osf * nombre d'échantillons de sortie
    return filtre_mise_en_forme->step(ArrayXcf::Zero((int) ceil((1.0 * nech) / osf)));
  }

  ArrayXcf step(const BitStream &bs)
  {
    // TODO : utiliser config.wf->génère_échantillons
    auto x = wf->génère_symboles(bs);
    auto x_symb = x;

    // Filtre de mise en forme, avec sur-échantillonnage intégré
    x = filtre_mise_en_forme->step(x);

    // Fait maintenant directement dans le filtre
    /*if(config.wf->filtre.type == SpecFiltreMiseEnForme::GAUSSIEN)
    {
      auto ma = tsd::filtrage::filtre_mg<cfloat,cdouble>(osf);
      x = ma->step(x);
    }*/

    ArrayXcf x_filtre = x;

    ArrayXf vfreqs, vphase;

    if(config.wf->est_fsk)
    {
      // df = 0.5 * h * fsymb
      // df sur un symbole = 0.5 * h / osf
      // 2 π / osf <=> h = 2
      // => θ = h * 2 * pi / osf / 2 = h * pi / osf

      auto Ω_max = (π * config.wf->index) / osf;
      // h = 2 -> Omega_max = 2 * pi / osf

      //msg("Ω max = {} degrés.", rad2deg(Ω_max));

      // h = 2 fd / fsymb = excursion / fsymb
      // => omega_max = 2%pi*fd = %pi * h * fsymb
      // => Omega_max = omega_max / fs
      //              = %pi * h / ovs;
      // Conversion phase -> IQ

      ArrayXf xr = x.real();
      // normalisation entre [-fmax,fmax]
      vfreqs = xr * (Ω_max / xr.abs().maxCoeff());

      // TODO : échantillon précédent !
      vphase = cumsum(vfreqs);
      x = polar(vphase);
    }

    if(fe1 != config.fe)
    {
      x = ra->step(x);
    }

    // Modulation RF
    if(config.fi != 0)
    {
      x *= ol->step(x.rows());
      if(config.sortie_reelle)
        x = x.real();
    }


    if(config.debug_actif)
    {
      msg("<h2>Modulation</h2>");

      {
        Figures f;
        f.subplot().plot(bs.array(), "|b", "Signal binaire");

        f.subplot().plot(x_symb, "|", "Symboles");

        f.subplot().plot_iq(x_symb, "ob", "Constellation");

        f.subplot().plot(x_filtre, "", "Pulse shaping");

        f.subplot().plot_iq(x_filtre, "ob-", "Constellation filtree");

        f.afficher("Modulation (1)");
      }

      {
        Figures f;
        f.subplot().plot(x_filtre, "", "Signal bande de base");

        f.subplot().plot_psd(x_filtre, fe1);

        if(config.wf->est_fsk)
        {
          f.subplot().plot(vfreqs, "", "Vfreqs");
          f.subplot().plot(vphase, "", "Phase (cumsum)");
        }

        f.subplot().plot(x, "-", "Signal RF");

        f.subplot().plot_psd(x, config.fe);
        f.afficher("Modulation (2)");
      }

    }

    return x;
  }
};


sptr<Modulateur> modulateur_création(const ModConfig &config)
{
  return std::make_shared<ModGen>(config);
}




}




