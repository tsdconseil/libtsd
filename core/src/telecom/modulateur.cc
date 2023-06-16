// Démodulation générique ((a priori) modulation numérique)

#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"

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

  bouléen valide = non;

  ModConfig config;

  sptr<FiltreGen<cfloat>> filtre_mise_en_forme;

  sptr<FormeOnde> forme_onde;

  void def_forme_onde(sptr<FormeOnde> fo)
  {
    this->forme_onde = fo;
  }

  ModGen(const ModConfig &config)
  {
    configure(config);
    valide = oui;
  }

  float delais() const
  {
    retourne latence;
  }

  void configure(const ModConfig &config)
  {
    latence = 0;
    this->config = config;

    forme_onde = config.forme_onde;

    si(!forme_onde)
      échec("Création modulateur : forme d'onde non spécifiée.");

    msg("<h3>Configuration modulateur</h3>");
    msg("Configuration : fe={} Hz, fi={} Hz, fsymb={} Hz.",
          config.fe, config.fi, config.fsymb);

    fe1 = config.fe;

    bouléen auto_adapation_rythme = non;

    si(auto_adapation_rythme)
    {
      si((fe1 > (16 * config.fsymb)) || (fmod(config.fe,config.fsymb) != 0))
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



    //si(fmod(config.fe,config.fsymb) != 0)
    //{
      //msg_erreur("mod_init: sample frequency must be a multiple of symbol frequency.");
      //retourne -1;
    //}

    entier ncoefs = config.ncoefs_filtre_mise_en_forme;

    filtre_mise_en_forme = config.forme_onde->filtre.filtre_mise_en_forme(ncoefs, osf);
    ol = source_ohc(config.fi / config.fe);

    {
      soit h = config.forme_onde->filtre.get_coefs(ncoefs, osf);

      //assertion(h.rows() == ncoefs);

      // Temps vers le milieu du premier bit transmis
      // Filtre = filtre_rif_ups, R = osf
      //latence = (h.rows() - 1) / 2.0f + osf - 1;
      latence = filtre_rif_ups_délais(h.rows(), osf);
      msg("Modulateur : calcule délais = {} (ncoefs1={}, ncoefs={}, osf={})", latence, ncoefs, h.rows(), osf);
    }

    si(config.debug_actif)
    {
      soit h = config.forme_onde->filtre.get_coefs(ncoefs, osf);
      soit f = plot_filtre(h, non, config.fe);
      f.afficher("Filtre de mise en forme");
    }

    this->config.forme_onde->cnt = 0;
  }

  Veccf flush(entier nech)
  {
    si(nech <= 0)
      retourne {};
    // Nombre d'échantillons d'entrée = 1/osf * nombre d'échantillons de sortie
    //retourne filtre_mise_en_forme->step(ArrayXcf::Zero((entier) ceil((1.0 * nech) / osf)));
    retourne step(Veccf::zeros((entier) ceil((1.0 * nech) / osf)));
  }

  // D'après les données binaires
  Veccf step(const BitStream &bs)
  {
    retourne step(forme_onde->génère_symboles(bs));
  }

    // D'après les symboles
  Veccf step(const Veccf &x_)
  {
    // Filtre de mise en forme, avec sur-échantillonnage intégré
    soit x = filtre_mise_en_forme->step(x_);
    soit x_filtre = x.clone();

    Vecf vfreqs, vphase;

    si(config.forme_onde->infos.est_fsk)
    {
      // df = 0.5 * h * fsymb
      // df sur un symbole = 0.5 * h / osf
      // 2 π / osf <=> h = 2
      // => θ = h * 2 * pi / osf / 2 = h * pi / osf

      soit Ω_max = (π * config.forme_onde->infos.index) / osf;
      // h = 2 -> Omega_max = 2 * pi / osf

      //msg("Ω max = {} degrés.", rad2deg(Ω_max));

      // h = 2 fd / fsymb = excursion / fsymb
      // => omega_max = 2%pi*fd = %pi * h * fsymb
      // => Omega_max = omega_max / fs
      //              = %pi * h / ovs;
      // Conversion phase -> IQ

      soit xr = real(x);

      soit denom = abs(xr).valeur_max();
      si(denom == 0)
        denom = 1.0f;

      // normalisation entre [-fmax,fmax]
      vfreqs = xr * (Ω_max / denom);

      // TODO : échantillon précédent !
      vphase = cumsum(vfreqs);
      x = polar(vphase);
    }

    si(fe1 != config.fe)
    {
      x = ra->step(x);
    }

    // Modulation RF
    si(config.fi != 0)
    {
      x *= ol->step(x.rows());
      si(config.sortie_reelle)
        x = real(x);
    }


    si(config.debug_actif)
    {
      msg("<h2>Modulation</h2>");

      {
        Figures f;
        f.subplot().plot(x_, "|", "Symboles");
        f.subplot().plot_iq(x_, "ob", "Constellation");
        f.subplot().plot(x_filtre, "", "Pulse shaping");
        f.subplot().plot_iq(x_filtre, "ob-", "Constellation filtree");
        f.afficher("Modulation (1)");
      }

      {
        Figures f;
        f.subplot().plot(x_filtre, "", "Signal bande de base");
        f.subplot().plot_psd(x_filtre, fe1);

        si(config.forme_onde->infos.est_fsk)
        {
          f.subplot().plot(vfreqs, "", "Vfreqs");
          f.subplot().plot(vphase, "", "Phase (cumsum)");
        }

        f.subplot().plot(x, "-", "Signal RF");
        f.subplot().plot_psd(x, config.fe);
        f.afficher("Modulation (2)");
      }

    }

    retourne x;
  }
};


sptr<Modulateur> modulateur_création(const ModConfig &config)
{
  retourne make_shared<ModGen>(config);
}




}




