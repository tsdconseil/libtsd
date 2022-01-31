#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"
#include <optional>


#define VERBOSE(AA)

using namespace std;
using namespace tsd::filtrage;
using namespace tsd::vue;

namespace tsd::telecom
{


// Config du recouvrement d'horloge
struct RecHorlogeConfig
{
  // Forme d'onde
  sptr<FormeOnde> forme_onde;

  // Facteur de sur-échantillonnage
  int osf;

  // Si faux, seule l'interpolation est active (pas de correction)
  bool actif = true;

  // Constante de temps du filtre de boucle
  float tc = 10;

  // Plots de débug
  bool debug_actif = false;

  int ncoefs_filtre_mise_en_forme = 0;
};



// Recouvrement d'horloge pour architecture basée sur la décision
struct RecHorloge
{
  // Interpolateur
  sptr<Interpolateur<cfloat>> itrp;

  // Ligne à retard pour l'interpolateur
  ArrayXcf fenetre_x;
  // Phase en cours
  float phase = 0, ph0 = 0;
  // Gain de correction
  float gain = 1.0f;
  // Filtre adapté
  sptr<FiltreGen<cfloat>> fa;
  // Configuration
  RecHorlogeConfig config;



  // Variable de stockage pour les plots
  float erreur, dec;

  // Compteur d'échantillons
  int cnt = 0;

  // Dernière sortie de l'interpolateur
  cfloat lyi = 0.0f;

  bool disable = false;

  void reset()
  {
    cnt = 0;

    // On suppose que le démodulateur est démarré au début d'un symbole
    // -> il faut commencer à sampler au milieu, soit 1/2 symbole plus loin
    phase = (config.osf / 2) + 1;

    fenetre_x  = ArrayXcf::Zero(itrp->K);

    // TODO : 5 * config.osf : arbitraire
    // TODO : reset direct de la ligne à retard
    ArrayXcf z = ArrayXcf::Zero(5 * config.osf);
    fa->step(z);
  }

  void configure(const RecHorlogeConfig &config)
  {
    tsd_assert_msg(config.osf > 0, "Clock rec : OSF invalide ({})", config.osf);

    this->config = config;

    gain = config.osf * (1 - std::exp(-1/(config.tc * config.osf)));

    //itrp = itrp_sinc<cfloat>(7, 0.4, "hn");
    //itrp = itrp_sinc<cfloat>(2*config.osf+1, 0.5, "hn");
    itrp = itrp_lineaire<cfloat>();

    fa = config.forme_onde->filtre.filtre_adapte(config.ncoefs_filtre_mise_en_forme, config.osf);

    reset();

    msg("rec horloge: osf = {}, npts itrp = {}. phase initiale = {}, tc = {}",
        config.osf, itrp->K, phase, config.tc);
  }

  inline void maj_fenetre(ArrayXcf &wnd, cfloat x)
  {
    auto n = wnd.rows();
    wnd.head(n-1) = wnd.tail(n-1).eval();
    wnd(n-1) = x;
  }

  // Filtrage adapté de n échantillons
  inline ArrayXcf step0_filtre_adapte(const ArrayXcf &x)
  {
    return fa->step(x);
  }

  // Traite un échantillon, après le filtrage adapté
  // Renvoi vrai si un échantillon est disponible en sortie
  inline std::optional<cfloat> step1_interpolation(cfloat x)
  {
    maj_fenetre(fenetre_x,  x);

    // Requiert: phase >= 1
    phase--;
    if(phase > 1)
      return {}; // pas de symbole à sortir

    // Requiert: phase >= 0
    if(phase < 0)
    {
      msg_erreur("clock rec : phase négative ({}). Incrément phase = {}", phase, config.osf);
      phase = 0;
    }

    // Ici on est à la fréquence de la TED
    auto yi  = itrp->step(fenetre_x, 0, phase);

    if(std::isnan(yi.real()) || std::isnan(yi.imag()))
      msg_erreur("Itrp : nan, fen itrp = {}.", fenetre_x);

    // Lecture au rythme de la TED, qui travaille à deux fois la fréquence symbole
    phase += ((float) config.osf) / 2;

    if(cnt == 1)
    {
      cnt = 0;
      return yi;
    }
    lyi = yi;
    cnt++;
    return {}; // Pas de symbole à sortir (on stocke juste yi = valeur intermédiaire entre deux symboles)
  }

  inline float cdot(cfloat a, cfloat b)
  {
    return real(a) * real(b) + imag(a) * imag(b);
  }

  inline void step2_maj_retard(cfloat yd0, cfloat yd1)
  {
    if((yd0 != yd1) && (config.actif))
    {
      erreur = cdot(yd1 - yd0, lyi - (yd0 + yd1) / 2.0f) / std::abs(yd1 - yd0);

      // Filtre IIR du premier ordre
      // mu est exprimé en : nombre de samples d'entrée
      // e : en multiple de la période symbole
      dec = gain * erreur;

      // Décalage maximum = 0.25 symboles
      dec = std::clamp(dec, -config.osf/4.0f, config.osf/4.0f);

      phase -= dec;
      ph0   -= dec;
    }
  }

};



struct DemodGen2: Démodulateur
{
  /** Facteur de sur-échantillonnage */
  float osf;

  /** Transposition en bande de base */
  sptr<FiltreGen<cfloat>> transpo;

  /** Configuration en cours */
  DemodConfig config;
  ModConfig modconfig;

  /** Recouvrement d'horloge */
  RecHorloge rec_horloge;

  /** Déphasage en cours */
  float θ = 0;

  /** Filtre de boucle pour la correction de phase */
  sptr<FiltreBoucle> lf;

  /** Symbole précédent (après décision) */
  cfloat lye = 0.0f;

  /** Compteur d'échantillons */
  int cnt = 0, cnt1 = 0;

  /** Oscillateur local pour la correction de phase, basé sur une LUT. */
  OLUT lut;




  sptr<FormeOnde::Ctx> ctx_fo;

  DemodGen2(const ModConfig &modconfig, const DemodConfig &config)
  {
    configure(modconfig, config);
  }

  float delais()
  {
    echec("TODO : DemodGen2::delais()");
    return 0;
  }

  void reset(int cnt)
  {
    this->cnt = cnt;
    cnt1 = osf - 1;
    θ = 0;
    lye   = 0;
    lf->reset();
    if(config.dec.clock_rec.actif)
      rec_horloge.reset();
    ctx_fo->reset();
  }




  int configure(const ModConfig &modconfig, const DemodConfig &config)
  {
    this->config    = config;
    this->modconfig = modconfig;
    auto fe = modconfig.fe, fsymb = modconfig.fsymb, fi = modconfig.fi;
    cnt   = 0;
    osf   = fe / fsymb;
    cnt1  = osf-1;

    tsd_assert_msg(modconfig.wf, "Démodulateur : la forme d'onde doit être renseignée.");

    msg("Configuration démod: fe={}, fsymb={} (OSF={}), excursion={}",
        fe, fsymb, osf, modconfig.wf->excursion());

    // Configuration du recouvrement d'horloge
    if(config.dec.clock_rec.actif)
    {
      RecHorlogeConfig rhmc;
      rhmc.forme_onde   = modconfig.wf;
      rhmc.osf          = osf;
      rhmc.debug_actif  = config.debug_actif;
      rhmc.tc           = config.dec.clock_rec.tc;
      rhmc.actif        = config.dec.clock_rec.actif;
      rhmc.ncoefs_filtre_mise_en_forme = modconfig.ncoefs_filtre_mise_en_forme;
      rec_horloge.configure(rhmc);
    }



    // Initialisation du filtre de boucle pour la correction de phase
    lf = filtre_boucle_ordre_2(config.dec.carrier_rec.BL, config.dec.carrier_rec.η);


    auto reste = fmod(fe, fsymb);
    if(abs(reste) > 1e-6 * fe)
    {
      msg_avert("demod_init: la fréquence d'échantillonnage ({}) doit être un multiple de la fréquence symbole ({}) -- reste = {}.", fe, fsymb, reste);
    }

    // Configuratio de la transposition en bande de base
    if(fi != 0)
    {
      TranspoBBConfig config_tbb;
      config_tbb.fi = fi / fe;
      transpo = transpo_bb<cfloat>(config_tbb);
    }

    ctx_fo = modconfig.wf->get_ctx(osf);

    reset(0);



    return 0;
  }


  /*void regle_horloge(float delais)
  {
    tsd_assert((delais >= 0) && (delais <= 1));
    ArrayXf h = itrp->coefs(delais);
    filtre_itrp = tsd::filtrage::filtre_rif<float,cfloat>(h);
  }*/


  void step(const ArrayXcf &x_, BitStream &bs, ArrayXXf &llr)
  {
    VERBOSE(msg(" demod: start...");)

    if(x_.rows() == 0)
      return;

    ArrayXcf x, x_dn, x_crr, x_mf, x_clk, x_agc, x_clk_crr;


    // (1) Transposition en bande de base
    if(modconfig.fi != 0)
      x_dn = transpo->step(x_);
    else
      x_dn = x_;

    // Autre solution : démodulation cohérente
    //if(config.wf->est_fsk)
    //{
    //  x_dn = discri->step(x_dn);
    //}

    VERBOSE(msg(" demod: rec horloge...");)

    // Filtrage adapté
    ArrayXcf xf;

    if((osf > 1) && (config.dec.clock_rec.actif))
      xf = rec_horloge.step0_filtre_adapte(x_dn);
    else
      xf = x_dn; // Pas de filtrage possible

    int n = xf.rows();

    //////////////////////////////////////////////////////////////////
    // Tableaux ci-dessous : utilisés uniquement pour le débug
    vector<int> si, idx;
    ArrayXcf v1, v2, v3, v8;
    ArrayXf v4, v5, v6, v7, v9, v10;
    if(config.debug_actif)
    {
      v1 = v2 = v3 = v8 = ArrayXcf::Zero(n);
      v4 = v5 = v6 = v7 = v9 = v10 = ArrayXf::Zero(n);
    }
    /////////////////////////////////////////////////////////////////

    VERBOSE(msg(" demod: boucle principale ({} échans)...", n);)
    for(auto i = 0; i < n; i++)
    {
      cnt++;

      // Correction de phase
      cfloat y = xf(i);

      if(config.dec.carrier_rec.actif)
        y *= lut.step(-θ); //* std::polar(1, -θ)

      // yi = valeur interpolée (si le bouléen f vaut vrai, sinon pas de valeur)
      cfloat yi;

      if(config.dec.clock_rec.actif)
      {
        auto interp = rec_horloge.step1_interpolation(y);

        if(config.debug_actif)
        {
          v1(i) = y;
          v6(i) = rec_horloge.phase * 100.0 / osf;
          v7(i) = rec_horloge.ph0 * 100.0 / osf;
        }

        // Pas de nouveau symbole, on continue
        if(!interp)
          continue;
        yi = *interp;
      }
      else
      {
        if(config.debug_actif)
        {
          v1(i) = y;
          v6(i) = 0;
          v7(i) = 0;
        }

        // TODO: en FSK, il ne faut pas faire ça...
        cnt1 = (cnt1 + 1) % ((int) osf);

        if(cnt1 == 0)
          yi  = y;
        else
          continue;
      }

      //config.wf->cnt = cnt-1;

      // Décision :
      //  - s = index du symbole le plus proche
      //  - ye = symbole I/Q correspondant
      auto [s, ye] = ctx_fo->step(yi);

      // Pas de symbole à sortir
      if(s == -1)
        continue;

      // Décision : s = index du symbole le plus proche
      //int s = config.wf->symbole_plus_proche(yi);

      // Enregistre le nouveau symbole décodé
      si.push_back(s);

      // Calcule le symbole I/Q correspondant
      //cfloat ye = config.wf->lis_symbole(s);

      // Maj niveau de bruit
      float bruit = std::norm(ye - yi); // Carré

      // Erreur de phase, basée sur la décision
      float erreur_phase = std::arg(yi * std::conj(ye));

      // Maj de l'adapteur de rythme
      if(config.dec.clock_rec.actif && (cnt >= 2))
        rec_horloge.step2_maj_retard(lye, ye);

      // Enregistre le dernier symbole décodé
      lye = ye;

      // Boucle de correction de phase
      if((cnt >= 2) && config.dec.carrier_rec.actif)
        θ = lf->step(erreur_phase);

      if(config.debug_actif)
      {
        idx.push_back(i);
        v8(i) = ye;
        v2(i) = yi;
        //v3(i) = dyi;
        v4(i) = erreur_phase;
        v5(i) = rec_horloge.erreur;
        v9(i) = rec_horloge.dec;
        v10(i) = std::sqrt(bruit);
      }
    }

    // Conversion index de symbole vers train binaire
    VERBOSE(msg(" demod: symboles -> train binaire...");)
    ArrayXi si2 = Eigen::Map<ArrayXi>(si.data(), si.size());
    symdemap_binaire(bs, si2, modconfig.wf->k);
    VERBOSE(msg("nb symboles décodés : {}, nb bits : {}", si.size(), bs.lon());)


    if(config.debug_actif)
    {
      // Plot transposition
      {
        Figures f;
        f.subplot().plot_iq(x_, ".b", "Signal entree (I/Q)");
        f.subplot().plot(x_);
        f.subplot().plot_psd(x_, modconfig.fe);
        f.subplot().plot_iq(x_dn, ".b", "Transposition");
        f.subplot().plot(x_dn);
        f.subplot().plot_psd(x_dn, modconfig.fe);// * ratio_ra);
        f.afficher("Démodulation / 1 - transposition");
      }

      // Plot filtrage adapté
      {
        Figures f;
        f.subplot().plot_iq(xf, ".b", "Filtre adapté (const)");
        f.subplot().plot(xf, "-", "Filtre adapté");
        f.afficher("Démodulation / 2 - Filtrage adapté");
      }

      // Plot correction de phase
      {
        Figures figs;

        auto f = figs.subplot();
        f.plot((180/π) * sousvec(v4, idx));
        f.titres("Erreur de phase instantannée", "Echantillons", "Erreur (degrés)");

        figs.subplot().plot(v10, "-r", "Erreur (décision - interpol)");

        figs.subplot().plot_iq(v1, ".b", "Corr. de porteuse (const)");

        figs.subplot().plot(v1, "-o", "Corr. de porteuse");
        figs.afficher("Démodulation / 3 - Correction de porteuse");
      }

      // Plot correction d'horloge
      if(config.dec.clock_rec.actif)
      {
        Figures figs;

        auto f = figs.subplot();
        f.plot(sousvec(v2, idx).real(), "b-o", "y interpolé");
        f.plot(sousvec(v8, idx).real(), "g-o", "y décision");
        //f.plot(sousvec(v3, idx).real(), "m-o", "dy interpolé");

        f = figs.subplot();
        f.plot(100.0 * sousvec(v5, idx), "r-o");
        f.titres("Erreur d'horloge instantannée", "Echantillons", "% période symbole");

        //f.subplot();
        //f.plot(100 * v9);
        //f.titres("Corrections de phase", "Echantillons", "% période symbole");

        f = figs.subplot();
        f.plot(v7);
        f.titres("Accu corrections de phase", "Echantillons", "% période symbole");

        f = figs.subplot();
        f.plot(v6);
        f.titres("Phase (wrapping)", "Echantillons", "% période symbole");

        f = figs.subplot();
        f.plot_iq(sousvec(v2, idx), ".b", "Corr. horloge (const)");
        f = figs.subplot();
        f.plot(sousvec(v2, idx), "o-", "Corr. horloge");

        figs.afficher("Démodulation / 4 - Correction d'horloge");
      }

      // Plot démapping
      {
        Figure f;
        f.plot(bs.array(), "|b", "Demapping");
        f.afficher("Démodulation / 5 - Démapping");
      }
      VERBOSE(msg("ok.");)
    }

  }
};


sptr<Démodulateur> demodulateur2(const ModConfig &modconfig, const DemodConfig &config)
{
  return std::make_shared<DemodGen2>(modconfig, config);
}
}

