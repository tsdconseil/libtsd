#include "tsd/tsd-all.hpp"
#include <optional>

#define VERBOSE(AA)

namespace tsd::telecom
{

using namespace std;

// Config du recouvrement d'horloge
struct RecHorlogeConfig
{
  // Forme d'onde
  sptr<FormeOnde> forme_onde;

  // Facteur de sur-échantillonnage
  entier osf;

  // si faux, seule l'interpolation est active (pas de correction)
  bouléen actif = oui;

  // Constante de temps du filtre de boucle
  float tc = 10;

  // Plots de débug
  bouléen debug_actif = non;

  entier ncoefs_filtre_mise_en_forme = 0;

  ItrpType itrp = ItrpType::LINEAIRE;
  entier itrp_lagrange_degré = 1;
};



// Recouvrement d'horloge pour architecture basée sur la décision
struct RecHorloge
{
  // Interpolateur
  sptr<Interpolateur<cfloat>> itrp;

  // Ligne à retard pour l'interpolateur
  Veccf fenetre_x;
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
  entier cnt = 0;

  // Dernière sortie de l'interpolateur
  cfloat lyi = 0.0f;

  bouléen disable = non;

  void reset()
  {
    cnt = 0;

    // On suppose que le démodulateur est démarré au début d'un symbole
    // -> il faut commencer à sampler au milieu, soit 1/2 symbole plus loin
    phase = (config.osf / 2) + 1;

    fenetre_x  = Veccf::zeros(itrp->K);

    // TODO : 5 * config.osf : arbitraire
    // TODO : reset direct de la ligne à retard
    soit z = Veccf::zeros(5 * config.osf);
    fa->step(z);
  }

  void configure(const RecHorlogeConfig &config)
  {
    tsd_assert_msg(config.osf > 0, "Clock rec : OSF invalide ({})", config.osf);

    this->config = config;


    gain = config.osf * lexp_tc_vers_coef(config.tc);

    //itrp = itrp_sinc<cfloat>(7, 0.4, "hn");
    //itrp = itrp_sinc<cfloat>(2*config.osf+1, 0.5, "hn");
    //itrp = itrp_lineaire<cfloat>();

    si(config.itrp == ItrpType::CSPLINE)
      itrp = itrp_cspline<cfloat>();
    sinon si(config.itrp == ItrpType::LINEAIRE)
      itrp = itrp_lineaire<cfloat>();
    sinon si(config.itrp == ItrpType::LAGRANGE)
      itrp = itrp_lagrange<cfloat>(config.itrp_lagrange_degré);
    sinon
    {
      echec("clock rec: itrp inconnu.");
    }

    fa = config.forme_onde->filtre.filtre_adapté(config.ncoefs_filtre_mise_en_forme, config.osf);

    reset();

    msg("rec horloge: osf = {}, npts itrp = {}. phase initiale = {}, tc = {} symboles, gain={}",
        config.osf, itrp->K, phase, config.tc, gain);
  }

  inline void maj_fenetre(Veccf &wnd, cfloat x)
  {
    soit n = wnd.rows();
    wnd.head(n-1) = wnd.tail(n-1);
    wnd(n-1) = x;
  }

  // Filtrage adapté de n échantillons
  inline Veccf step0_filtre_adapte(const Veccf &x)
  {
    retourne fa->step(x);
  }

  // Traite un échantillon, après le filtrage adapté
  // Renvoi vrai si un échantillon est disponible en sortie
  inline std::optional<cfloat> step1_interpolation(cfloat x)
  {
    maj_fenetre(fenetre_x,  x);

    // Requiert: phase >= 1
    phase--;
    si(phase > 1)
      retourne {}; // pas de symbole à sortir

    // Requiert: phase >= 0
    si(phase < 0)
    {
      msg_erreur("clock rec : phase négative ({}). Incrément phase = {}", phase, config.osf);
      phase = 0;
    }

    // Ici on est à la fréquence de la TED
    soit yi  = itrp->step(fenetre_x, 0, phase);

    si(std::isnan(yi.real()) || std::isnan(yi.imag()))
      msg_erreur("Itrp : nan");//, fen itrp = {}.", fenetre_x);

    // Lecture au rythme de la TED, qui travaille à deux fois la fréquence symbole
    phase += ((float) config.osf) / 2;

    si(cnt == 1)
    {
      cnt = 0;
      retourne yi;
    }
    lyi = yi;
    cnt++;
    retourne {}; // Pas de symbole à sortir (on stocke juste yi = valeur intermédiaire entre deux symboles)
  }

  inline float cdot(cfloat a, cfloat b)
  {
    retourne real(a) * real(b) + imag(a) * imag(b);
  }

  inline void step2_maj_retard(cfloat yd0, cfloat yd1)
  {
    si((yd0 != yd1) && (config.actif))
    {
      erreur = cdot(yd1 - yd0, lyi - (yd0 + yd1) / 2.0f) / abs(yd1 - yd0);

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
  entier cnt = 0, cnt1 = 0;

  /** Oscillateur local pour la correction de phase, basé sur une LUT. */
  OLUT lut;

  float gain_cag = 1;

  float alpha_cag = 1;

  sptr<FormeOnde::Ctx> ctx_fo;

  DemodGen2(const ModConfig &modconfig, const DemodConfig &config)
  {
    configure(modconfig, config);
  }

  float delais()
  {
    echec("TODO : DemodGen2::delais()");
    retourne 0;
  }

  void reset(entier cnt)
  {
    this->cnt = cnt;
    cnt1      = osf - 1;
    θ         = 0;
    gain_cag  = 1;
    lye       = 0;
    lf->reset();
    si(config.dec.clock_rec.actif)
      rec_horloge.reset();
    ctx_fo->reset();
  }




  entier configure(const ModConfig &modconfig, const DemodConfig &config)
  {
    this->config    = config;
    this->modconfig = modconfig;
    soit fe     = modconfig.fe,
         fsymb  = modconfig.fsymb,
         fi     = modconfig.fi;
    cnt   = 0;
    osf   = fe / fsymb;
    cnt1  = osf - 1;

    tsd_assert_msg(modconfig.forme_onde, "Démodulateur : la forme d'onde doit être renseignée.");

    msg("Configuration démod: fe={}, fsymb={} (OSF={}), excursion={}",
        fe, fsymb, osf, modconfig.forme_onde->excursion());

    // Configuration du recouvrement d'horloge
    si(config.dec.clock_rec.actif)
    {
      RecHorlogeConfig rhmc;
      rhmc.itrp                         = config.dec.clock_rec.itrp;
      rhmc.itrp_lagrange_degré          = config.dec.clock_rec.itrp_lagrange_degré;
      rhmc.forme_onde                   = modconfig.forme_onde;
      rhmc.osf                          = osf;
      rhmc.debug_actif                  = config.debug_actif;
      rhmc.tc                           = config.dec.clock_rec.tc;
      rhmc.actif                        = config.dec.clock_rec.actif;
      rhmc.ncoefs_filtre_mise_en_forme  = modconfig.ncoefs_filtre_mise_en_forme;
      rec_horloge.configure(rhmc);
    }

    alpha_cag = tsd::filtrage::lexp_tc_vers_coef(config.dec.cag.tc);


    // Initialisation du filtre de boucle pour la correction de phase
    lf = filtre_boucle_ordre_2(config.dec.carrier_rec.BL, config.dec.carrier_rec.η);


    soit reste = fmod(fe, fsymb);
    si(abs(reste) > 1e-6 * fe)
    {
      msg_avert("demod_init: la fréquence d'échantillonnage ({}) doit être un multiple de la fréquence symbole ({}) -- reste = {}.", fe, fsymb, reste);
    }

    // Configuratio de la transposition en bande de base
    si(fi != 0)
    {
      TranspoBBConfig config_tbb;
      config_tbb.fi = fi / fe;
      transpo = transpo_bb<cfloat>(config_tbb);
    }

    ctx_fo = this->modconfig.forme_onde->get_ctx(osf);

    reset(0);



    retourne 0;
  }


  void step(const Veccf &x_, BitStream &bs, Tabf &llr)
  {
    VERBOSE(msg(" demod: start...");)

    si(x_.rows() == 0)
      retourne;

    Veccf x, x_dn, x_crr, x_mf, x_clk, x_agc, x_clk_crr, xf;


    // (1) Transposition en bande de base
    si(modconfig.fi != 0)
      x_dn = transpo->step(x_);
    sinon
      x_dn = x_;

    // Autre solution : démodulation cohérente
    //si(config.wf->est_fsk)
    //{
    //  x_dn = discri->step(x_dn);
    //}

    VERBOSE(msg(" demod: rec horloge...");)

    // Filtrage adapté
    si((osf > 1) && config.dec.clock_rec.actif && config.dec.fa_actif)
      xf = rec_horloge.step0_filtre_adapte(x_dn);
    sinon
      xf = x_dn; // Pas de filtrage possible

    soit n = xf.rows(), ids = 0;

    //////////////////////////////////////////////////////////////////
    // Tableaux ci-dessous : utilisés uniquement pour le débug
    vector<int32_t> si_, idx;
    Veccf v1, v2, v3, v8, v2b, v14;
    Vecf v4, v5, v6, v7, v9, v10, v11, v12, v13;
    si(config.debug_actif)
    {
      // TODO : utiliser un tableau !
      v1  = Veccf::zeros(n);
      v2  = Veccf::zeros(n);
      v2b = Veccf::zeros(n);
      v3  = Veccf::zeros(n);
      v8  = Veccf::zeros(n);
      v14 = Veccf::zeros(n);
      v4  = Vecf::zeros(n);
      v5  = Vecf::zeros(n);
      v6  = Vecf::zeros(n);
      v7  = Vecf::zeros(n);
      v9  = Vecf::zeros(n);
      v10 = Vecf::zeros(n);
      v11 = Vecf::zeros(n);
      v12 = Vecf::zeros(n);
      v13 = Vecf::zeros(n);

    }
    /////////////////////////////////////////////////////////////////

    VERBOSE(msg(" demod: boucle principale ({} échans)...", n);)
    pour(auto i = 0; i < n; i++)
    {
      cnt++;

      // Correction de phase
      cfloat y = xf(i);

      si(config.dec.carrier_rec.actif)
        y *= lut.step(-θ); //* std::polar(1, -θ)

      si(config.debug_actif)
      {
        v14(i) = y;
      }

      // Eventuellement, un étage de CAG
      si(config.dec.cag.actif)
        y *= gain_cag;



      // yi = valeur interpolée (si le bouléen f vaut vrai, sinon pas de valeur)
      cfloat yi;

      si(config.dec.clock_rec.actif)
      {
        soit interp = rec_horloge.step1_interpolation(y);

        si(config.debug_actif)
        {
          v1(i) = y;
          v6(i) = rec_horloge.phase * 100.0 / osf;
          v7(i) = rec_horloge.ph0 * 100.0 / osf;
        }

        // Pas de nouveau symbole, on continue
        si(!interp)
          continue;
        yi = *interp;
      }
      sinon
      {
        si(config.debug_actif)
        {
          v1(i) = y;
          v6(i) = 0;
          v7(i) = 0;
        }

        // TODO: en FSK, il ne faut pas faire ça...
        cnt1 = (cnt1 + 1) % ((entier) osf);

        si(cnt1 == 0)
          yi  = y;
        sinon
          continue;
      }

      si(config.debug_actif)
        v2b(i) = yi;

      // Décision :
      //  - s = index du symbole le plus proche
      //  - ye = symbole I/Q correspondant
      soit [s, ye] = ctx_fo->step(yi);

      // Pas de symbole à sortir
      si(s == -1)
        continue;

      // Enregistre le nouveau symbole décodé
      si_.push_back(s);

      si(config.dec.cag.actif)
      {
        float erreur_gain = abs(yi) / abs(ye);
        si(config.debug_actif)
          v12(i) = erreur_gain - 1;
        gain_cag = (1 - alpha_cag) * gain_cag + alpha_cag * (1/erreur_gain);
      }

      // Maj niveau de bruit
      float bruit = std::norm(ye - yi); // Carré

      // Erreur de phase, basée sur la décision
      float erreur_phase = std::arg(yi * conj(ye));

      // Maj de l'adapteur de rythme
      si(config.dec.clock_rec.actif && (cnt >= 2))
        rec_horloge.step2_maj_retard(lye, ye);

      // Enregistre le dernier symbole décodé
      lye = ye;

      // Boucle de correction de phase
      si((cnt >= 2) && config.dec.carrier_rec.actif)
        θ = lf->step(erreur_phase);

      si(config.debug_actif)
      {
        idx.push_back(i);


        float SNR_lin = abs(ye) / sqrt(bruit);
        float SNR = mag2db(SNR_lin);

        v13(ids)  = SNR;
        v11(ids)  = gain_cag;
        v8(i)     = ye;
        v2(i)     = yi;
        v2b(i)    = yi;
        v4(i)     = erreur_phase;
        v5(i)     = rec_horloge.erreur;
        v9(i)     = rec_horloge.dec;
        v10(ids)  = sqrt(bruit);

        ids++;
      }
    }

    // Conversion index de symbole vers train binaire
    VERBOSE(msg(" demod: symboles -> train binaire...");)
    Veci si2(si_.size());// = Eigen::Map<ArrayXi>(si.data(), si.size());
    memcpy(si2.data(), si_.data(), si_.size() * sizeof(int32_t));
    symdemap_binaire(bs, si2, modconfig.forme_onde->infos.k);
    VERBOSE(msg("nb symboles décodés : {}, nb bits : {}", si_.size(), bs.lon());)


    si(config.debug_actif)
    {
      v10.conservativeResize(ids);
      v11.conservativeResize(ids);
      v13.conservativeResize(ids);

      // Plot transposition
      {
        Figures f;
        f.subplot().plot_iq(x_, ".b", "Signal d'entrée (I/Q)");
        f.subplot().plot(x_);
        f.subplot().plot_psd(x_, modconfig.fe, "", "PSD");

        si(modconfig.fi != 0)
        {
          f.subplot().plot_iq(x_dn, ".b", "Transposition");
          f.subplot().plot(x_dn);
          f.subplot().plot_psd(x_dn, modconfig.fe);// * ratio_ra);
          f.afficher("Démodulation / 1 - transposition");
        }
        sinon
          f.afficher("Démodulation / 1 - entrée");

      }

      // Plot filtrage adapté
      {
        Figures f;
        f.subplot().plot_iq(xf, ".b", "Filtre adapté (const)");
        f.subplot().plot(xf, "-", "Filtre adapté");
        f.afficher("Démodulation / 2 - Filtrage adapté");
      }

      {
        Figures f;
        f.subplot().plot(real(sousvec(v2b, idx)), "a-o", "v2b");
        f.subplot().plot(real(sousvec(v2, idx)), "b-o", "y interpolé");
        f.subplot().plot(real(sousvec(v8, idx)), "g-o", "y décision");
        f.afficher("Démodulation / 2' - Décision");
      }

      // CAG
      {
        Figures f;

        f.subplot().plot(v12, "r-o", "Erreur CAG");
        f.subplot().plot(v11, "b-o", "Gain CAG");
        f.subplot().plot(v13, "b-o", "SNR (dB)");

        f.afficher("Démodulation / 3 - CAG", {1500,1500});
      }

      // Plot correction de phase
      {
        Figures figs;

        soit f = figs.subplot();
        f.plot(sousvec(v4, idx) * (180/π));
        f.titres("Erreur de phase instantanée", "Echantillons", "Erreur (degrés)");

        figs.subplot().plot(v10, "-r", "Erreur (décision - interpol)");

        figs.subplot().plot_iq(v1, ".b", "Corr. de porteuse (const)");

        figs.subplot().plot(v14, "-o", "Corr. de porteuse");
        figs.subplot().plot(v1, "-o", "Corr. de porteuse et CAG");
        figs.afficher("Démodulation / 3 - Correction de porteuse");
      }



      // Plot correction d'horloge
      si(config.dec.clock_rec.actif)
      {
        Figures figs;

        soit f = figs.subplot();
        f.plot(real(sousvec(v2, idx)), "b-o", "y interpolé");
        f.plot(real(sousvec(v8, idx)), "g-o", "y décision");
        //f.plot(sousvec(v3, idx).real(), "m-o", "dy interpolé");

        f = figs.subplot();
        f.plot(sousvec(v5, idx) * 100, "r-o");
        f.titres("Erreur d'horloge instantanée", "Echantillons", "% période symbole");

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
  retourne std::make_shared<DemodGen2>(modconfig, config);
}
}

