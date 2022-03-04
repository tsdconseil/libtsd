#include "tsd/telecom.hpp"
#include "tsd/fourier.hpp"
#include "tsd/vue.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/moniteur-cpu.hpp"

#include <algorithm>

using namespace std;
using namespace tsd::filtrage;
using namespace tsd::fourier;
using namespace tsd::vue;

#define VERB(AA)
//AA



namespace tsd::telecom
{


////////////////////////////////////////////////////////////////////////////////
// Flux de données :                                                          //
// step(x) -> tampon -> step_Ne(x : Ne échantillons) -> OLA ou détection RIF  //
//                                                   -> step_demod(x)         //
////////////////////////////////////////////////////////////////////////////////



struct RécepteurImpl: Récepteur
{
  //////////////////////////////////////////////////////////////////
  // Interpolateur passe-tout, pour corriger l'écart d'horloge
  // Calcul des coefficients
  sptr<InterpolateurRIF<cfloat>> itrp;
  // Filtre
  sptr<FiltreGen<cfloat>> filtre_itrp, fa;
  //////////////////////////////////////////////////////////////////

  sptr<FiltreGen<cfloat,float>> discri;

  // Tampon pour avoir des paquets de dimension Ne (--> step_b)
  sptr<SinkGen<cfloat>> tampon;

  // Configuration
  RécepteurConfig config;

  // Démodulateur cohérent
  sptr<Démodulateur> demod;

  // Détecteur d'en-tête
  sptr<Detecteur> détecteur;

  sptr<FiltreGen<cfloat>> decim;

  sptr<FormeOnde> wf;

  float osf;

  bool demod_en_cours = false;

  int dim_motif = 0;

  vector<Detection> fingers;

  RécepteurTrame trame; // Trame en cours de réception

  int Ne = 0;

  ArrayXcf last;

  // Compteur d'échantillons, remis à zéro à chaque appel externe de la fonction step,
  // et incrémenté de Ne échantillons à chaque appel interne de step_Ne.
  int cnt_x = 0;

  vector<RécepteurTrame> trames;

  bool tampon_vide = true;

  // Nombre de bits / symboles
  int nbits_par_symbole_entete, nbits_par_symbole_données;

  MoniteurCpu mon_ola{"recepteur/ola"},
              mon_demod{"recepteur/demod"},
              mon_misc{"recepteur/misc"};

  // Délais interpolateur + filtre adapté
  int delais_interpolateur;

  // Nombre d'échantillons supplémentaires en début de motif
  // (vaut 0 ou 0.5, suivant si le délais du modulateur est entier ou non)
  // Eg. si δt_modulateur == 0, le premier échantillon du motif est le milieu du premier symbole
  // sinon 0.5 symboles avant.
  // TODO : prendre plutôt un échantillon avant
  float δt_modulateur = 0;

  float δt_interpolateur = 0;

  MoniteursStats moniteurs()
  {
    auto res = détecteur->moniteurs();
    res.ajoute(mon_ola);
    res.ajoute(mon_demod);
    res.ajoute(mon_misc);
    return res;
  }


  // Fonction avec          OSF=3, filtre NRZ
  // Ne fonctionne pas avec OSF=3, filtre SRRC


  std::tuple<int,float> calc_retard()
  {
    ArrayXf cf = config.format.modulation.wf->filtre.get_coefs(
        config.format.modulation.ncoefs_filtre_mise_en_forme, osf);

    // osf/2 pour aller au milieu du premier bit
    float d = rif_delais(config.ncoefs_interpolateur)
            + rif_delais(cf.rows())
            + 1
            + osf / 2.0f;


    int     delais_interpolateur = floor(d);
    float δt_interpolateur   = d - delais_interpolateur;

    return {delais_interpolateur, δt_interpolateur};
  }

  void regle_delais(float δ)
  {
    //fa = config.format.modulation.wf->filtre.filtre_adapte(osf);

    δ = 1 - δ;

    VERB(msg("Régle délais : δ={}", δ);)

    tsd_assert(itrp);

    δ = std::clamp(δ, 0.0f, 1.0f);

    ArrayXf h1 = itrp->coefs(δ);
    ArrayXf h2 = config.format.modulation.wf->filtre.get_coefs(
        config.format.modulation.ncoefs_filtre_mise_en_forme, osf);

    if(config.debug_actif)
    {
      analyse_filtre(h2, 1.0f).afficher("Filtre adapté");
      msg("Création d'un filtre d'interpolation, délais = {} échantillons ({} symboles)", δ, δ / osf);
      analyse_filtre(h1, 1.0f).afficher("Filtre d'interpolation");
    }

    filtre_itrp = filtre_rif<float,cfloat>(h1);
    fa          = filtre_rif<float,cfloat>(h2);

    std::tie(delais_interpolateur, δt_interpolateur) = calc_retard();

    VERB(msg("Nb coefs : interpolateur = {}, fa = {}", h1.rows(), h2.rows());)
    VERB(msg("Délais interpolateur: {} + {} échans.", delais_interpolateur, δt_interpolateur));
  }


  RécepteurEtat get_etat()
  {
    RécepteurEtat res;

    res.Ne = Ne;

    return res;
  }



  int configure(const RécepteurConfig &rc)
  {
    config = rc;

    const auto &conf_mod = config.format.modulation;

    this->wf    = conf_mod.wf;
    auto fe     = conf_mod.fe;
    auto fsymb  = conf_mod.fsymb;

    auto fo_entete = config.format.modulation.wf;

    if(config.format.fo_entete)
    {
      fo_entete = config.format.fo_entete;
      msg_majeur("Récepteur: fo entete = {}",  *fo_entete);
      msg_majeur("Récepteur: fo données = {}", *wf);
    }

    if((fe <= 0)
        || (fsymb <= 0)
        || (modulo(fe, fsymb) != 0))
      echec("Récepteur : fréquences invalides (fe={} Hz, fsymb={} Hz).", fe, fsymb);

    if(wf->infos.est_fsk)
    {
      discri = discriminateur_fm();
    }

    demod_en_cours = false;
    mon_ola.reset();
    mon_demod.reset();
    fingers.clear();
    tampon_vide = true;
    trames.clear();
    trame.bs.clear();
    last.resize(0);

    //ModConfig config_mod;

    // Pas besoin de recouvrement d'horloge,
    // on est déjà calé
    /*config.config_demod.dec.clock_rec.actif = false;
    {
      // La décimation à 1SPS est faite ici
      config.config_demod.fsymb = fsymb;
      config.config_demod.fe    = fsymb;
      config.config_demod.fi    = 0;
    }*/
    // En FSK, pas de décimation ici
    //if(wf->est_fsk)
    //  config.config_demod.fe  = config.fe;

    ModConfig config_mod_entete = conf_mod;
    config_mod_entete.wf = fo_entete;
    auto mod = modulateur_création(config_mod_entete);


    ModConfig config_mod_données = conf_mod;
    config_mod_données.fsymb  = fsymb;
    config_mod_données.fe     = fsymb;
    config_mod_données.fi     = 0;
    osf = conf_mod.fe / conf_mod.fsymb; // Toujours 1 !
    config.config_demod.dec.clock_rec.actif = false;
    demod = démodulateur_création(config_mod_données, config.config_demod);

    // Nombre de bits / symboles
    nbits_par_symbole_données = wf->infos.k;
    nbits_par_symbole_entete  = fo_entete->infos.k;
    float df  = mod->delais();

    msg_majeur("Récepteur / modulateur : df = {}", df);

    // df est le délais vers le milieu du premier bit transmis.
    // il faut donc démarrer 1/2 échantillon avant pour avoir aussi le début du bit
    df -= osf / 2.0f;

    // Pb si NRZ, OSF = 2 : df = - 0.5 !!!!

    // En effet, df(premier) = 0.5

    int   di  = (int) floor(df);
    δt_modulateur = df - di;
    msg_majeur("Récepteur / modulateur : df = {}, di = {}", df, di);

    BitStream entete = config.format.entete;
    entete.pad_mult(fo_entete->infos.k);

    // Pour flusher le modulateur
    BitStream et2 = entete;
    for(auto i = 0; i < (di + 8) * nbits_par_symbole_entete; i++)
      et2.push(0);

    ArrayXcf motif = mod->step(et2);
    auto e = motif.abs2().mean();

    // Peut arriver avec un OSF de 2, et filtre NRZ :
    // dans ce cas la latence du modulateur est de 0.5, et le premier bit commence à -0.5
    if(di < 0)
    {
      motif = vconcat(ArrayXcf::Zero(1), motif);
      di++;
      df++;
    }

    tsd_assert(di >= 0);

    if(wf->infos.est_fsk)
    {
      motif = discri->step(motif);
      // Conversion angle -> valeur entre -1 et 1
      motif /= (2 * π * wf->excursion() / 2) / (config_mod_entete.fe / config_mod_entete.fsymb);
    }

    float fc = 0.5;
    if(osf > 1)
      // Moyenne entre 0.5 et 0.5 / osf
      fc = 0.45;//0.5 * (0.5 / osf + 0.5);

    VERB(msg("Récepteur : fréquence de coupure interpolateur : {} (osf = {})", fc, osf);)

    itrp = tsd::filtrage::itrp_sinc<cfloat>({config.ncoefs_interpolateur, 1024, fc, "hn"});


    int nsymbs = (entete.lon() + nbits_par_symbole_entete - 1) / nbits_par_symbole_entete;
    // M = nombre d'échantillons de l'en-tête
    int M = nsymbs * osf;

    VERB(msg("Récepteur : {} bits / symbole, osf={}, nsymbs={} -> M = {} échans.", nbits_par_symbole_entete, osf, nsymbs, M);)


    DetecteurConfig config_detecteur;

    float C;
    int Nf, Nz;
    ola_complexité_optimise(M, C, Nf, Nz, Ne);

    msg("Récepteur : calcul auto Ne optimal : M={} --> Ne={},Nz={},Nf=Ne+Nz={},C={} FLOPS/ech", M, Ne, Nz, Nf, C);
    config_detecteur.Ne = Ne;
    config_detecteur.calculer_signal_correlation = config.callback_corr ? true : false;
    config_detecteur.mode = DetecteurConfig::MODE_OLA;
    //ola_config.mode = DetecteurConfig::MODE_RIF;



    // if(config_mod.wf->M != 2)
      // echec("TODO : récepteur / gestion délais M != 2");

    msg("Latence modulateur : {} (floor: {}), osf : {}", df, di, osf);
    msg("Longueur d'en-tête avec flush : {} échans.", motif.rows());
    msg("Longueur à extraire : {} échans @ {}.", M, di);
    msg("Découpage en blocs de {} échan.", Ne);

    tsd_assert_msg(di + M <= motif.rows(),
        "Latence modulateur : di = {}, nb échantillons théo en-tête : M = {}, nb échan générés : {}", di, M, motif.rows());

    config_detecteur.motif  = motif.segment(di, M);

    if(config.debug_actif)
    {
      Figures f;
      f.subplot().plot(config.format.entete.array(), "|b", "En-tête binaire");
      f.subplot().plot(motif);
      f.gcf().titre("En-tête modulé (complet)");
      f.subplot().plot(config_detecteur.motif);
      f.gcf().titre("En-tête modulé (segment)");
      f.afficher(fmt::format("Récepteur : en-tête détection (e = {}).", e));
    }

    if(config.debug_actif)
    {
      auto [lags, xc] = xcorrb(config_detecteur.motif, config_detecteur.motif);
      Figure f;
      f.plot(lags, xc.abs(), "", "Auto-corrélation motif (biaisée)");
      f.afficher();
    }


    // TODO: renommer en M
    dim_motif               = config_detecteur.motif.rows();
    config_detecteur.seuil  = config.seuil;
    config_detecteur.gere_detection = [&](const Detection &det)
      {
        if(det.SNR_dB >= config.SNR_mini)
        {
          VERB(msg("Récepteur, {}", det);)
          fingers.push_back(det);
        }
      };

    config_detecteur.debug_actif = config.debug_actif;
    détecteur = détecteur_création(config_detecteur);

    tampon = tampon_création<cfloat>(Ne,
        [&](const Vecteur<cfloat> &x)
        {
          step_Ne(x);
        });

    last = ArrayXcf::Zero(Ne);
    return 0;
  }








  vector<RécepteurTrame> step(const ArrayXcf &x)
  {
    VERB(msg("récepteur : step({} échans)", x.rows()));
    cnt_x = 0;
    if(tampon_vide && (x.rows() == Ne))
    {
      step_Ne(x);
    }
    else
    {
      tampon_vide = false;
      tampon->step(x);
    }
    auto res = trames;
    trames.clear();
    return res;
  }

  // A voir si on peut pas s'en passer (du fait que la dimension doive être Ne)
  void step_Ne(const ArrayXcf &x_)
  {
    tsd_assert(x_.rows() == Ne);

    VERB(msg("Récepteur : step Ne={}, cnt_x = {}", Ne, cnt_x));

    ArrayXcf x;

    if(wf->infos.est_fsk)
    {
      // Ceci introduit un retard de 1 échantillon
      x = discri->step(x_);

      if(config.debug_actif)
      {
        Figures f;
        f.subplot().plot(x_, "", "Avant discri FM.");
        f.subplot().plot(x, "", "Après discri FM.");
        f.afficher("Récepteur : discri FM");
      }

    }
    else
      x = x_;

    fingers.clear();

    mon_ola.commence_op();
    ArrayXf corr = détecteur->step(x);
    mon_ola.fin_op();

    // idx dans [0,Ne-1] - delais_corr
    // Avec delais_corr = Ne  en mode FFT
    //                    M-1 en mode RIF

    if(config.callback_corr)
      config.callback_corr(corr);

    if(demod_en_cours)
    {
      step_demod(x);
      if(trame.bs.lon() == config.format.nbits)
      {
        trames.push_back(trame);
        trame.bs.clear();
      }
    }

    mon_misc.commence_op();

    for(auto &finger: fingers)
    {
      trame.det = finger;

      trame.bs.clear();
      trame.x  = ArrayXcf();
      trame.x1 = ArrayXcf();
      cnt_demod = 0;

      // Position relative au buffer appelant
      trame.det.position      += cnt_x;
      trame.det.position_prec += cnt_x;
      trame.det.position_prec += δt_modulateur;

      VERB(msg("Récepteur, après offset cnt (cnt_x={}): {}", cnt_x, trame.det));

      auto fo_entete = config.format.modulation.wf;
      if(config.format.fo_entete)
        fo_entete = config.format.fo_entete;

      int nb_bits_par_symbole = fo_entete->infos.k;
      int nb_symb_entete = config.format.entete.lon() / nb_bits_par_symbole;

      trame.EbN0 = pow2db((db2pow(finger.SNR_dB) * osf) / nb_bits_par_symbole);

      demod_en_cours = true;

      VERB(msg("nb_bits_par_symbole = {}, nb bits en tete = {}, nb_symb_entete = {}",
          nb_bits_par_symbole, config.format.entete.lon(), nb_symb_entete);)

      // Démarrage du démodulateur en lui disant le nombre de symboles dans l'en-tête
      // (nécessaire pour avoir l'état adéquat en π/4-QPSK)
      demod->reset(nb_symb_entete);
      decim = decimateur<cfloat>(osf);

      // Reset des filtres
      if(filtre_itrp)
      {
        ArrayXcf z = ArrayXcf::Zero(config.ncoefs_interpolateur);
        filtre_itrp->step(z);
        fa->step(z);
      }

      tsd_assert(abs(finger.position_prec - finger.position) < 1);

      // Idem regle_delais()
      std::tie(delais_interpolateur, δt_interpolateur) = calc_retard();

      //msg("Récepteur : pos avant: prec={}, int={}", finger.position_prec, finger.position);
      finger.position_prec += δt_modulateur;
      finger.position_prec += δt_interpolateur;
      while(finger.position_prec < finger.position)
        finger.position--;
      while(finger.position_prec >= finger.position + 1)
        finger.position++;

      auto δ = finger.position_prec - finger.position;

      //msg("Récepteur : pos après: prec={}, int={}", finger.position_prec, finger.position);
      VERB(msg_majeur("Récepteur : réglage délais : δt_modulateur={}, δt_interpolateur={}.", δt_modulateur, δt_interpolateur);)

      //regle_delais(1-δ);
      regle_delais(δ);

      // Index à partir du buffer précédent

      // finger.position_prec = index à partir du buffer en cours
      // (peut être négatif pour indiquer un début sur la fin du buffer précédent)

      int nspl_strict = (config.format.nbits * osf + nbits_par_symbole_données-1) / nbits_par_symbole_données;
      // + de quoi flusher les filtres
      int nspl = nspl_strict + 4 * osf + delais_interpolateur;

      // Entre 0 et Ne-1
      // Position relative au début du tampon précédent
      int idx = finger.position + Ne;

      // Saute le motif
      idx += dim_motif;

      // A cause du discriminateur, qui induit un délais de 2 échantillons
      /*if(wf->est_fsk)
      {
        idx += 2;
        delais_interpolateur -= 2;
      }*/

      //msg("nspl strict = {} (nbits = {}), idx={}, Ne={}, osf={}, nbits_par_symbole_données={}", nspl_strict, config.format.nbits, idx, Ne, osf, nbits_par_symbole_données);

      ArrayXcf y;

      // Tout sur le deuxième tampon ?
      if(idx >= Ne)
      {
        if(nspl > Ne - (idx - Ne))
          nspl = Ne - (idx - Ne);

        tsd_assert((idx - Ne >= 0) && (idx - Ne + nspl <= Ne));

        y = x.segment(idx - Ne, nspl);
      }
      else
      {
        // Nombre d'échantillons à prendre sur le tampon précédent
        int nspl1 = min(nspl, Ne - idx);

        // Tous les échantillons peuvent être pris du tampon précédent
        if(nspl1 == nspl)
        {
          tsd_assert((idx >= 0) && (idx + nspl <= Ne));
          y = last.segment(idx, nspl);
        }
        else
        {
          // Check débordement tampon suivant
          if(nspl - nspl1 > Ne)
            nspl = Ne + nspl1;

          tsd_assert(nspl > 0);
          tsd_assert(nspl1 > 0);
          tsd_assert(nspl1 <= Ne);
          tsd_assert((nspl >= nspl1) && (nspl - nspl1 <= Ne));

          y.resize(nspl);
          y.head(nspl1) = last.tail(nspl1);
          y.tail(nspl - nspl1) = x.head(nspl - nspl1);
        }
      }

      // int demod_start = finger.position_prec + dim_motif + x.rows();
      // tsd_assert(demod_start >= 0);
      // ArrayXcf y = last.tail(last.rows() - demod_start);
      // y = tsd::vconcat(y, x);
      // TODO : ajouter la suite à partir de x

      if(config.debug_actif && (last.rows() > 0))
      {
        Figures fs;
        auto f = fs.subplot();
        f.plot(last);
        f.plot((float) (idx - dim_motif), 0.0f, "ro");
        f.titre("Buffer complet");

        f = fs.subplot();
        f.plot(last.segment(idx-dim_motif, std::min(nspl_strict + dim_motif, Ne - (idx-dim_motif))));

        if(idx - dim_motif + 1 < last.rows())
        {
          f.plot((float) (0), real(last(idx-dim_motif)), "ro");
          f.plot((float) (δ), real(last((int) round(idx-dim_motif+δ))), "bs");
        }
        if(idx + 1 < last.rows())
        {
          f.plot((float) (dim_motif), real(last(idx)), "mo");
          f.plot((float) (dim_motif+δ), real(last((int) round(idx+δ))), "bs");
        }
        f.titre("Zoom");


        f = fs.subplot();
        f.plot(nspl_strict >= y.rows() ? y : y.head(nspl_strict));

        //f.plot(last.segment(finger.position + x.rows(), 500));
        f.titre("Sous-buffer : données à démoduler");
        fs.afficher("Données reçues");
      }

      // Il faut ajouter en tête de y de quoi faire démarrer le filtre d'interpolation

      mon_misc.fin_op();

      step_demod(y);


      mon_misc.commence_op();

      if(trame.bs.lon() == config.format.nbits)
      {
        trames.push_back(trame);
        trame.bs.clear();
      }
      else if(&finger != &(fingers.back()))
      {
        msg_avert("Pas de quoi finir de démoduler dans ce buffer, il manque {} bits / {} (un deuxième finger écrase le précédent).",
            config.format.nbits - trame.bs.lon(), config.format.nbits);
        msg("  score finger en cours : {} (SNR = {:.1f} dB), dernier finger : {} (SNR = {:.1f} dB)",
            finger.score, finger.SNR_dB, fingers.back().score, fingers.back().SNR_dB);
        int nspl_strict = (config.format.nbits * osf + nbits_par_symbole_données-1) / nbits_par_symbole_données;
        msg("  néchantillons dispos : {}, nbéchantillons nécessaires : {} (nbits={}, osf={}, nbits_par_symbole={} bits/symbole).", Ne - idx, nspl_strict, config.format.nbits, osf, nbits_par_symbole_données);
      }
    }

    last = x;
    cnt_x += Ne;

    mon_misc.fin_op();
  }

  int cnt_demod = 0;
  void step_demod(const ArrayXcf &x)
  {
    mon_demod.commence_op();

    BitStream bs;
    ArrayXXf llr;

    trame.x = vconcat(trame.x, x);

    int nb_echans_theo = osf * (config.format.nbits + nbits_par_symbole_données - 1) / nbits_par_symbole_données;

    if(trame.x.rows() > nb_echans_theo)
      trame.x = trame.x.head(nb_echans_theo).eval();

    ArrayXcf x1 = x * std::polar(1.0f / trame.det.gain, -trame.det.θ);

    // Le premier échantillon de x1 est exactement le milieu du premier symbole

    ArrayXcf y = filtre_itrp->step(x1);
    ArrayXcf y1 = fa->step(y);

    // TODO : intégrer ensemble, dans un seul filtre polyphase,
    // l'interpolateur, le filtre adapté et la décimation.

    // osf/2 pour ce placer au milieu du premier symbole
    ArrayXcf y2, y3;

    if(delais_interpolateur == 0)
      y2 = y1;
    else if(delais_interpolateur >= y1.rows())
    {
      delais_interpolateur -= y1.rows(); // ???
    }
    else
    {

      VERB(msg("Recalage signal : délais = {}", delais_interpolateur));

      y2 = y1.tail(y1.rows() - delais_interpolateur);
      delais_interpolateur = 0; // ???
    }

    // Maintenant, on ne garde que une échantillon tous les R
    //  sauf si la waveform en veut plus
    //if(wf->est_fsk)
    //  y3 = y2;
    //else
    y3 = decim->step(y2);


    //msg("TX1/a: {} elms, min coeff = {}", trame.x1.rows(), trame.x1.rows() > 0 ? trame.x1.abs().minCoeff() : 0);

    trame.x1 = vconcat(trame.x1, y3).eval();


    //msg("TX1/b: {} elms, min coeff = {}", trame.x1.rows(), trame.x1.rows() > 0 ? trame.x1.abs().minCoeff() : 0);

    if(trame.x1.rows() > nb_echans_theo / osf)
      trame.x1 = trame.x1.head(nb_echans_theo / osf).eval();


    /*if(trame.x1.abs().minCoeff() > 1e10)
    {
      msg_avert("PB REC: det={}", trame.det);
      msg("x.mincoef={},max={}", x.abs().minCoeff(),x.abs().maxCoeff());
      msg("x1.mincoef={},max={}", x1.abs().minCoeff(),x1.abs().maxCoeff());
      msg("y: {}, {}", y.abs().minCoeff(), y.abs().maxCoeff());
      msg("y1: {}, {}", y1.abs().minCoeff(), y1.abs().maxCoeff());
      msg("y2: {}, {}", y2.abs().minCoeff(), y2.abs().maxCoeff());
      msg("y3: {} elmts, {}, {}", y3.rows(), y3.abs().minCoeff(), y3.abs().maxCoeff());
      msg("trame.x1: {}, {}", trame.x1.abs().minCoeff(), trame.x1.abs().maxCoeff());
    }*/

    if(config.debug_actif)
    {

      //msg("Correction de porteuse : phase={}°, gain={:.2f}", rad2deg(trame.det.θ), trame.det.gain);
      Figures fs;
      fs.subplot().plot(x, "", "Signal à démoduler (x)");
      fs.subplot().plot(x1, "", "Correction de gain et phase (x1)");
      fs.subplot().plot(y, "", "Après interpolateur (y)");
      fs.subplot().plot(y1, "", "Après filtre adapté (y1)");
      fs.subplot().plot(y2, "", "Décalage début premier symbole (y2)");
      fs.subplot().plot(y3, "", "Décimation (y3)");
      fs.afficher("Interpolation");

      {
        Figure f;
        f.plot(y2);
        for(auto i = 0; i < y2.rows(); i += osf)
          f.plot((float) i, real(y2(i)), "ro");
        f.afficher("Timing");
      }

      {
        Figure f;
        ArrayXcf y3((int) floor(y2.rows() / osf));
        for(auto i = 0; i < y2.rows(); i += osf)
          if(i / osf < y3.rows())
            y3(i / osf) = y2(i);
        f.plot_iq(y3, "bo");
        f.afficher("Constellation après timing");
      }
    }


    demod->step(y3, bs, llr);

    // TODO: optim
    for(auto i = 0; (i < bs.lon()) && (trame.bs.lon() < config.format.nbits); i++)
      trame.bs.push(bs[i]);

    if(trame.bs.lon() == config.format.nbits)
    {
      demod_en_cours = false;
    }

    mon_demod.fin_op();
  }

};

sptr<Récepteur> récepteur_création(const RécepteurConfig &rc)
{
  auto res = make_shared<RécepteurImpl>();
  res->configure(rc);
  return res;
}



}
