#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"


#define VERBOSE(AA)

using namespace std;
using namespace tsd::filtrage;
using namespace tsd::vue;

namespace tsd::telecom {



struct DemodGen: Démodulateur
{
  sptr<FiltreGen<cfloat,float>> discri;

  // OSF avant et après transposition
  float osf, osf2, osf3;
  //unsigned int osf;
  sptr<FiltreGen<cfloat>> ra, ra_dp;
  sptr<Filtre<cfloat, cfloat, PLLConfig>> carrier_rec;
  sptr<FiltreGen<cfloat>> psf, dn, clock_rec;
  sptr<FiltreGen<float>> filtre_rssi, filtre_rssi2;
  Ped ped;
  sptr<FiltreGen<cfloat>> filtre_image;

  // Ra après transposition
  float ratio_ra = 1.0;

  // Ra après discri polaire (FSK seulement)
  float ratio_ra_dp = 1.0;

  bouléen valide = non;


  DemodConfig config;
  ModConfig modconfig;

  float delais()
  {
    échec("TODO : DemodGen::delais()");
    retourne 0;
  }

  DemodGen(const ModConfig &modconfig, const DemodConfig &config)
  {
    configure(modconfig, config);
    valide = oui;
  }

  void configure(const ModConfig &modconfig, const DemodConfig &config)
  {
    this->config = config;
    this->modconfig = modconfig;

    soit &fo = modconfig.forme_onde;

    discri = discriminateur_fm();

    si(!fo)
      échec("Démodulateur : la forme d'onde doit être renseignée.");

    si(fmod(modconfig.fe,modconfig.fsymb) != 0)
    {
      msg_avert("Démodulateur : sample frequency must be a multiple of symbol frequency ?");
    }

    this->osf2 = this->osf = modconfig.fe / modconfig.fsymb;

    ratio_ra = 1.0f;

    float ratio_osf_ideal = 8;

    si(modconfig.fe > ratio_osf_ideal * fo->excursion() * modconfig.fsymb)
    {
      ratio_ra = fo->excursion() * modconfig.fsymb * ratio_osf_ideal / modconfig.fe;
      ra = tsd::filtrage::filtre_reechan<cfloat>(ratio_ra);
      //osf2 = ratio_osf_ideal;
      osf2 = ratio_ra * modconfig.fe / modconfig.fsymb;

      msg("Décimation active ({}) après transposition (osf avant = {}, osf après = {}).",
          ratio_ra, osf, osf2);
    }
    osf3 = osf2;

    si(fo->infos.est_fsk && (modconfig.fe * ratio_ra > ratio_osf_ideal * modconfig.fsymb))
    {
      ratio_ra_dp = modconfig.fsymb * ratio_osf_ideal / (modconfig.fe * ratio_ra);
      ra_dp = tsd::filtrage::filtre_reechan<cfloat>(ratio_ra_dp);
      osf3 = ratio_ra * ratio_ra_dp * modconfig.fe / modconfig.fsymb;

      msg("Décimation active ({}) après discri polaire (osf avant = {}, osf après = {}).",
          ratio_ra_dp, osf2, osf3);
    }

    psf = fo->filtre.filtre_adapté(modconfig.ncoefs_filtre_mise_en_forme, osf3);
    TranspoBBConfig config_tbb;
    config_tbb.fi = modconfig.fi / modconfig.fe;
    dn = transpo_bb<cfloat>(config_tbb);

    soit ted  = ted_init(config.ndec.clock_rec.ted);

    assertion(ted);

    sptr<Interpolateur<cfloat>> itrp;
    si(config.ndec.clock_rec.itrp == ItrpType::CSPLINE)
      itrp = itrp_cspline<cfloat>();
    sinon si(config.ndec.clock_rec.itrp == ItrpType::LINEAIRE)
      itrp = itrp_lineaire<cfloat>();
    sinon si(config.ndec.clock_rec.itrp == ItrpType::LAGRANGE)
      itrp = itrp_lagrange<cfloat>(config.ndec.clock_rec.itrp_lagrange_degre);
    sinon
      échec("itrp inconnu.");



    ClockRecConfig crecconfig;
    crecconfig.ted          = ted;
    crecconfig.itrp         = itrp;
    crecconfig.osf          = osf3;
    crecconfig.tc           = config.ndec.clock_rec.tc;
    crecconfig.debug_actif  = config.debug_actif;
    crecconfig.h_fa         = modconfig.forme_onde->filtre.get_coefs(modconfig.ncoefs_filtre_mise_en_forme, osf3);

    msg("crecconfig.tc = {}", crecconfig.tc);

    si(config.ndec.clock_rec.mode_ml)
      clock_rec = clock_rec2_init(crecconfig);
    sinon
      clock_rec = clock_rec_init(crecconfig);

    filtre_rssi = filtre_lexp<float>(lexp_tc_vers_coef(config.ndec.tc_rssi_coarse * osf3)); // 10 symboles
    filtre_rssi2 = filtre_lexp<float>(lexp_tc_vers_coef(config.ndec.tc_rssi_fine)); // 3 symboles


    si(!fo->infos.est_fsk)
    {
      ped = ped_init(config.ndec.carrier_rec.ped, modconfig.forme_onde);


      PLLConfig crr_config;
      crr_config.ped = ped_init(PedType::AUTO, modconfig.forme_onde);
      crr_config.loop_filter_order = 2;
      msg_avert("TODO: config loop filter.");
      //crr_config.
      //crr_config.lf  = lf2_init(config.carrier_rec.BL,config.carrier_rec.η);
      crr_config.debug = config.debug_actif;

      //carrier_rec = carrier_rec_init(crr_config);
      carrier_rec = cpll_création(crr_config);
      msg("Démod init : fe = {} Hz, fi = {} Hz, fsymb = {} Hz.",
          modconfig.fe, modconfig.fi, modconfig.fsymb);
    }

    valide = oui;
  }

  /*template<typename T, typename T1> auto operator()(const std::shared_ptr<T> p)
    {
          delete[] p;
  }*/

  void reset(entier cnt)
  {
    // TODO
  }

  void step(const Veccf &x_, BitStream &bs, Tabf &llr)
  {
    soit &wf = modconfig.forme_onde;

    si(!valide)
      retourne;

    msg(" demod: start...");

    //y.resize(0);

    //ArrayXf coarse_rssi;
    Veccf x, x_dn, x_crr, x_mf, x_clk, x_agc, x_clk_crr;


    // (1) Transposition en bande de base
    si(modconfig.fi != 0)
    {
      msg(" demod: transpo...");
      x_dn = dn->step(x_);
    }
    sinon
    {
      //msg("Downconcersion : fi = 0.");
      x_dn = x_;//x_.eval();
    }


    si(ra)
    {
      msg(" demod: ra...");
      x_dn = ra->step(x_dn);
    }

    ////////////////////////////////
    // Channel filter ?
    // (no, hypothesis: white noise)
    // TODO : devrait être configurable
    ////////////////////////////////

    ///////////////////////////////////////////////
    /// FSK incoherent demodulation
    ///  --> No need pour carrier recovery
    /// But, matched filter in the frequency domain ?
    ///////////////////////////////////////////////

    // (1') Coarse RSSI estimation
    //msg("CRSSI : %d rows.", coarse_rssi.rows());
    //msg("CRSSI : %d rows, min = %f, max = %f", coarse_rssi.rows(), coarse_rssi.minCoeff(), coarse_rssi.maxCoeff());


    // PB SI PED = basé sur la décision :
    // oblige à faire toute la suite échantillon par échantillon...
    // Attention PED QAM : nécessite une étape d'accrochage dans tous les cas
    // Ou alors, d'abord clock rec, puis seulement phase rec.

#   if 0
    // (2) Carrier recovery
    si(!config.forme_onde->est_fsk)
      x_crr = carrier_rec->step(x_dn);
    sinon
    {
      // Incoherent demodulation
      auto n = x_dn.rows();
      x_crr = imag(x_dn.tail(n-1) * conj(x_dn.head(n-1)));
      coarse_rssi = coarse_rssi.head(n-1);
    }
#   endif


    //ArrayXcf x_dn2 = x_dn;

    //assert(!x_dn2.hasNaN());

    // En FSK, pas de carrier rec.
    si(wf->infos.est_fsk) // TODO : tester mode différentiel plutôt
    {
      // Incoherent demodulation
      //auto n = x_dn.rows();
      //si(n == 0)
        //retourne;

      x_dn = discri->step(x_dn);


      //assert(coarse_rssi.rows() == n);
      //coarse_rssi = coarse_rssi.head(n-1).eval();
      //assert(coarse_rssi.rows() == x_dn2.rows());

      assertion(!x_dn.hasNaN());
      si(ra_dp)
      {
        x_dn = ra_dp->step(x_dn);
      }

      assertion(!x_dn.hasNaN());
    }

    //msg(" demod: coarse rssi...");
    //coarse_rssi = filtre_rssi->step(abs(x_dn));

    // Attention ici, la clock rec ML intégre le filtre adapté

    si(config.ndec.clock_rec.mode_ml)
    {
      x_mf = x_dn;

      x_crr = carrier_rec->step(x_mf);

      si(osf3 == 1)
        x_clk = x_crr; // TODO : mueller & mueller devrait pouvoir fonctionner avec un osf de 1
      sinon
        x_clk = clock_rec->step(x_crr);// / (coarse_rssi + 1e-20f));

      x_clk_crr = x_clk;
    }
    sinon
    {
      msg(" demod: psf...");
      // (3) Matched filter
      x_mf = psf->step(x_dn);
      msg(" demod: clkrec...");
      si(osf3 == 1)
        x_clk = x_mf; // TODO : mueller & mueller devrait pouvoir fonctionner avec un osf de 1
      sinon
        x_clk = clock_rec->step(x_mf);// / (coarse_rssi + 1e-20f));

      msg(" demod: crrec...");
      si(!wf->infos.est_fsk) // TODO : tester mode différentiel plutôt
        x_crr = carrier_rec->step(x_clk);
      sinon
        x_crr = x_clk;

      x_clk_crr = x_crr;
    }

    // (3) Clock recovery
    // TODO : ne faire la division que pour le signal d'erreur !
    // --> Il faut pouvoir passer en paramètre le coarse_rssi...


    //msg("clock rec : {} échan entrée -> {} échan sortie.", x_mf.rows(), x_clk.rows());


    // TODO: à mettre ailleurs
    si(wf->infos.est_psk && (wf->infos.M == 4))
    {
      // La carrier rec fait en sorte que la phase soit nulle
      // (ce qui n'est pas le cas de la convention choisie pour QPSK)
      x_crr *= std::polar(1.0f, -π_f/4);
      x_clk_crr *= std::polar(1.0f, -π_f/4);
    }

    msg(" demod: rssi...");
    // (4) Fine RSSI estimation
    Vecf ax = abs(x_clk_crr);
    // PB avec les modulations qui ne sont pas à amplitude constante !
    si(wf->infos.est_psk || wf->infos.est_fsk)
    {
      ax = filtre_rssi2->step(ax);
    }
    sinon
    {
      // TRES BRICOLE !!!!
 //      auto m = ax.mean();
 //      idx = find(ax > m / 2);
 //      ax = mean(ax(idx)) * ones(length(ax),1);
      ax = filtre_rssi2->step(ax); // TEMPORAIRE, TODO
    }

    // (5) AGC
    msg(" demod: agc...");
    x_agc = x_clk_crr / (ax + 1e-10f);

    msg(" demod: demapping...");
    // (6) Demapping
    wf->decode_symboles(bs, x_agc);


    // Ou alors ne pas faire le démapping ici ?

    si(config.debug_actif)
    {
      {
        Figures f;
        f.subplot().plot_iq(x_, ".b", "Signal entree (I/Q)");
        f.subplot().plot(x_);
        f.subplot().plot_psd(x_, modconfig.fe);

        f.subplot().plot_iq(x_dn, ".b", "Downconversion");
        f.subplot().plot(x_dn);
        f.subplot().plot_psd(x_dn, modconfig.fe * ratio_ra);
        f.afficher("Démodulation / 1");
      }

      si(wf->infos.est_fsk)
      {
        Figure f;
        f.plot(x_dn, "", "Discri. polaire, fech={:.1f} Hz, OSF={}",
            modconfig.fsymb * osf3, (entier) osf3);
        f.afficher("Démodulation / fsk");
      }



      {
        Figures f;

        f.subplot().plot_iq(x_mf, "-b", "Matched filter");
        f.subplot().plot(x_mf, "", "matched filter");
        f.subplot().plot_iq(x_clk, ".b", "Clock recovery");
        f.subplot().plot(x_clk, "", "clk rec");
        f.subplot().plot_iq(x_crr, ".b", "After carrier rec");
        f.subplot().plot(x_crr, "", "after carrier rec");
        f.subplot().plot(ax, "-b", "Fine RSSI");
        f.subplot().plot_iq(x_agc, ".b", "AGC");
        f.subplot().plot(x_agc, "", "AGC");
        f.subplot().plot(bs.array(), "|b", "Demapping");
        f.afficher("Démodulation / 2");
      }
      msg("ok.");
    }

  }
};


extern sptr<Démodulateur> demodulateur2(const ModConfig &modconfig, const DemodConfig &config);

sptr<Démodulateur> démodulateur_création(const ModConfig &modconfig, const DemodConfig &config)
{
  si(config.architecture == DemodConfig::ARCHI_SANS_DECISION)
    retourne std::make_shared<DemodGen>(modconfig, config);
  sinon
    retourne demodulateur2(modconfig, config);
}
}
