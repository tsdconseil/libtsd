#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include "tsd/telecom.hpp"
#include "tsd/telecom/carrier-rec.hpp"
#include <iostream>


// TOUT CA (sauf PLL peut-être) doit être invisible pour les TP

using namespace tsd::vue;

namespace tsd::telecom {

using namespace tsd::filtrage;


struct FMDiscri: FiltreGen<cfloat, float>
{
  entier cnt = 0;
  cfloat x0 = 0.0f, x1 = 0.0f, x2 = 0.0f;
  cfloat last = 0.0f;

  /*void configure_impl(const FMDiscriConfig &c)
  {
    //config = c;
    x0 = x1 = x2 = 0.0f;
  }*/
  void step(const Vecteur<cfloat> &x, Vecteur<float> &y)
  {
    // D'après https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
    entier n = x.rows();

    y.resize(n);

    // PB si x0 = x1 = 0, x2 != 0

    pour(auto i = 0; i < n; i++)
    {
      x2 = x1;
      x1 = x0;
      x0 = x(i);

      si((cnt <= 2) || (x1 == 0.0f)) // TODO : test de nullité bancal
        y(i) = 0;
      sinon
        y(i) = (real(x1) * (imag(x0) - imag(x2)) - imag(x1) * (real(x0) - real(x2))) / (2 * std::norm(x1));

      y(i) = std::clamp(y(i), -π_f, π_f);

      cnt++;
    }

    // Attention, le discri génère une forme de filtrage !!!

#   if 1
    pour(auto i = 0; i < n; i++)
    {
      // Discrimination polaire
      y(i) = std::arg(x(i) * conj(last));
      last = x(i);
    }
#   endif
  }
};



sptr<FiltreGen<cfloat,float>> discriminateur_fm()
{
  retourne std::make_shared<FMDiscri>();
}





struct ModulateurAM: Filtre<float, float, AMConfig>
{
  sptr<FiltreGen<float>> ra;
  sptr<SourceGen<cfloat>> ol;
  sptr<Filtre<float, cfloat, HilbertTransformeurConfig>> hilbert;

  ModulateurAM()
  {
    configure(config);
  }
  void configure_impl(const AMConfig &c)
  {
    ra = filtre_reechan<float>(c.fe_rf / c.fe_sortie);
    ol = source_ohc(c.f_rf / c.fe_rf);

    si(c.est_BLU())
    {
      hilbert = hilbert_transformeur();
      hilbert->configure({255, "hn"});
    }
  }
  void step(const Vecteur<float> &x1, Vecteur<float> &y)
  {
    soit &config = Configurable<AMConfig>::config;

    // (1) Sur-échantillonnage
    ra->step(x1, y);

    // https://www.mathworks.com/help/signal/examples/single-sideband-modulation-via-the-hilbert-transform.html


    si(config.mode == AMConfig::Mode::DSB)
    {
      soit mxcoef = abs(y).valeur_max();

      // (2) Transposition
      y = 1 + config.indice * y / mxcoef;
      y *= real(ol->step(y.rows()));
    }
    sinon si (config.mode == AMConfig::Mode::DSB_SUPPRESSED_CARRIER)
    {
      y *= real(ol->step(y.rows()));
    }
    sinon si(config.est_BLU())
    {
      // https://www.mathworks.com/help/signal/examples/single-sideband-modulation-via-the-hilbert-transform.html
      // Calcul du signal analytique
      soit z    = hilbert->step(y);
      soit xol  = ol->step(z.rows());
      float signe = config.mode == AMConfig::Mode::LSB ? -1 : 1;
      y = real(z) * real(xol) + signe * imag(z) * imag(xol);
    }


    si(config.debug_actif)
    {
      Figures fig;
      fig.subplot().plot(x1, "", "Signal audio");
      fig.subplot().plot_psd(x1, config.fe_sortie);
      fig.subplot().plot(y, "", "Signal AM");
      fig.subplot().plot_psd(y, config.fe_rf);
      fig.afficher(format("Modulateur AM - frf = {:.1f} Hz", config.f_rf));
    }


  }

};


struct DemodulateurAM: Filtre<cfloat, float, AMConfig>
{
  sptr<SourceGen<cfloat>> ol;
  // pour modes LSB, USB
  sptr<FiltreGen<float>> retard, hilbert, ra, filtre_lp, filtre_hp;

  DemodulateurAM()
  {
    configure(Configurable<AMConfig>::config);
  }
  void configure_impl(const AMConfig &c)
  {
    soit &config = Configurable<AMConfig>::config;
    config = c;
    ra = filtre_reechan<float>(config.fe_sortie / config.fe_rf);
    ol = source_ohc(-config.f_rf / config.fe_rf);

    soit fc_high = config.fcut_audio_high / config.fe_rf;
    soit fc_low  = config.fcut_audio_low / config.fe_rf;

    soit hlp = design_riia(4, "pb", "butt", fc_high);
    soit hhp = design_riia(4, "pb", "butt", fc_low);
        //butterworth_lp(fc_high, 2);
    //auto hhp = butterworth_lp(fc_low,  2);

    msg("FC low = {}, FC high = {}", fc_low, fc_high);

    //msg("Coefs LP :\n {}\n {}", hlp.numer, hlp.denom);

    //analyse_filtre(hlp, config.fe_rf);
    //analyse_filtre(hhp, config.fe_rf);

    filtre_lp = filtre_sois<float>(hlp);
    filtre_hp = filtre_sois<float>(hhp);

    si(config.est_BLU())
    {
      hilbert = filtre_rif<float>(design_rif_hilbert(127));
      retard = ligne_a_retard<float>(127/2);
    }

    //fech_sortie = config.fe_rf;
  }
  void step(const Veccf &x, Vecf &y)
  {
    soit &config = Configurable<AMConfig>::config;
    soit n = x.rows();
    soit xol = ol->step(n);
    soit rftc = x * xol;

    Vecf y1;


    // DSB ou DSB sans porteuse
    si(!config.est_BLU())
    {
      // Idem méthode Wikipédia,
      // mais modifiée pour ne pas être dépendant de la phase
      // relative de l'OL par rapport à la porteuse
      soit rft = real(rftc);

      rft = filtre_lp->step(rft);
      rft -= filtre_hp->step(rft);
      ra->step(rft, y);
  #   if 0
      sinon
      {
      // Méthode telle que décrite dans Wikipédia
      // https://fr.wikipedia.org/wiki/Modulation_d%27amplitude#D%C3%A9modulation_par_DSP

      ArrayXf ol = impl->ol.step(n).real();
      ArrayXf rft = rf * ol;

      Figure f;

      f.subplot(221);
      f.plot(rft, "b-", "Transpo");

      // Filtrage passe-bande
      rft = impl->filtre_lp.step(rft);

      f.subplot(222);
      f.plot(rft, "b-", "Passe-bas");

      msg("RFT out : {} elems, min = {}, max = {}.", rft.rows(), rft.minCoeff(), rft.maxCoeff());

      rft -= impl->filtre_hp.step(rft);

      f.subplot(223);
      f.plot(rft, "b-", "Passe-haut");

      // (3) Sous-échantillonnage
      audio = impl->ra.step(rft);

      f.subplot(224);
      f.plot(audio, "b-", "Sous-ech");

      stdo << f;
      }
      //f.enregistrer("./build/ex/analogique/entier.jpg");
  #   endif
    }

    // BLU
    sinon
    {
      // Filtre de Hilbert :
      // sin(t) -> -cos(t)
      // cos(t) -> sin(t)
      // C'est à dire en fréquentiel :
      // cos(wt) = 0.5 * (exp(iwt) + exp(-iwt))
      // sin(wt) = 0.5/i * (exp(iwt) - exp(-iwt)) = 0.5 * i * (-exp(iwt) + exp(-iwt))
      // exp(iwt) --> -i * signe(w) * exp(iwt)

      // Imaginons qu'on aie un signal cos(wt) = 0.5 * (exp(iwt) + exp(-iwt))
      // --> f de hilbert --> - i * 0.5 * (exp(iwt) - exp(-iwt))
      //                  --> sin(wt)

      // TRF de hilbert : T(x) = x + i h(x)
      // T(exp(iwt)) = exp(iwt) + signe(w) * exp(iwt)
      //  --> suppression des fréquences négatives

      // On veut :
      //  - (1) Supprimer la partie inférieure du spectre (ou supérieur pour LSB)
      //  - (2) Puis dupliquer la partie supérieure

      // (1) TRF de hilbert :
      //   y1 = x + i h(x)
      // (2) Duplication :
      //   y2 = reel(y1)

      // => y2 = reel(x) - imag(h(x))
      //       = reel(x) - h(imag(x))

      soit ix = hilbert->step(imag(rftc));
      soit rx = retard->step(real(rftc));
      y1 = rx - ix;





      si(config.debug_actif)
      {
        Figures f;

        soit s = f.subplot();
        s.plot_psd(x, config.fe_rf);
        s.titres("PSD signal entree");
        s = f.subplot();
        s.plot_psd(rftc, config.fe_rf);
        s.titres("PSD apres ol");
        s = f.subplot(223);
        s.plot_psd(y1, config.fe_rf);
        s.titres("PSD signal reconstruit");
        f.afficher();
      }

      // Première idée naïve :
      // Laisse passer uniquement les fréquences positives
      // impl->hilbert.step(rftc); // x2 = retard(x) + i * h(y)
      // Attention, ce n'est pas un filtre qui travaille indiferrement sur Real / Imag
      // (ce qui est le cas pour les filtres habituelles, e.g. leur réponse est symétrique).
      // audio = rftc.real();
      //       = retard(x) - h(imag(rftc)) équivalent mais plus coûteux.
    }

    // Post - traitements audio
    // Filtrage passe-bande
    soit y_lp = filtre_lp->step(y1);

    soit y_hp = filtre_hp->step(y1);

    assertion(y_lp.rows() == y_hp.rows());

    soit y_bp = y_lp - y_hp;


    si(config.debug_actif)
    {
      Figures f;

      f.subplot().plot_psd(y1, config.fe_rf, "", "PSD signal reconstruit");
      f.subplot().plot_psd(y_lp, config.fe_rf, "", "PSD LP");
      f.subplot().plot_psd(y_hp, config.fe_rf, "", "PSD HP");
      f.subplot().plot_psd(y_bp, config.fe_rf, "", "PSD BP");
      f.afficher("dbg-am-demod-filtrage.png");
    }


    y = y_bp;


    /*
    // (1) Transposition
    ArrayXXcf rft;
    impl->osc.step(rf, rft);

    // (2) Filtage passe-bas
    // A FAIRE

    // (3) Sous-échantillonnage
    ArrayXXcf y;
    impl->ra.step(rft, y);

    // (4) On prends la partie réelle ?
    audio = y.real();
    */
  }
  //void demodulation(const ArrayXf &rf, ArrayXf &audio);
};

sptr<Filtre<float, float, AMConfig>> modulateurAM()
{
  retourne std::make_shared<ModulateurAM>();
}
sptr<Filtre<cfloat, float, AMConfig>> demodulateurAM()
{
  retourne std::make_shared<DemodulateurAM>();
}













struct FMDemod: Filtre<cfloat, cfloat, FMDemodConfig>
{
  // 1er filtre  : wideband
  // 2eme filtre : après discri
  // 3eme filtre : audio
  sptr<FiltreGen<cfloat>> rif_wb;
  sptr<FiltreGen<float>>  rif_fm, rif_audio_L, rif_audio_R;
  float fech2 = 0;
  // Ratio de décimation
  entier R = 1;

  sptr<FiltreGen<cfloat,float>> discri = discriminateur_fm();

  sptr<Filtre<float, float, RPLLConfig>> pll = rpll_création(RPLLConfig());


  /*sptr<Filtre<cfloat, float, DemodConfig>>*/sptr<Démodulateur> demod_rds;

  Tabi H;
  Veci offsets[4];
  Veci rdec;

  void configure_impl(const FMDemodConfig &cfg)
  {
    soit &config = Configurable<FMDemodConfig>::config;
    config = cfg;


    rdec.setZero(26 * 4 * 2, 1);

    //H.resize(26, 10);


    H = vconcat(Veci::zeros(10*10), Veci::valeurs(
        {
              1, 0, 1, 1, 0, 1, 1, 1, 0, 0,
              0, 1, 0, 1, 1, 0, 1, 1, 1, 0,
              0, 0, 1, 0, 1, 1, 0, 1, 1, 1,
              1, 0, 1, 0, 0, 0, 0, 1, 1, 1,
              1, 1, 1, 0, 0, 1, 1, 1, 1, 1,
              1, 1, 0, 0, 0, 1, 0, 0, 1, 1,
              1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
              1, 1, 0, 1, 1, 1, 0, 1, 1, 0,
              0, 1, 1, 0, 1, 1, 1, 0, 1, 1,
              1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
              1, 1, 1, 1, 0, 1, 1, 1, 0, 0,
              0, 1, 1, 1, 1, 0, 1, 1, 1, 0,
              0, 0, 1, 1, 1, 1, 0, 1, 1, 1,
              1, 0, 1, 0, 1, 0, 0, 1, 1, 1,
              1, 1, 1, 0, 0, 0, 1, 1, 1, 1,
              1, 1, 0, 0, 0, 1, 1, 0, 1, 1
        }
        )).reshape(26,10);


    //H.topRows(10) = Eigen::MatrixXi::Identity();
    /*H.bottomRows(16)
     << 1, 0, 1, 1, 0, 1, 1, 1, 0, 0,
        0, 1, 0, 1, 1, 0, 1, 1, 1, 0,
        0, 0, 1, 0, 1, 1, 0, 1, 1, 1,
        1, 0, 1, 0, 0, 0, 0, 1, 1, 1,
        1, 1, 1, 0, 0, 1, 1, 1, 1, 1,
        1, 1, 0, 0, 0, 1, 0, 0, 1, 1,
        1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
        1, 1, 0, 1, 1, 1, 0, 1, 1, 0,
        0, 1, 1, 0, 1, 1, 1, 0, 1, 1,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        1, 1, 1, 1, 0, 1, 1, 1, 0, 0,
        0, 1, 1, 1, 1, 0, 1, 1, 1, 0,
        0, 0, 1, 1, 1, 1, 0, 1, 1, 1,
        1, 0, 1, 0, 1, 0, 0, 1, 1, 1,
        1, 1, 1, 0, 0, 0, 1, 1, 1, 1,
        1, 1, 0, 0, 0, 1, 1, 0, 1, 1;*/

    //pour(soit k = 0; k < 4; k++)
    //  offsets[k].resize(10);
    offsets[0] = Veci::valeurs({1, 1, 1, 1, 0, 1, 1, 0, 0, 0});
    offsets[1] = Veci::valeurs({1, 1, 1, 1, 0, 1, 0, 1, 0, 0});
    offsets[2] = Veci::valeurs({1, 0, 0, 1, 0, 1, 1, 1, 0, 0});
    offsets[3] = Veci::valeurs({1, 1, 0, 1, 0, 1, 0, 0, 0, 0});


   // Filtrage +/- 100 kHz
    soit coefs = tsd::filtrage::design_rif_cs(255, 0.2, 100e3/config.fe);
    rif_wb = tsd::filtrage::filtre_rif<float, cfloat>(coefs);

    // Filtrage après discri
    soit coefs2 = tsd::filtrage::design_rif_cs(255, 0.1, 65e3/config.fe);
    rif_fm = tsd::filtrage::filtre_rif<float, float>(coefs2);

    // Ré-échantillonnage
    fech2 = config.fe;
    R = 1;
    si(config.fe >= 200)
    {
      R = (entier) floor(config.fe / (65e3 * 2));
      fech2 = config.fe / R;
    }
    msg("Configuration démod FM : fe = {:.1f} Hz, fe2 = {:.1f} Hz (R = 1/{}).", config.fe, fech2, R);

    ModConfig mc;
    DemodConfig dc;
    dc.debug_actif = config.genere_img_debug;
    mc.fe     = fech2;
    mc.fi     = 57e3;
    mc.fsymb  = 1187.5 * 2; // Encodage de Manchester
    mc.forme_onde     = forme_onde_bpsk();
    demod_rds = démodulateur_création(mc, dc);

    // Filtrage audio
    soit coefs3 = tsd::filtrage::design_rif_cs(255, 0.1, 15e3/fech2);
    rif_audio_L = tsd::filtrage::filtre_rif<float, float>(coefs3);
    rif_audio_R = tsd::filtrage::filtre_rif<float, float>(coefs3);

    //discri->configure({config.fe, config.fe});

    RPLLConfig pll_config;
    pll_config.pll_interne.freq = 19.0e3f/fech2;
    pll_config.pll_interne.bp   = 250.0f/fech2;
    pll->configure(pll_config);

    //fech_sortie = fech2;

    si(config.genere_img_debug)
    {
      Figures f;
      f.subplot().plot(coefs, "b-", "Filtre wb");
      f.subplot().plot(coefs2, "b-", "Filtre fm");
      f.subplot().plot(coefs3, "b-", "Filtre audio");
      f.afficher("fm-demod-coefs.png");
    }
  }

  void step(const Veccf &x, Veccf &y)
  {
    soit &config = Configurable<FMDemodConfig>::config;
    soit y1 = rif_wb->step(x);
    soit a = discri->step(y1);
    a /= 75e3; // Excursion de 75 kHz

    soit a2 = rif_fm->step(a);
    soit a3 = sousech(a2, R);

    //ArrayXf x_rds = demod_rds->step(ArrayXcf(a3));
    BitStream bs;
    Tabf llr;
    demod_rds->step(a3, bs, llr);
    // Manchester : 2 position à essayer à chaque fois
    // Donc on travaille sur des blocs de 26 * 4 * 2 bits à chaque fois


    soit N = 26 * 4 * 2;
    pour(auto i = 0; i < bs.lon(); i++)
    {
      // Mise à jour de la fenêtre
      rdec.head(N - 1) = rdec.tail(N -1).eval();
      rdec(N-1) = bs[i];//x_rds(i) > 0 ? 1 : 0;

      // Décodage manchester
      Tabi dec(1, 26 * 4);
      pour(auto j = 0; j < N / 2; j++)
        dec(j) = rdec(2*j) == rdec(2*j+1) ? 0 : 1;

      soit nerrs = 0u;
      pour(auto j = 0; j < 4; j++)
      {
        Tabi syndrome = dec.block(0, j * 26, 1, 26) * H;
        //mod2(syndrome);
        pour(auto c = 0; c < 26; c++)
          syndrome(0, c) = syndrome(0, c) % 2;
        //nerrs += (syndrome - offsets[j]);//.abs().sum(); // TODO
      }
      si(nerrs < 5)
      {
        msg("RDS nerrs = {}", nerrs);
        uint16_t blk[4] = {0};

        string station;

        //pour(auto k = 0u; k < 4; k++)
        //  blk[k] = mati_vers_binaire(dec.block(0, j * 26, 1, 16));
        si((blk[1] & 0x0f) == 0)
        {
          char st[5];
          st[0] = blk[2] & 0xff;
          st[1] = (blk[2] >> 8) & 0xff;
          st[2] = blk[3] & 0xff;
          st[3] = (blk[3] >> 8) & 0xff;
          st[4] = 0;
          station = string(st);
        }

        msg("station = [{}]", station);
      }
    }



    // On a donc L+R à 0 Hz - 15 kHz,
    // un signal CW de référence à 19 kHz
    // et L-R à 38 kHz [23 - 53 kHz]
    // Et le signal RDS à 57 kHz
    soit x_LpR = rif_audio_L->step(a3);
    //audio::wav_enregistre("./build/out-mono.wav", fech2, x_LpR);


    // Pb ici : il y a un décalage de phase entre x_LmR et x_LpR
    // --> Trouver la phase de la porteuse grâce au signal à 38k
    // Il vaut mieux le faire avec une PLL, plus général
    soit cw = pll->step(a3);

    //ArrayXcf osc = OLH_signal(a3.rows(), 38e3 / fech);
    soit x_LmR_non_filtre = square(cw) * a3;
    soit x_LmR = rif_audio_R->step(x_LmR_non_filtre);

    soit L = x_LpR + x_LmR;
    soit R = x_LpR - x_LmR;

    y.resize(L.rows());
    y.set_real(L);
    y.set_imag(R);

    //audio::wav_enregistre_stereo("./build/out-stereo.wav", fech2, z);






    si(config.genere_img_debug)
    {

      {
        Figures f;
        f.subplot().plot_psd(y, 1e6, "", "PSD avant filtrage");
        f.subplot().plot_psd(y1, 1e6, "", "PSD apres filtrage");
        f.afficher("2-fig-psd");
      }
      {
        Figures f;
        f.subplot().plot_psd(a, 1e6, "", "PSD avant filtrage");
        f.subplot().plot_psd(a2, 1e6, "", "PSD apres filtrage");
        f.subplot().plot_psd(a3, fech2, "", "PSD apres sous-ech");
        f.enregistrer("3-fig-psd-demod");
      }
      {
        Figures f;
        f.subplot().plot_psd(a3, fech2, "", "PSD LpR avant filtrage");
        f.subplot().plot_psd(x_LpR, fech2, "", "PSD LpR apres filtrage");
        f.enregistrer("4-fig-mono");
      }
      {
        Figures f;
        f.subplot(421).plot(cw, "b-", "Signal CW regenere");
        f.subplot(422).plot_psd(cw, fech2, "", "PSD CW");
        f.subplot(423).plot_psd(a3, fech2, "", "Signal complet");
        f.subplot(424).plot_psd(x_LmR_non_filtre, fech2, "", "PSD LmR non filtre");
        f.subplot(425).plot_psd(x_LmR, fech2, "", "PSD LmR filtre");
        f.afficher("5-fig-pll");
      }
      {
        Figures f;
        f.subplot().plot(L, "b-", "Signal L");
        f.subplot().plot(R, "b-", "Signal R");
        f.afficher("6-fig-stereo");
      }
    }
  }
};

sptr<Filtre<cfloat, cfloat, FMDemodConfig>> demodulateurFM()
{
  retourne std::make_shared<FMDemod>();
}


}







