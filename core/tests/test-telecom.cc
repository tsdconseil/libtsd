#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"
#include "tsd/figure.hpp"
#include "tsd/fourier.hpp"
#include "tsd/tests.hpp"
#include <cstdio>

using namespace tsd;
using namespace tsd::telecom;
using namespace tsd::vue;
using namespace tsd::filtrage;


static void verifie_delais(const ArrayXf &x, int pos_attendue, const std::string &desc)
{
  //tsd_assert_msg((x.minCoeff() == 0) && (x.maxCoeff() == 1) && (x.sum() == 1), "vérification délais : valeur invalide (v = {})", x.transpose());
  int pos;
  x.maxCoeff(&pos);
  msg("Vérification délais {} : pos attendue = {}, pos = {}", desc, pos_attendue, pos);
  tsd_assert_msg(pos == pos_attendue, "vérification délais {} : erreur (pos attendue = {}, pos = {})", desc, pos_attendue, pos);
}


int test_bitstream()
{
  BitStream bs = BitStream::zéros(5);
  tsd_assert(bs.lon() == 5);
  for(auto i = 0; i < 5; i++)
  {
    tsd_assert(!bs[i]);
  }

  bs = BitStream::uns(13);
  tsd_assert(bs.lon() == 13);
  for(auto i = 0; i < 13; i++)
  {
    tsd_assert(bs[i]);
  }

  bs.set(4, false);
  tsd_assert(!bs[4]);

  bs.set(4, true);
  BitStream bs2 = bs + BitStream::zéros(7);
  tsd_assert(bs2.lon() == 20);
  for(auto i = 0; i < 13; i++)
  {
    tsd_assert(bs2[i]);
  }
  for(auto i = 14; i < 20; i++)
  {
    tsd_assert(!bs2[i]);
  }


  return 0;
}



void test_forme_onde(sptr<FormeOnde> fo)
{
  msg("Test forme d'onde [{}]...", fo->desc());


  BitStream bs("0000111101010101"), bs2;

  while((bs.lon() % fo->infos.k) != 0)
    bs.push(0);

  ArrayXcf x = fo->génère_symboles(bs);
  fo->decode_symboles(bs2, x);

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(bs.array(), "|b", "Bitstream");
    f.subplot().plot(x.real(), "ob", "Symboles (re)");
    f.gcf().plot(x.imag(), "og", "Symboles (im)");
    f.subplot().plot(bs2.array(), "|b", "Bitstream décodé");
    f.afficher();
  }

  tsd_assert_msg(bs.dst_Hamming(bs2) == 0, "Echec test forme d'onde.");
}


int test_formes_ondes()
{
  msg_majeur("Test des formes d'ondes");
  auto lst = {forme_onde_bpsk(), forme_onde_qpsk(), forme_onde_π4_qpsk(), forme_onde_psk(8), forme_onde_qam(16)};

  for(auto fo: lst)
    test_forme_onde(fo);

  return 0;
}





int test_delais_filtres()
{
  msg_majeur("Test délais des filtres.");


  int nc = 7;
  ArrayXf h = design_rif_fen(7, "lp", 0.3, "hn");
  ArrayXf x = ArrayXf::Zero(15);
  x(0) = 1;
  auto f11 = tsd::filtrage::filtre_rif<float,float>(h);
  ArrayXf y = f11->step(x);

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(h, "|bo", "h");
    f.subplot().plot(x, "|go", "x");
    f.subplot().plot(y, "|mo", "y1 [filtrage simple]");
    f.afficher();
  }

  verifie_delais(x,  0, "x");
  verifie_delais(y, rif_delais(nc),   "filtrage simple");

  for(auto R: {1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 31, 32})
  {
    auto filtre = tsd::filtrage::filtre_rif_ups<float,float>(h, R);
    ArrayXf y = filtre->step(x);
    if(tests_debug_actif)
    {
      Figure f;
      f.plot(y, "|mo", "upsampling R = {}", R);
      f.afficher();
    }
    verifie_delais(y,  filtre_rif_ups_délais(nc, R), fmt::format("upsampling R = {}", R));
  }

  return 0;
}


void test_discri_fm()
{
  // Signal test
  ArrayXf x = siggauss(1000) / 10;

  // Modulation FM
  ArrayXcf y = polar(cumsum(x));

  // Démod
  auto discri = discriminateur_fm();
  ArrayXf z = discri->step(y);

  float score;
  ArrayXcf x2 = x;
  ArrayXcf z2 = z;
  auto d = tsd::fourier::estimation_delais_entier(x2, z2, score);

  msg("Délais détecté : {} (corr = {})", d, score);

  //int n = z.rows();
  //z.head(n-2) = z.tail(n-2).eval();

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x, "", "Signal test");
    f.subplot().plot(y, "", "Signal modulé FM");
    f.subplot().plot(z, "", "Sortie discriminateur");
    f.subplot().plot(z - x, "-r", "Erreur");
    f.afficher();
  }

  ArrayXf err = z - x;

  float emax = err.abs().maxCoeff();
  msg("Erreur discri = {}", emax);
  if(emax > 1e-3)
    echec("Disci fm : trop d'erreur.");

}


// Attention, ne marque qu'en NRZ pour l'instant !!!!
//ArrayXcf tr_module(const BitStream &bs, sptr<FormeOnde> fo, int osf)
//{
  //ArrayXcf x = 2 * (bs.vers_array() - 0.5);
  //ArrayXcf x = fo->gene_symboles(bs);
  //assert(x.rows() == bs.lon());
  //msg("tr module : ")
  //if(fo->filtre.)
  //return sah(x, osf);
//}


struct TestRecepteurConfig
{
  int osf = 1;
  sptr<FormeOnde> fo = forme_onde_bpsk();
  float SNR_min = 0;
  float SNR_max = 10;
  bool  avec_delais = false;
  bool  avec_plot = false;
  int   nb_rep = 9;
  int   ncoefs      = 63;//15;
  bool  check_errs  = true;
  bool  check_trames = true;
  bool  premiere_trame_SNR_inf = true;
  int   nbits = -1;
  float carrier_rec_bl = 0.001;
  float seuil = 0.6;
  int   ncoefs_filtre_mise_en_forme = 63;

  struct
  {
    bool actif = false;
  } cocanal;


  int nreg_mls = -1;
};


struct TestRecepteurRes
{
  ArrayXf SNR;
  ArrayXf ber; // NaN si trame non reçue
  ArrayXf ber_theo;
  ArrayXf EbN0;
  ArrayXf err_phase, err_pos;

  bool succès = true;
};

TestRecepteurRes test_recepteur_unit(const TestRecepteurConfig &config)
{
  TestRecepteurRes res;
  srand(0x124DF531);

  if(config.avec_plot)
  {
    stdo.flush();
    stdo.def_dossier_sortie(fmt::format("./build/test-log/recepteur-osf-{}-{}-{}", config.osf, *(config.fo), config.avec_delais ? "delais-frac" : "delais-entier"));
  }

  msg_majeur("\n\n\n\n\nTest récepteur, OSF = {}, {}, filtre = {}, {}...\n\n\n\n\n", config.osf, *(config.fo), config.fo->filtre, config.avec_delais ? "delais-frac" : "delais-entier");
  RécepteurConfig rc;

  int nreg_mls = 7 + config.fo->infos.k;

  if(config.nreg_mls >= 0)
    nreg_mls = config.nreg_mls;
  rc.format.entete  = code_mls(nreg_mls); // 127 bits

  msg("Dimension de l'en-tête : {} bits.", rc.format.entete.lon());

  // Note : le padding de l'en-tête (multiple de k)
  // est fait automatiquement dans l'émetteur et dans le récepteur

  ArrayXf full_corr;

  int nbits = config.nbits == -1 ? 64 : config.nbits;

  auto fo = config.fo;

  rc.BS                                             = 512;
  rc.format.nbits                                   = nbits;
  //rc.format.fo_entete                               = fo;
  rc.seuil                                          = config.seuil;
  rc.debug_actif                                    = config.avec_plot;
  rc.config_demod.dec.carrier_rec.actif             = false;
  rc.config_demod.dec.clock_rec.actif               = false;
  rc.format.modulation.fsymb                        = 1e3;
  rc.format.modulation.fe                           = rc.format.modulation.fsymb * config.osf;
  rc.format.modulation.fi                           = 0;
  rc.format.modulation.wf                           = fo;
  rc.config_demod.debug_actif                       = config.avec_plot;
  rc.config_demod.dec.carrier_rec.BL                = config.carrier_rec_bl;
  rc.SNR_mini                                       = -10;
  rc.format.modulation.ncoefs_filtre_mise_en_forme  = config.ncoefs_filtre_mise_en_forme;
  if(fo->filtre.type == SpecFiltreMiseEnForme::NRZ)
    rc.format.modulation.ncoefs_filtre_mise_en_forme = config.osf;
  //rc.ncoefs_interpolateur                           = 63;
  rc.callback_corr = [&](const ArrayXf &corr)
  {
    full_corr = vconcat(full_corr, corr);
  };
  auto rec = récepteur_création(rc);

  ÉmetteurConfig ec;
  ec.format       = rc.format;
  ec.debug_actif  = config.avec_plot;

  auto émet = émetteur_création(ec);

  BitStream data;

  if(config.nbits == -1)
  {
    BitStream data1 = BitStream("1010101011110000");
    data = data1;
    for(auto i = 0; i < 3; i++)
      data += data1;
  }
  else
  {
    data = BitStream::rand(config.nbits);
  }

# if 0
  float retard;
  BitStream bs = rc.format.entete + data;
  ArrayXcf xd = config.fo->génère_échantillons(bs, config.ncoefs, config.osf, retard);
# endif

  float retard = émet->retard();

  ArrayXcf xd = émet->step(data);

  int nd = xd.rows();

  retard -= config.osf / 2.0f; // Début du symbole

  msg(" *** retard modulation = {}", retard);

  int Nt = 24 * 1024;

  /*if(Nt < ((nd + 1024) * config.nb_rep + 1024))
  {
    Nt = (nd + 1024) * config.nb_rep + 1024;
  }*/

  // 1023 * 4 / 3 = 1300
  auto re = rec->get_etat();
  int nb_echans_padding = re.Ne;//(nbits_entete * config.osf) / config.fo->k;

  msg(" Nb échan padding = {}", nb_echans_padding);


  int dim_ideale = (nd + 1024) * config.nb_rep + 2 * nb_echans_padding;

  // Ajout de zéros de manière à flusher le détecteur
  Nt = std::max(Nt, dim_ideale);


  // (8192 est la dimension des blocs d'entrée)
  //if(Nt < ((nd + 1024) * config.nb_rep + 8192 + 2048))
  //{
    //Nt = (nd + 1024) * config.nb_rep + 8192 + 2048;
  //}

  ArrayXcf x = 0.0001 * randn(Nt);

  int nb_rep = config.nb_rep;//9;

  ArrayXf SNR = linspace(config.SNR_min, std::max(config.SNR_max, config.SNR_min), nb_rep).reverse();

  if(config.premiere_trame_SNR_inf)
    SNR(0) = 200; // Test de référence, (presque) pas de bruit sur la première trame
  std::vector<float> pos(nb_rep);

  res.SNR = SNR;
  res.ber = res.ber_theo = res.EbN0 = ArrayXf::Constant(nb_rep, 1);//;//-1);
  res.err_phase = res.err_pos = ArrayXf::Constant(nb_rep, std::nan(""));

  int offset = 100;
  for(auto i = 0; i < nb_rep; i++)
  {


    tsd_assert(offset + nd <= x.rows());

    //float σ = pow(10, - SNR(i) / 20) / sqrt(2.0f);
    // Pareil
    float σ = sqrt((1/db2pow(SNR(i)))) / sqrt(2.0f);



    float τ = 0;
    if(config.avec_delais)
      τ = 0.5;//randu(1)(0);

    pos[i] = offset + retard + τ;


    msg("Placement trame[{}] @ offset {}, SNR={:.1f} dB, σ={}.", i, pos[i], SNR(i), σ);

    ArrayXcf xde;

    if(config.avec_delais)
    {
      xde = tsd::fourier::delais(xd, τ);

      if(config.avec_plot)
      {
        if(tests_debug_actif)
        {
          Figures f;
          f.subplot().plot(xd, "", "Avant retard");
          f.subplot().plot(xde, "", "Après retard τ={}", τ);
          f.afficher("Insertion : avant / après retard");
        }
      }
    }
    else
      xde = xd;

    float puissance_trame = xde.abs2().mean();

    σ *= sqrt(puissance_trame);

    float SNR_emp = 10 * log10(xde.abs2().mean() / (2*σ*σ));
    msg("Energie trame : {}, σ² : {}, SNR empirique = {:.1f} dB.", xde.abs2().mean(), 2*σ*σ, SNR_emp);


    ArrayXcf xr(nd);
    xr.real() = σ * randn(nd);
    xr.imag() = σ * randn(nd);

    float σ_emp     = sqrt(xr.abs2().mean());
    float σ_emp_re  = sqrt(xr.real().abs2().mean());
    float σ_emp_im  = sqrt(xr.imag().abs2().mean());

    msg("σ cible = {} (/sqrt(2)={}), σ emp = {}, σ emp re = {}, σ emp im = {}",
        σ * sqrt(2.0f), σ, σ_emp, σ_emp_re, σ_emp_im);

    ArrayXcf co = ArrayXcf::Zero(nd);

    if(config.cocanal.actif)
    {
      // Tests Pascal :
      //  - QPSK ou pi/4 QPSK : 156.25 ks/s (312 kbit/s)
      //  - ASK 40 KSPS, PRBS9
      //     -> soit 40/156.25 de fsymb
      //  - Décalage ~ 8 kHz
      //     -> soit 8/156.25 de fsymb

      auto fo = forme_onde_ask(2);

      ModConfig cfg_ask = rc.format.modulation;
      cfg_ask.fi    = rc.format.modulation.fsymb * 8.0f / 156.25f;
      cfg_ask.fsymb = rc.format.modulation.fsymb * 40.0f / 156.25f;
      cfg_ask.sortie_reelle = false;

      auto mod = modulateur_création(cfg_ask);

      int nbits = nd * rc.format.modulation.fe / (rc.format.modulation.fsymb / 10) + 128;
      BitStream bs = BitStream::rand(nbits);
      ArrayXcf y = mod->step(bs);
      tsd_assert(y.rows() >= nd);

      co = y.head(nd) / 2;

      if(config.avec_plot)
      {
        {
          Figure f;
          f.plot(co);
          f.titre("Modulation ASK parasite");
          f.afficher();
        }
        {
          Figure f;
          f.plot_iq(co, "ob");
          f.titre("Modulation ASK parasite (constel)");
          f.afficher();
        }
      }
    }


    x.segment(offset, nd) = xde + xr + co;


    if(config.avec_plot)
    {
      ArrayXcf xe = x.segment(offset, nd);
      Figures f;
      f.subplot().plot(xde, "", "Avant bruit");
      f.subplot().plot(xe, "", "Après bruit σ={:.3f}, SNR={:.1f} dB", σ, SNR(i));
      f.afficher("Insertion : avant / après retard");

      {
        Figure f;
        f.plot_iq(xde, "ob", "Avant bruit");
        f.afficher("Constellation, avant bruit");
      }
      {
        Figure f;
        f.plot_iq(xe, "ob", "Après bruit");
        f.afficher("Constellation, avant bruit");
      }
    }

    if(std::abs(SNR_emp - SNR(i)) >= 0.5)
    {
      echec("SNR empirique hors borne.");
    }

    offset += nd;
    int ecart = (rand() & 511) + 128 * 4;
    x.segment(offset, ecart) = bruit_awgn(ArrayXcf::Zero(ecart), σ);
    offset += ecart;
  }

  //x = x.head(offset + 8 * 1024).eval(); // laisse quelques échantillons à la fin pour flusher

  if(config.avec_plot)
  {
    Figures f;
    f.subplot().plot(rc.format.entete.array(), "", "en-tête");
    f.subplot().plot(data.array(), "", "data");
    f.subplot().plot(x, "", "x");
    f.afficher(fmt::format("Test récepteur, osf={}, mod={}", config.osf, *(config.fo)));
  }

  auto trames = rec->step(x);

  msg("Nb trames décodées : {}", trames.size());
  //tsd_assert(!trames.empty());


  if(config.avec_plot)
  {
    Figures fg;
    fg.subplot().plot(x, "", "Vue globale");

    for(auto j = 0; ((j < (int) trames.size()) && (j < nb_rep)); j++)
    {
      auto &t = trames[j];
      fg.gcf().plot((float) t.det.position, 1, "sr");
    }

    fg.subplot().plot(full_corr, "", "Corrélation");
    fg.afficher();
  }

  for(auto j = 0; j < nb_rep; j++)
  {
    float EsN0 = pow2db(config.osf * db2pow(SNR(j)));
    float EbN0 = EsN0 - 10 * log10(fo->infos.k);
    float ber_theo = fo->ber(EbN0);
    res.ber_theo(j) = ber_theo;
    res.EbN0(j)     = EbN0;
  }

  for(auto j = 0; ((j < (int) trames.size()) && (j < nb_rep)); j++)
  {
    msg("Vérification trames #{}...", j);
    auto &t = trames[j];
    msg("  détection : {}", t.det);
    //msg("  données décodées : {}", t.bs);


    if(config.avec_plot)
    {
      Figures f;
      f.subplot().plot(t.bs.array(), "|b", "Données reçues");
      f.subplot().plot(data.array(), "|g", "Données envoyées");
      f.afficher(fmt::format("Trame décodée vs envoyée (#{})", j));
    }
    if(data.lon() != t.bs.lon())
    {
      msg_avert("Nombre de bits reçu invalide ({} vs {} attendus)", t.bs.lon(), data.lon());
      res.succès = false;
      return res;
    }
    auto dst = data.dst_Hamming(t.bs);

    float ber = ((float) dst) / t.bs.lon();

    res.ber(j) = ber;

    float EsN0 = pow2db(config.osf * db2pow(SNR(j)));
    float EbN0 = EsN0 - 10 * log10(fo->infos.k);
    float ber_theo = fo->ber(EbN0);
    msg("  SNR = {:.1f} dB, EsN0 = {:.1f} dB, EbN0 = {:.1f} dB, err=\033[32m{}\033[0m, ber=\033[32m{:.2e}\033[0m, ber theo={:.2e}", SNR(j), EsN0, EbN0, dst, ber, ber_theo);
    auto ep = t.det.position_prec - pos[j];
    msg("  Erreur position = {}.", ep);


    res.err_pos(j) = ep;
    res.err_phase(j) = t.det.θ;

    res.EbN0(j)     = EbN0;

    if(config.check_trames)
    {
      if(std::abs(ep) >= 0.2)
      {
        msg_avert("Position invalide (attendue : {}, détectée : {}).", pos[j], t.det.position_prec);
        res.succès = false;
        return res;
      }
    }

    float nb_err_attendues = ber_theo * data.lon();

    // Tolérance car le BER théorique ne prends pas en compte
    // l'erreur de phase et de gain
    if(config.check_errs && (dst > 5 + nb_err_attendues * 3.0))
    {
      msg_avert(" BER trop élevé (nb err attendues = {:.2f}, nb err détectées = {}).", nb_err_attendues, dst);
      res.succès = false;
      return res;
    }
  }

  if(config.check_trames)
  {
    if((int) trames.size() != nb_rep)
    {
      msg_avert("nombre de trames reçues invalide ({} attendues, {} reçues).", nb_rep, trames.size());
      res.succès = false;
      return res;
    }
  }

  if(config.avec_plot)
    stdo.flush();

  res.EbN0 = res.EbN0.reverse().eval();
  res.SNR = res.SNR.reverse().eval();
  res.ber = res.ber.reverse().eval();
  res.ber_theo = res.ber_theo.reverse().eval();
  res.err_phase = res.err_phase.reverse().eval();
  res.err_pos = res.err_pos.reverse().eval();

  return res;
}


int bench_recepteur_a()
{
  auto wf = forme_onde_qpsk();
  wf->filtre = SpecFiltreMiseEnForme::srrc(0.25);


  /*{
    test_recepteur_unit(
    {
      .osf                    = 4,
      .fo                     = wf,
      .SNR_min                = -4,
      .SNR_max                = 8,
      .avec_delais            = false,
      .avec_plot              = true,
      .nb_rep                 = 13,
      .check_errs             = false,
      .check_trames           = false,
      .premiere_trame_SNR_inf = false,
      .nbits                  = 32,
      .carrier_rec_bl         = 0,
      .seuil                  = 0.4,
    });
  }
  return 0;*/



  auto res = test_recepteur_unit(
  {
    .osf                    = 4,
    .fo                     = wf,
    .SNR_min                = -4,
    .SNR_max                = 8,
    .avec_delais            = false,
    .avec_plot              = false,
    .nb_rep                 = 13,
    .check_errs             = false,
    .check_trames           = false,
    .premiere_trame_SNR_inf = false,
    .nbits                  = 16 * 1024,
    .carrier_rec_bl         = 0,
    .seuil                  = 0.4,
    .ncoefs_filtre_mise_en_forme  = 63
  });

  {
    Figure f;
    f.plot(res.SNR, res.EbN0, "b-o", "Eb/N0 = f(SNR)");
    f.afficher("BPSK - Eb/N0");
  }

  {
    Figure f;
    f.axes().def_echelle("lin", "log");
    f.plot(res.EbN0, res.ber_theo, "g-o", "Ber théorique");
    f.plot(res.EbN0, res.ber, "b-o", "Ber simulé");
    f.titres("BPSK - BER = f(Eb/N0)", "Eb/N0", "ber");
    f.def_rdi({res.EbN0.minCoeff() - 1, 1e-6, res.EbN0.maxCoeff() + 1, 1.0});
    f.afficher("BPSK - BER");
  }
  {
    Figures f;
    f.subplot().plot(res.SNR, res.err_phase * (180 / π), "r-o", "Erreur de phase (degrés)");
    f.gcf().titres("", "SNR (dB)");
    f.subplot().plot(res.SNR, res.err_pos, "r-o", "Erreur de position");
    f.gcf().titres("", "SNR (dB)");
    f.afficher("BPSK - Erreurs détection");
  }





  return 0;
}



int bench_recepteur()
{
  bench_recepteur_a();
  //return 0;
  auto lst_m =
  {
      forme_onde_bpsk(),
      forme_onde_qpsk(),
      forme_onde_π4_qpsk(),
      forme_onde_psk(8),
      forme_onde_qam(16),
      forme_onde_fsk(4, 1.0f), // Index de 1.0f en 4FSK, pour meilleure discrimination
  };
  //auto lst_f  = {SpecFiltreMiseEnForme::nrz(), SpecFiltreMiseEnForme::srrc(0.5), SpecFiltreMiseEnForme::srrc(0.2), SpecFiltreMiseEnForme::gaussien(0.8)};
  int osf     = 4;

  FILE *fo = fopen("./build/test-log/bench-recepteur.txt", "wt");

  bool first = true;

  for(auto m: lst_m)
  {
    m->filtre = m->infos.est_fsk ? SpecFiltreMiseEnForme::gaussien(2.0) : SpecFiltreMiseEnForme::srrc(0.25);
    auto res = test_recepteur_unit(
    {
      .osf                    = osf,
      .fo                     = m,
      .SNR_min                = -4,
      .SNR_max                = 16,
      .avec_delais            = false,
      .avec_plot              = false,
      .nb_rep                 = 11,
      .check_errs             = false,
      .check_trames           = false,
      .premiere_trame_SNR_inf = false,
      .nbits                  = 16 * 1024,
      .carrier_rec_bl         = 0,
      .ncoefs_filtre_mise_en_forme  = 63
    });

    if(first)
    {
      std::string s = "Modulation ; SNR";
      for(auto i = 0; i < res.SNR.rows(); i++)
        s += fmt::format(" ; {:.1f} dB", res.SNR(i));
      fprintf(fo, "%s\n", s.c_str());
      first = false;
    }

    std::string s = fmt::format("{}", *m);

    s += " ; Eb/N0";
    for(auto i = 0; i < res.ber.rows(); i++)
      s += fmt::format(" ; {:.1f} dB", res.EbN0(i));
    s += "\n";

    s += " ; Ber théo";
    for(auto i = 0; i < res.ber.rows(); i++)
      s += fmt::format(" ; {:.4f} %", 100 * res.ber_theo(i));
    s += "\n";

    s += " ; Ber simulé";
    for(auto i = 0; i < res.ber.rows(); i++)
    {
      if(res.ber(i) >= 0)
        s += fmt::format(" ; {:.4f} %", 100 * res.ber(i));
      else
        s += " ; n/d";
    }
    fprintf(fo, "%s\n", s.c_str());
  }

  fclose(fo);

  return 0;
}

// Code = MLS 31 bits
// BPSK
// Signal 1.25 MHz
// sous-échantilloné à 250 kHz après filtre adapté NRZ
int test_recepteur()
{
  auto lst_m =
  {
      forme_onde_bpsk(),
      forme_onde_qpsk(),
      forme_onde_π4_qpsk(),
      forme_onde_psk(8),
      forme_onde_fsk(),
      forme_onde_fsk(4, 1.0f), // Index de 1.0f en 4FSK, pour meilleure discrimination
      forme_onde_qam(4),
      forme_onde_qam(16)
  };



# if 0
  // Vérif amplitudes
  {
    TestRecepteurConfig cfg;
    cfg.osf = 4;
    cfg.fo = forme_onde_qpsk();
    //cfg.fo->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    cfg.SNR_min = 4;
    cfg.nb_rep = 2;
    cfg.avec_plot = true;
    test_recepteur_unit(cfg);
  }


  {
    TestRecepteurConfig cfg;
    cfg.osf = 4;
    cfg.fo = forme_onde_π4_qpsk();
    //cfg.fo->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    cfg.SNR_min = 4;
    cfg.nb_rep = 2;
    cfg.avec_plot = true;
    test_recepteur_unit(cfg);
  }
# endif

# if 0
  // Vérif cocanal
  {
    TestRecepteurConfig cfg;
    cfg.osf = 8;
    cfg.fo = forme_onde_qpsk();
    cfg.fo->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    cfg.SNR_min = 12;
    cfg.nb_rep = 2;
    cfg.avec_plot = true;
    cfg.cocanal.actif = true;
    cfg.nreg_mls = 5;
    cfg.nbits     = 512;
    cfg.check_trames = false;
    //cfg.
    test_recepteur_unit(cfg);


    cfg.fo = forme_onde_π4_qpsk();
    cfg.fo->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    test_recepteur_unit(cfg);
  }
  exit(0);
# endif

  // Ceci fonctionnne !!!!
  //tsd_assert(test_recepteur_unit({.osf = 3, .fo = waveform_bpsk(), .SNR_min = 4, .avec_delais = false, .avec_plot = true, .nb_rep = 1}).succès);
  //tsd_assert(test_recepteur_unit({.osf = 3, .fo = waveform_bpsk(), .SNR_min = 4, .avec_delais = true, .avec_plot = true, .nb_rep = 1}).succès);

  /*{
    auto m = waveform_bpsk();
    m->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    tsd_assert(test_recepteur_unit({.osf = 3, .fo = m, .SNR_min = 4, .avec_delais = false, .avec_plot = true, .nb_rep = 1}).succès);
  }*/
  /*{
    auto fo = waveform_fsk(4);
    fo->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    test_recepteur_unit({.osf = 4, .fo = fo, .SNR_min = 15, .avec_delais = false, .avec_plot = true});
  }*/

  /*
  test_recepteur_unit({.osf = 4, .fo = waveform_bpsk(), .SNR_min = 4, .avec_delais = false, .avec_plot = true, .nb_rep = 1});
  test_recepteur_unit({.osf = 4, .fo = waveform_fsk(),  .SNR_min = 10, .avec_delais = false, .avec_plot = true, .nb_rep = 1});
  auto fo = waveform_fsk();
  fo->filtre = SpecFiltreMiseEnForme::gaussien(0.8);
  test_recepteur_unit({.osf = 4, .fo =fo,  .SNR_min = 10, .avec_delais = false, .avec_plot = true, .nb_rep = 1});*/
  //test_recepteur_unit({.osf = 2, .fo = waveform_psk(8), .SNR_min = 4, .avec_delais = false, .avec_plot = true});

  //auto lst_f = {SpecFiltreMiseEnForme::nrz(), SpecFiltreMiseEnForme::srrc(0.5)};
  auto lst_f = {SpecFiltreMiseEnForme::srrc(0.5), SpecFiltreMiseEnForme::nrz()};

  for(auto m: lst_m)
  {
    for(auto f: lst_f)
    {
      m->filtre = f;
      for(auto osf: {4, 2, 3, 4, 5, 6})
        for(auto avec_delais: {false, true})
        {
          float SNR_min = 4;
          if(m->infos.est_fsk)
            SNR_min = 15;
          auto res = test_recepteur_unit({.osf = osf, .fo = m, .SNR_min = SNR_min, .avec_delais = avec_delais, .avec_plot = false});
          if(!res.succès)
          {
            test_recepteur_unit({.osf = osf, .fo = m, .SNR_min = SNR_min, .avec_delais = avec_delais, .avec_plot = true});

            msg_erreur("Le test unitaire du récepteur a échoué : osf={}, fo={}", osf, *m);

            return -1;
          }
        }
    }
  }

  stdo.def_dossier_sortie("./build/test-log/recepteur");
  return 0;
}


int test_filtre_adapte()
{

  for(auto i = 0; i < 2; i++)
  {
    auto sf = SpecFiltreMiseEnForme::srrc(0.5);
    int osf = 2;

    auto fmef = sf.filtre_mise_en_forme(15, osf);
    auto fa   = sf.filtre_adapte(15, osf);

    int n = 30;
    ArrayXcf x = ArrayXcf::Zero(n);
    if(i == 0)
    {
      for(auto i = 0; i + 1 < n; i += 2)
      {
        x(i)   = -1;
        x(i+1) = 1;
      }
    }
    else
      x(n/2) = 1;

    ArrayXcf y = fmef->step(x);
    ArrayXcf z = fa->step(y);


    auto itrp = tsd::filtrage::itrp_sinc<cfloat>({15, 1024, 0.45, "hn"});

    ArrayXf c0 = itrp->coefs(0);
    ArrayXf c1 = itrp->coefs(0.5);

    if(tests_debug_actif)
    {
      {
        Figure f;
        f.plot(c0);
        f.afficher("Coef itrp délais 0");
      }

      {
        Figure f;
        f.plot(c1);
        f.afficher("Coef itrp délais 0.5");
      }
    }

    auto fc0 = filtre_rif<float, cfloat>(c0);
    auto fc1 = filtre_rif<float, cfloat>(c1);

    ArrayXcf z0 = fc0->step(z);
    ArrayXcf z1 = fc1->step(z);

    if(tests_debug_actif)
    {
      Figures f;
      f.subplot().plot(x, "ob-", "x : bitstream");
      f.subplot().plot(y, "og-", "y : filtre de mise en forme");
      f.subplot().plot(z, "om-", "z : filtrage adapté");
      f.subplot().plot(z0, "or-", "z0 : interpolation délais = 0");
      f.subplot().plot(z1, "or-", "z1 : interpolation délais = 0.5");
      f.afficher();
    }
  }
  return 0;
}




// TODO !!!
void test_fsk()
{
  // 16 bits aléatoires
  auto bs = randstream(16);

  ModConfig cfg;
  cfg.wf            = forme_onde_fsk(2, 4);//2, index, filt, BT);
  cfg.debug_actif   = true;
  cfg.fe            = 10e3;
  cfg.fi            = 500;
  cfg.fsymb         = 0.1e3;
  cfg.sortie_reelle = false;

  // Fréquence d'échantillonnage  = 10 kHz
  // Fréquence intermédiaire      = 500 Hz
  // Fréquence symbole            = 0.1 kHz
  // -> Déviation = fsymb * h / 2 = 200 Hz
  auto mod = modulateur_création(cfg);

  ArrayXcf x = mod->step(bs);

  auto discri = discriminateur_fm();
  ArrayXf y = discri->step(x);

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(bs.array(), "hb", "Train binaire");
    f.subplot().plot(x, "", "Modulation FSK");
    f.subplot().plot(y.tail(y.rows() * 0.9), "", "Discri");
    f.afficher("Test FSK");
  }
}

int test_demod()
{
  msg_majeur("Tests démodulation...");

  auto filtres = {SpecFiltreMiseEnForme::nrz(), SpecFiltreMiseEnForme::srrc(0.4)};
  auto Ms      = {2, 4, 8};

  std::vector<sptr<FormeOnde>> wfs;

  // TODO : ne fonctionne plus !!!
  /*wfs.push_back(waveform_fsk(2, 0.5,  SpecFiltreMiseEnForme::nrz()));
  wfs.push_back(waveform_fsk(4, 0.5,  SpecFiltreMiseEnForme::nrz()));
  wfs.push_back(waveform_fsk(2, 0.5,  SpecFiltreMiseEnForme::gaussien(1.5)));
  wfs.push_back(waveform_fsk(2, 1,    SpecFiltreMiseEnForme::nrz()));
  wfs.push_back(waveform_fsk(2, 0.7,  SpecFiltreMiseEnForme::nrz()));
  wfs.push_back(waveform_fsk(2, 2,    SpecFiltreMiseEnForme::nrz()));*/

  for(auto filtre: filtres)
  {
    wfs.push_back(forme_onde_qam(16, filtre));
    for(auto M: Ms)
      wfs.push_back(forme_onde_psk(M, filtre));
  }

  for(auto wf: wfs)
  {
    msg_majeur("Test démodulation {}", *wf);
    stdo.printf(fmt::format("<h2>Test démodulation {}</h2>", *wf));

    ModConfig modcfg;
    DemodConfig cfg;
    modcfg.wf          = wf;//waveform_psk(M, filtre);
    modcfg.fe          = 100e3;
    modcfg.fi          = 0;
    modcfg.fsymb       = 20e3;
    modcfg.debug_actif = false;
    cfg.debug_actif    = false;
    auto demod         = démodulateur_création(modcfg, cfg);

    // Create a modulator to simulate RF signal
    auto mod = modulateur_création(modcfg);

    auto bs = randstream(1000);
    ArrayXcf x = mod->step(bs);

    x = bruit_awgn(x, 0.01);

    // Proceed to demodulation
    BitStream bs2;
    ArrayXXf llr;
    demod->step(x, bs2, llr);

    // Comment assurer la synchro d'horloge ?

    /*Figure f;
    f.subplot();
    f.plot(bs.vers_array(), "hb", "Train binaire");
    f.subplot();
    f.plot(x, "", "Données modulées");
    f.subplot();
    f.plot(bs2.vers_array(), "hb", "Train binaire décodé");
    f.afficher("bits entrées / sorties");*/

    // Ignore les 300 premiers bits (convergence des boucles)
    bs = BitStream(bs.array().tail(700));
    bs2 = BitStream(bs2.array().tail(700));

    // Alignment of the two bit vectors (and phase ambiguity resolution)
    auto res = cmp_bits_psk(bs, bs2, wf->infos.k);
    //auto res = cmp_bits(bs, bs2);

    if(tests_debug_actif)
    {
      Figures f;
      f.subplot().plot(res.b0, "hb", "Train binaire émis");
      f.subplot().plot(res.b1, "hg", "Train binaire décodé");
      f.subplot().plot(res.b1 - res.b0, "hr", "Erreurs (nb errs = {}, ber = {:.1e})", res.nerr, res.ber);
      f.afficher("Bits entrées / sorties");
    }



    // Verification
    if(res.nerr != 0)
      echec("Bits erronés non attendus !");
  }
  return 0;
}



void test_sah()
{
  ArrayXf x(3);
  x << 0, 1, 2;
  ArrayXf y = sah(x, 3);
  ArrayXf yref(9);
  yref << 0, 0, 0, 1, 1, 1, 2, 2, 2;
  auto err = (y - yref).abs().maxCoeff();
  tsd_assert_msg(err == 0, "Erreur SAH");
  msg("sah ok.");
}


static void test_emetteur()
{
  RécepteurConfig rc;
  rc.format.entete = BitStream::rand(127);
  rc.format.modulation.wf     = forme_onde_bpsk();
  rc.format.modulation.fe     = 1e6;
  rc.format.modulation.fsymb  = 1e5;
  rc.format.nbits = 256;
  rc.debug_actif = true;

  ÉmetteurConfig ec;
  ec.debug_actif = true;
  ec.format = rc.format;

  msg("Création émetteur...");
  auto eme = émetteur_création(ec);
  msg("Création récepteur...");
  auto rec = récepteur_création(rc);

  auto bs = BitStream::rand(256);

  msg("Emission...");
  ArrayXcf x = eme->step(bs);

  x = vconcat(ArrayXcf::Zero(1000), x);
  x = vconcat(x, ArrayXcf::Zero(40e3));

  msg("Réception...");
  auto tr = rec->step(x);

  msg("Reçu : {} trames.", tr.size());

  tsd_assert(tr.size() == 1u);
  tsd_assert(tr[0].bs == bs);
}

int test_telecom()
{
  test_emetteur();
  test_filtre_adapte();
  test_discri_fm();
  test_fsk();
  test_demod();
  test_sah();
  return 0;
}
