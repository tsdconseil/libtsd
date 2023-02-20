#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


static void test_filtre_boucle_ordre_1()
{
  msg_majeur("Test filtre de boucle d'ordre 1...");

  // Constante de temps
  soit τ = 5;
  soit flt = filtre_boucle_ordre_1(τ);

  // Vérifie la constante de temps
  soit n = 100;
  soit y = Vecf::zeros(n);
  pour(auto i = 1; i < n; i++)
    y(i) = flt->step(1 - y(i-1));

  Figure f;
  f.plot(y);
  f.afficher("Test filtre boucle ordre 1");

  msg("y(5) = {}", y(5));

  // Attendu : 63 % à τ
  tsd_assert(abs(y(5) - 0.632) < 1e-3);

  msg("ok.");
}


static void verifie_delais(const Vecf &x, entier pos_attendue, const std::string &desc)
{
  soit pos = x.index_max();
  msg("Vérification délais {} : pos attendue = {}, pos = {}", desc, pos_attendue, pos);
  tsd_assert_msg(pos == pos_attendue, "vérification délais {} : erreur (pos attendue = {}, pos = {})", desc, pos_attendue, pos);
}


entier test_bitstream()
{
  BitStream bs = BitStream::zéros(5);
  tsd_assert(bs.lon() == 5);
  pour(auto i = 0; i < 5; i++)
  {
    tsd_assert(!bs[i]);
  }

  bs = BitStream::uns(13);
  tsd_assert(bs.lon() == 13);
  pour(auto i = 0; i < 13; i++)
  {
    tsd_assert(bs[i]);
  }

  bs.set(4, non);
  tsd_assert(!bs[4]);

  bs.set(4, oui);
  BitStream bs2 = bs + BitStream::zéros(7);
  tsd_assert(bs2.lon() == 20);
  pour(auto i = 0; i < 13; i++)
  {
    tsd_assert(bs2[i]);
  }
  pour(auto i = 14; i < 20; i++)
  {
    tsd_assert(!bs2[i]);
  }


  retourne 0;
}



void test_forme_onde(sptr<FormeOnde> fo)
{
  msg("Test forme d'onde [{}]...", fo->desc());


  BitStream bs("0000111101010101"), bs2;

  tantque((bs.lon() % fo->infos.k) != 0)
    bs.push(0);

  soit x = fo->génère_symboles(bs);
  fo->decode_symboles(bs2, x);

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(bs.array(), "|b", "Bitstream");
    f.subplot().plot(real(x), "ob", "Symboles (re)");
    f.gcf().plot(imag(x), "og", "Symboles (im)");
    f.subplot().plot(bs2.array(), "|b", "Bitstream décodé");
    f.afficher();
  }

  tsd_assert_msg(bs.dst_Hamming(bs2) == 0, "Echec test forme d'onde.");
}


entier test_formes_ondes()
{
  msg_majeur("Test des formes d'ondes");
  soit lst = {forme_onde_bpsk(), forme_onde_qpsk(), forme_onde_π4_qpsk(), forme_onde_psk(8), forme_onde_qam(16)};

  pour(auto fo: lst)
    test_forme_onde(fo);

  retourne 0;
}





entier test_delais_filtres()
{
  msg_majeur("Test délais des filtres.");

  entier nc = 7;
  soit h = design_rif_fen(7, "lp", 0.3, "hn");
  soit x = Vecf::zeros(15);
  x(0) = 1;
  soit f11 = tsd::filtrage::filtre_rif<float,float>(h);
  soit y = f11->step(x);

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(h, "|bo", "h");
    f.subplot().plot(x, "|go", "x");
    f.subplot().plot(y, "|mo", "y1 [filtrage simple]");
    f.afficher();
  }

  verifie_delais(x,  0, "x");
  verifie_delais(y, rif_delais(nc),   "filtrage simple");

  pour(auto R: {1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 31, 32})
  {
    soit filtre = tsd::filtrage::filtre_rif_ups<float,float>(h, R);
    soit y = filtre->step(x);
    si(tests_debug_actif)
    {
      Figure f;
      f.plot(y, "|mo", "upsampling R = {}", R);
      f.afficher();
    }
    verifie_delais(y,  filtre_rif_ups_délais(nc, R), fmt::format("upsampling R = {}", R));
  }

  retourne 0;
}


void test_discri_fm()
{
  // Signal test
  soit x = siggauss(1000) / 10;

  // Modulation FM
  soit y = polar(cumsum(x));

  // Démod
  soit discri = discriminateur_fm();
  soit z = discri->step(y);


  soit x2 = x.clone();
  soit z2 = z.clone();
  soit [d, score] = tsd::fourier::estimation_délais(x2, z2);

  msg("Délais détecté : {} (corr = {})", d, score);

  //entier n = z.rows();
  //z.head(n-2) = z.tail(n-2).eval();

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(x, "", "Signal test");
    f.subplot().plot(y, "", "Signal modulé FM");
    f.subplot().plot(z, "", "Sortie discriminateur");
    f.subplot().plot(z - x, "-r", "Erreur");
    f.afficher();
  }

  soit err = z - x;

  soit emax = abs(err).valeur_max();
  msg("Erreur discri = {}", emax);
  si(emax > 1e-3)
    echec("Disci fm : trop d'erreur.");

}


// Attention, ne marque qu'en NRZ pour l'instant !!!!
//ArrayXcf tr_module(const BitStream &bs, sptr<FormeOnde> fo, entier osf)
//{
  //ArrayXcf x = 2 * (bs.vers_array() - 0.5);
  //ArrayXcf x = fo->gene_symboles(bs);
  //assert(x.rows() == bs.lon());
  //msg("tr module : ")
  //si(fo->filtre.)
  //retourne sah(x, osf);
//}


struct TestRecepteurConfig
{
  entier osf = 1;
  sptr<FormeOnde> fo = forme_onde_bpsk();
  float SNR_min = 0;
  float SNR_max = 10;
  bouléen  avec_delais = non;
  bouléen  avec_plot = non;
  entier   nb_rep = 9;
  entier   ncoefs      = 63;//15;
  bouléen  check_errs  = oui;
  bouléen  check_trames = oui;
  bouléen  premiere_trame_SNR_inf = oui;
  entier   nbits = -1;
  float carrier_rec_bl = 0.001;
  float seuil = 0.6;
  entier   ncoefs_filtre_mise_en_forme = 63;

  struct
  {
    bouléen actif = non;
  } cocanal;


  entier nreg_mls = -1;
};


struct TestRecepteurRes
{
  Vecf SNR;
  Vecf ber; // NaN si trame non reçue
  Vecf ber_theo;
  Vecf EbN0;
  Vecf err_phase, err_pos;

  bouléen succès = oui;
};

TestRecepteurRes test_recepteur_unit(const TestRecepteurConfig &config)
{
  TestRecepteurRes res;
  srand(0x124DF531);
  generateur_aleatoire.seed(0x124DF531);

  si(config.avec_plot)
  {
    stdo.flush();
    stdo.def_dossier_sortie(fmt::format("./build/test-log/recepteur-osf-{}-{}-{}", config.osf, *(config.fo), config.avec_delais ? "delais-frac" : "delais-entier"));
  }


  msg_majeur("\n\n\n\n\nTest récepteur, OSF = {}, {}, filtre = {}, {}...\n\n\n\n\n", config.osf, *(config.fo), config.fo->filtre, config.avec_delais ? "delais-frac" : "delais-entier");
  RécepteurConfig rc;

  si(config.avec_plot)
  {
    soit x = config.fo->constellation();
    Figure f;
    f.plot_iq(x, "oa");
    f.afficher(fmt::format("Constellation {}", config.fo->desc_courte()));
  }

  entier nreg_mls = 7 + config.fo->infos.k;

  si(config.nreg_mls >= 0)
    nreg_mls = config.nreg_mls;
  rc.format.entete  = code_mls(nreg_mls); // 127 bits

  // Note : le padding de l'en-tête (multiple de k)
  // est fait automatiquement dans l'émetteur et dans le récepteur

  Vecf full_corr;

  entier nbits = config.nbits == -1 ? 64 : config.nbits;

  msg("Dimension de l'en-tête : {} bits, data : {} bits.", rc.format.entete.lon(), nbits);

  soit fo = config.fo;

  rc.BS                                             = 512;
  rc.format.nbits                                   = nbits;
  //rc.format.fo_entete                               = fo;
  rc.seuil                                          = config.seuil;
  rc.debug_actif                                    = config.avec_plot;
  rc.config_demod.dec.carrier_rec.actif             = non;
  rc.config_demod.dec.clock_rec.actif               = non;
  rc.format.modulation.fsymb                        = 1e3;
  rc.format.modulation.fe                           = rc.format.modulation.fsymb * config.osf;
  rc.format.modulation.fi                           = 0;
  rc.format.modulation.forme_onde                           = fo;
  rc.config_demod.debug_actif                       = config.avec_plot;
  rc.config_demod.dec.carrier_rec.BL                = config.carrier_rec_bl;
  rc.SNR_mini                                       = -10;
  rc.format.modulation.ncoefs_filtre_mise_en_forme  = config.ncoefs_filtre_mise_en_forme;
  si(fo->filtre.type == SpecFiltreMiseEnForme::NRZ)
    rc.format.modulation.ncoefs_filtre_mise_en_forme = config.osf;
  //rc.ncoefs_interpolateur                           = 63;
  rc.callback_corr = [&](const Vecf &corr)
  {
    full_corr = vconcat(full_corr, corr);
  };
  soit rec = récepteur_création(rc);

  ÉmetteurConfig ec;
  ec.format       = rc.format;
  ec.debug_actif  = config.avec_plot;

  soit émet = émetteur_création(ec);

  BitStream data;

  si(config.nbits == -1)
  {
    BitStream data1 = BitStream("1010101011110000");
    data = data1;
    pour(auto i = 0; i < 3; i++)
      data += data1;
  }
  sinon
  {
    data = BitStream::rand(config.nbits);
  }

# if 0
  float retard;
  BitStream bs = rc.format.entete + data;
  ArrayXcf xd = config.fo->génère_échantillons(bs, config.ncoefs, config.osf, retard);
# endif

  float retard = émet->retard();

  soit xd = émet->step(data);

  entier nd = xd.rows();

  retard -= config.osf / 2.0f; // Début du symbole

  msg(" *** retard modulation = {}", retard);

  entier Nt = 24 * 1024;

  /*si(Nt < ((nd + 1024) * config.nb_rep + 1024))
  {
    Nt = (nd + 1024) * config.nb_rep + 1024;
  }*/

  // 1023 * 4 / 3 = 1300
  soit re = rec->get_etat();
  entier nb_echans_padding = re.Ne;//(nbits_entete * config.osf) / config.fo->k;

  msg(" Nb échan padding = {}", nb_echans_padding);


  entier dim_ideale = (nd + 1024) * config.nb_rep + 2 * nb_echans_padding;

  // Ajout de zéros de manière à flusher le détecteur
  Nt = max(Nt, dim_ideale);


  // (8192 est la dimension des blocs d'entrée)
  //si(Nt < ((nd + 1024) * config.nb_rep + 8192 + 2048))
  //{
    //Nt = (nd + 1024) * config.nb_rep + 8192 + 2048;
  //}

  Veccf x = 0.0001 * randn(Nt);

  entier nb_rep = config.nb_rep;//9;

  soit SNR = linspace(config.SNR_min, max(config.SNR_max, config.SNR_min), nb_rep).reverse();

  si(config.premiere_trame_SNR_inf)
    SNR(0) = 200; // Test de référence, (presque) pas de bruit sur la première trame
  std::vector<float> pos(nb_rep);

  res.SNR = SNR;
  res.ber = res.ber_theo = res.EbN0 = Vecf::constant(nb_rep, 1);//;//-1);
  res.err_phase = res.err_pos = Vecf::constant(nb_rep, std::nan(""));

  entier offset = 100;
  pour(auto i = 0; i < nb_rep; i++)
  {


    tsd_assert(offset + nd <= x.rows());

    //float σ = pow(10, - SNR(i) / 20) / sqrt(2.0f);
    // Pareil
    float σ = sqrt((1/db2pow(SNR(i)))) / sqrt(2.0f);



    float τ = 0;
    si(config.avec_delais)
      τ = 0.5;//randu(1)(0);

    pos[i] = offset + retard + τ;


    msg("Placement trame[{}] @ offset {}, SNR={:.1f} dB, σ={}.", i, pos[i], SNR(i), σ);

    Veccf xde;

    si(config.avec_delais)
    {
      xde = tsd::fourier::délais(xd, τ);

      si(config.avec_plot)
      {
        si(tests_debug_actif)
        {
          Figures f;
          f.subplot().plot(xd, "", "Avant retard");
          f.subplot().plot(xde, "", "Après retard τ={}", τ);
          f.afficher("Insertion : avant / après retard");
        }
      }
    }
    sinon
      xde = xd;


    float puissance_trame = abs2(xde).moyenne();
    σ *= sqrt(puissance_trame);
    float SNR_emp = 10 * log10(puissance_trame / (2*σ*σ));
    msg("Energie trame : {}, σ² : {}, SNR empirique = {:.1f} dB.", puissance_trame, 2*σ*σ, SNR_emp);


    Veccf xr(nd);
    xr.set_real(randn(nd) * σ);
    xr.set_imag(randn(nd) * σ);

    float σ_emp     = sqrt(abs2(xr).moyenne());
    float σ_emp_re  = sqrt(abs2(real(xr)).moyenne());
    float σ_emp_im  = sqrt(abs2(imag(xr)).moyenne());

    msg("σ cible = {} (/sqrt(2)={}), σ emp = {}, σ emp re = {}, σ emp im = {}",
        σ * sqrt(2.0f), σ, σ_emp, σ_emp_re, σ_emp_im);

    soit co = Veccf::zeros(nd);

    si(config.cocanal.actif)
    {
      // Tests Pascal :
      //  - QPSK ou pi/4 QPSK : 156.25 ks/s (312 kbit/s)
      //  - ASK 40 KSPS, PRBS9
      //     -> soit 40/156.25 de fsymb
      //  - Décalage ~ 8 kHz
      //     -> soit 8/156.25 de fsymb

      soit fo = forme_onde_ask(2);

      ModConfig cfg_ask = rc.format.modulation;
      cfg_ask.fi    = rc.format.modulation.fsymb * 8.0f / 156.25f;
      cfg_ask.fsymb = rc.format.modulation.fsymb * 40.0f / 156.25f;
      cfg_ask.sortie_reelle = non;

      soit mod = modulateur_création(cfg_ask);

      entier nbits = nd * rc.format.modulation.fe / (rc.format.modulation.fsymb / 10) + 128;
      BitStream bs = BitStream::rand(nbits);
      soit y = mod->step(bs);
      tsd_assert(y.rows() >= nd);

      co = y.head(nd) / 2;

      si(config.avec_plot)
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


    si(config.avec_plot)
    {
      soit xe = x.segment(offset, nd);
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

    si(abs(SNR_emp - SNR(i)) >= 0.5)
    {
      echec("SNR empirique hors borne.");
    }

    offset += nd;
    entier ecart = (rand() & 511) + 128 * 4;
    x.segment(offset, ecart) = bruit_awgn(Veccf::zeros(ecart), σ);
    offset += ecart;
  }

  //x = x.head(offset + 8 * 1024).eval(); // laisse quelques échantillons à la fin pour flusher

  si(config.avec_plot)
  {
    Figures f;
    f.subplot().plot(rc.format.entete.array(), "", "en-tête");
    f.subplot().plot(data.array(), "", "data");
    f.subplot().plot(x, "", "x");
    f.afficher(fmt::format("Test récepteur, osf={}, mod={}", config.osf, *(config.fo)));
  }

  soit trames = rec->step(x);

  msg("Nb trames décodées : {}", trames.size());
  //tsd_assert(!trames.empty());


  si(config.avec_plot)
  {
    Figures fg;
    fg.subplot().plot(x, "", "Vue globale");

    pour(auto j = 0; ((j < (entier) trames.size()) && (j < nb_rep)); j++)
    {
      soit &t = trames[j];
      fg.gcf().plot((float) t.det.position, 1, "sr");
    }

    fg.subplot().plot(full_corr, "", "Corrélation");
    fg.afficher();
  }

  pour(auto j = 0; j < nb_rep; j++)
  {
    float EsN0 = pow2db(config.osf * db2pow(SNR(j)));
    float EbN0 = EsN0 - 10 * log10(fo->infos.k);
    float ber_theo = fo->ber(EbN0);
    res.ber_theo(j) = ber_theo;
    res.EbN0(j)     = EbN0;
  }

  pour(auto j = 0; ((j < (entier) trames.size()) && (j < nb_rep)); j++)
  {
    msg("Vérification trames #{}...", j);
    soit &t = trames[j];
    msg("  détection : {}", t.det);
    //msg("  données décodées : {}", t.bs);


    si(config.avec_plot)
    {
      Figures f;
      f.subplot().plot(t.bs.array(), "o|b", "Données reçues");
      f.subplot().plot(data.array(), "o|g", "Données envoyées");
      f.subplot().plot(data.array() - t.bs.array(), "o|r", "Erreurs");
      f.afficher(fmt::format("Trame décodée vs envoyée (#{})", j));
    }
    si(data.lon() != t.bs.lon())
    {
      msg_avert("Nombre de bits reçu invalide ({} vs {} attendus)", t.bs.lon(), data.lon());
      res.succès = non;
      retourne res;
    }
    soit dst = data.dst_Hamming(t.bs);

    float ber = ((float) dst) / t.bs.lon();

    res.ber(j) = ber;

    float EsN0 = pow2db(config.osf * db2pow(SNR(j)));
    float EbN0 = EsN0 - 10 * log10(fo->infos.k);
    float ber_theo = fo->ber(EbN0);
    msg("  SNR = {:.1f} dB, EsN0 = {:.1f} dB, EbN0 = {:.1f} dB, err=\033[32m{}\033[0m, ber=\033[32m{:.2e}\033[0m, ber theo={:.2e}", SNR(j), EsN0, EbN0, dst, ber, ber_theo);
    soit ep = t.det.position_prec - pos[j];
    msg("  Erreur position = {}.", ep);


    res.err_pos(j) = ep;
    res.err_phase(j) = t.det.θ;

    res.EbN0(j)     = EbN0;

    si(config.check_trames)
    {
      si(abs(ep) >= 0.2)
      {
        msg_avert("Position invalide (attendue : {}, détectée : {}).", pos[j], t.det.position_prec);
        res.succès = non;
        retourne res;
      }
    }

    float nb_err_attendues = ber_theo * data.lon();

    msg("  Nb erreurs = {}, attendues = {}", dst, nb_err_attendues);

    // Tolérance car le BER théorique ne prends pas en compte
    // l'erreur de phase et de gain
    si(config.check_errs && (dst > 5 + nb_err_attendues * 3.0))
    {
      msg_avert(" BER trop élevé (nb err attendues = {:.2f}, nb err détectées = {}).", nb_err_attendues, dst);
      res.succès = non;
      retourne res;
    }
  }

  si(config.check_trames)
  {
    si((entier) trames.size() != nb_rep)
    {
      msg_avert("nombre de trames reçues invalide ({} attendues, {} reçues).", nb_rep, trames.size());
      res.succès = non;
      retourne res;
    }
  }

  //si(config.avec_plot)
  //  stdo.flush();

  res.EbN0 = res.EbN0.reverse().eval();
  res.SNR = res.SNR.reverse().eval();
  res.ber = res.ber.reverse().eval();
  res.ber_theo = res.ber_theo.reverse().eval();
  res.err_phase = res.err_phase.reverse().eval();
  res.err_pos = res.err_pos.reverse().eval();

  retourne res;
}


entier bench_recepteur_a()
{
  soit wf = forme_onde_qpsk();
  wf->filtre = SpecFiltreMiseEnForme::rcs(0.25);


  /*{
    test_recepteur_unit(
    {
      .osf                    = 4,
      .fo                     = wf,
      .SNR_min                = -4,
      .SNR_max                = 8,
      .avec_delais            = non,
      .avec_plot              = oui,
      .nb_rep                 = 13,
      .check_errs             = non,
      .check_trames           = non,
      .premiere_trame_SNR_inf = non,
      .nbits                  = 32,
      .carrier_rec_bl         = 0,
      .seuil                  = 0.4,
    });
  }
  retourne 0;*/



  soit res = test_recepteur_unit(
  {
    .osf                    = 4,
    .fo                     = wf,
    .SNR_min                = -4,
    .SNR_max                = 8,
    .avec_delais            = non,
    .avec_plot              = non,
    .nb_rep                 = 13,
    .check_errs             = non,
    .check_trames           = non,
    .premiere_trame_SNR_inf = non,
    .nbits                  = 16 * 1024,
    .carrier_rec_bl         = 0,
    .seuil                  = 0.4,
    .ncoefs_filtre_mise_en_forme  = 63
  });

  si(tests_debug_actif)
  {
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
      f.def_rdi({res.EbN0.minCoeff() - 1, 1e-6, res.EbN0.valeur_max() + 1, 1.0});
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
  }




  retourne 0;
}



entier bench_recepteur()
{
  bench_recepteur_a();
  //retourne 0;
  soit lst_m =
  {
      forme_onde_bpsk(),
      forme_onde_qpsk(),
      forme_onde_π4_qpsk(),
      forme_onde_psk(8),
      forme_onde_qam(16),
      forme_onde_fsk(4, 1.0f), // Index de 1.0f en 4FSK, pour meilleure discrimination
  };
  //soit lst_f  = {SpecFiltreMiseEnForme::nrz(), SpecFiltreMiseEnForme::srrc(0.5), SpecFiltreMiseEnForme::srrc(0.2), SpecFiltreMiseEnForme::gaussien(0.8)};
  entier osf     = 4;

  FILE *fo = fopen("./build/test-log/bench-recepteur.txt", "wt");

  bouléen first = oui;

  pour(auto m: lst_m)
  {
    m->filtre = m->infos.est_fsk ? SpecFiltreMiseEnForme::gaussien(2.0) : SpecFiltreMiseEnForme::rcs(0.25);
    soit res = test_recepteur_unit(
    {
      .osf                    = osf,
      .fo                     = m,
      .SNR_min                = -4,
      .SNR_max                = 16,
      .avec_delais            = non,
      .avec_plot              = non,
      .nb_rep                 = 11,
      .check_errs             = non,
      .check_trames           = non,
      .premiere_trame_SNR_inf = non,
      .nbits                  = 16 * 1024,
      .carrier_rec_bl         = 0,
      .ncoefs_filtre_mise_en_forme  = 63
    });

    si(first)
    {
      std::string s = "Modulation ; SNR";
      pour(auto i = 0; i < res.SNR.rows(); i++)
        s += fmt::format(" ; {:.1f} dB", res.SNR(i));
      fprintf(fo, "%s\n", s.c_str());
      first = non;
    }

    std::string s = fmt::format("{}", *m);

    s += " ; Eb/N0";
    pour(auto i = 0; i < res.ber.rows(); i++)
      s += fmt::format(" ; {:.1f} dB", res.EbN0(i));
    s += "\n";

    s += " ; Ber théo";
    pour(auto i = 0; i < res.ber.rows(); i++)
      s += fmt::format(" ; {:.4f} %", 100 * res.ber_theo(i));
    s += "\n";

    s += " ; Ber simulé";
    pour(auto i = 0; i < res.ber.rows(); i++)
    {
      si(res.ber(i) >= 0)
        s += fmt::format(" ; {:.4f} %", 100 * res.ber(i));
      sinon
        s += " ; n/d";
    }
    fprintf(fo, "%s\n", s.c_str());
  }

  fclose(fo);

  retourne 0;
}

// Code = MLS 31 bits
// BPSK
// Signal 1.25 MHz
// sous-échantilloné à 250 kHz après filtre adapté NRZ
entier test_recepteur()
{
  soit lst_m =
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
    cfg.avec_plot = oui;
    test_recepteur_unit(cfg);
  }


  {
    TestRecepteurConfig cfg;
    cfg.osf = 4;
    cfg.fo = forme_onde_π4_qpsk();
    //cfg.fo->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    cfg.SNR_min = 4;
    cfg.nb_rep = 2;
    cfg.avec_plot = oui;
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
    cfg.avec_plot = oui;
    cfg.cocanal.actif = oui;
    cfg.nreg_mls = 5;
    cfg.nbits     = 512;
    cfg.check_trames = non;
    //cfg.
    test_recepteur_unit(cfg);


    cfg.fo = forme_onde_π4_qpsk();
    cfg.fo->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    test_recepteur_unit(cfg);
  }
  exit(0);
# endif

  // Ceci fonctionnne !!!!
  //tsd_assert(test_recepteur_unit({.osf = 3, .fo = waveform_bpsk(), .SNR_min = 4, .avec_delais = non, .avec_plot = oui, .nb_rep = 1}).succès);
  //tsd_assert(test_recepteur_unit({.osf = 3, .fo = waveform_bpsk(), .SNR_min = 4, .avec_delais = oui, .avec_plot = oui, .nb_rep = 1}).succès);

  /*{
    soit m = waveform_bpsk();
    m->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    tsd_assert(test_recepteur_unit({.osf = 3, .fo = m, .SNR_min = 4, .avec_delais = non, .avec_plot = oui, .nb_rep = 1}).succès);
  }*/
  /*{
    soit fo = waveform_fsk(4);
    fo->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    test_recepteur_unit({.osf = 4, .fo = fo, .SNR_min = 15, .avec_delais = non, .avec_plot = oui});
  }*/

  /*
  test_recepteur_unit({.osf = 4, .fo = waveform_bpsk(), .SNR_min = 4, .avec_delais = non, .avec_plot = oui, .nb_rep = 1});
  test_recepteur_unit({.osf = 4, .fo = waveform_fsk(),  .SNR_min = 10, .avec_delais = non, .avec_plot = oui, .nb_rep = 1});
  soit fo = waveform_fsk();
  fo->filtre = SpecFiltreMiseEnForme::gaussien(0.8);
  test_recepteur_unit({.osf = 4, .fo =fo,  .SNR_min = 10, .avec_delais = non, .avec_plot = oui, .nb_rep = 1});*/
  //test_recepteur_unit({.osf = 2, .fo = waveform_psk(8), .SNR_min = 4, .avec_delais = non, .avec_plot = oui});

  //soit lst_f = {SpecFiltreMiseEnForme::nrz(), SpecFiltreMiseEnForme::rcs(0.5)};
  soit lst_f = {SpecFiltreMiseEnForme::rcs(0.5), SpecFiltreMiseEnForme::nrz()};

  pour(auto m: lst_m)
  {



    pour(auto f: lst_f)
    {
      m->filtre = f;
      pour(auto osf: {4, 2, 3, 4, 5, 6})
        pour(auto avec_delais: {non, oui})
        {
          float SNR_min = 4;
          si(m->infos.est_fsk)
            SNR_min = 15;
          soit res = test_recepteur_unit({.osf = osf, .fo = m, .SNR_min = SNR_min, .avec_delais = avec_delais, .avec_plot = non});
          stdo.def_dossier_sortie("./build/test-log/recepteur");
          si(!res.succès)
          {
            test_recepteur_unit({.osf = osf, .fo = m, .SNR_min = SNR_min, .avec_delais = avec_delais, .avec_plot = oui});
            echec("Le test unitaire du récepteur a échoué : osf={}, fo={}", osf, *m);
            retourne -1;
          }
        }
    }
  }


  retourne 0;
}


entier test_filtre_adapte()
{

  pour(auto i = 0; i < 2; i++)
  {
    soit sf = SpecFiltreMiseEnForme::rcs(0.5);
    entier osf = 2;

    soit fmef = sf.filtre_mise_en_forme(15, osf);
    soit fa   = sf.filtre_adapté(15, osf);

    entier n = 30;
    soit x = Veccf::zeros(n);
    si(i == 0)
    {
      pour(auto i = 0; i + 1 < n; i += 2)
      {
        x(i)   = -1;
        x(i+1) = 1;
      }
    }
    sinon
      x(n/2) = 1;

    soit y = fmef->step(x);
    soit z = fa->step(y);


    soit itrp = tsd::filtrage::itrp_sinc<cfloat>({15, 1024, 0.45, "hn"});

    soit c0 = itrp->coefs(0);
    soit c1 = itrp->coefs(0.5);

    si(tests_debug_actif)
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

    soit fc0 = filtre_rif<float, cfloat>(c0),
         fc1 = filtre_rif<float, cfloat>(c1);

    soit z0 = fc0->step(z), z1 = fc1->step(z);

    si(tests_debug_actif)
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
  retourne 0;
}




// TODO !!!
void test_fsk()
{
  // 16 bits aléatoires
  soit bs = randstream(16);

  ModConfig cfg;
  cfg.forme_onde            = forme_onde_fsk(2, 4);//2, index, filt, BT);
  cfg.debug_actif   = oui;
  cfg.fe            = 10e3;
  cfg.fi            = 500;
  cfg.fsymb         = 0.1e3;
  cfg.sortie_reelle = non;

  // Fréquence d'échantillonnage  = 10 kHz
  // Fréquence intermédiaire      = 500 Hz
  // Fréquence symbole            = 0.1 kHz
  // -> Déviation = fsymb * h / 2 = 200 Hz
  soit mod = modulateur_création(cfg);

  soit x = mod->step(bs);

  soit discri = discriminateur_fm();
  soit y = discri->step(x);

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot().plot(bs.array(), "hb", "Train binaire");
    f.subplot().plot(x, "", "Modulation FSK");
    f.subplot().plot(y.tail(y.rows() * 0.9), "", "Discri");
    f.afficher("Test FSK");
  }
}

entier test_demod()
{
  msg_majeur("Tests démodulation...");

  soit filtres = {SpecFiltreMiseEnForme::nrz(), SpecFiltreMiseEnForme::rcs(0.4)};
  soit Ms      = {2, 4, 8};

  std::vector<sptr<FormeOnde>> wfs;

  // TODO : ne fonctionne plus !!!
  /*wfs.push_back(waveform_fsk(2, 0.5,  SpecFiltreMiseEnForme::nrz()));
  wfs.push_back(waveform_fsk(4, 0.5,  SpecFiltreMiseEnForme::nrz()));
  wfs.push_back(waveform_fsk(2, 0.5,  SpecFiltreMiseEnForme::gaussien(1.5)));
  wfs.push_back(waveform_fsk(2, 1,    SpecFiltreMiseEnForme::nrz()));
  wfs.push_back(waveform_fsk(2, 0.7,  SpecFiltreMiseEnForme::nrz()));
  wfs.push_back(waveform_fsk(2, 2,    SpecFiltreMiseEnForme::nrz()));*/

  pour(auto filtre: filtres)
  {
    pour(auto M: Ms)
      wfs.push_back(forme_onde_psk(M, filtre));
    wfs.push_back(forme_onde_qam(16, filtre));
  }

  pour(auto wf: wfs)
  {
    msg_majeur("Test démodulation {}", *wf);
    stdo.printf(fmt::format("<h2>Test démodulation {}</h2>", *wf));

    ModConfig modcfg;
    DemodConfig cfg;
    modcfg.forme_onde  = wf;
    modcfg.fe          = 100e3;
    modcfg.fi          = 0;
    modcfg.fsymb       = 20e3;
    modcfg.debug_actif = tests_debug_actif;
    cfg.debug_actif    = tests_debug_actif;

    cfg.dec.clock_rec.tc = 10;
    cfg.dec.cag.actif = non;

    soit demod         = démodulateur_création(modcfg, cfg);

    // Create a modulator to simulate RF signal
    soit mod = modulateur_création(modcfg);

    entier nbits = 200;//1000

    soit bs = randstream(nbits);
    soit x = mod->step(bs);

    x = bruit_awgn(x, 0.001);

    // Proceed to demodulation
    BitStream bs2;
    Tabf llr;
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

    // Ignore les premiers bits (convergence des boucles)
    bs  = BitStream(bs.array().tail((nbits*7)/10));
    bs2 = BitStream(bs2.array().tail((nbits*7)/10));

    // Alignment of the two bit vectors (and phase ambiguity resolution)
    soit res = cmp_bits_psk(bs, bs2, wf->infos.k);

    si(tests_debug_actif)
    {
      Figures f;
      f.subplot().plot(res.b0, "hb", "Train binaire émis");
      f.subplot().plot(res.b1, "hg", "Train binaire décodé");
      f.subplot().plot(res.b1 - res.b0, "hr", "Erreurs (nb errs = {}, ber = {:.1e})", res.nerr, res.ber);
      f.afficher("Bits entrées / sorties");
    }



    // Verification
    si(res.nerr != 0)
      echec("Bits erronés non attendus !");
  }
  retourne 0;
}



void test_sah()
{
  soit x = Vecf::valeurs({0, 1, 2});
  soit y = sah(x, 3);
  soit yref = Vecf::valeurs({0, 0, 0, 1, 1, 1, 2, 2, 2});
  soit err = abs(y - yref).valeur_max();
  tsd_assert_msg(err == 0, "Erreur SAH");
  msg("sah ok.");
}


static void test_émetteur()
{
  msg_majeur("Test émetteur...");
  RécepteurConfig rc;
  rc.format.entete            = BitStream::rand(127);
  rc.format.modulation.forme_onde     = forme_onde_bpsk();
  rc.format.modulation.fe     = 1e6;
  rc.format.modulation.fsymb  = 1e5;
  rc.format.nbits             = 256;
  rc.debug_actif              = oui;

  ÉmetteurConfig ec;
  ec.debug_actif = oui;
  ec.format = rc.format;

  msg("Création émetteur...");
  soit eme = émetteur_création(ec);
  msg("Création récepteur...");
  soit rec = récepteur_création(rc);

  soit bs = BitStream::rand(256);

  msg("Emission...");
  soit x = eme->step(bs);

  x = vconcat(Veccf::zeros(1000), x);
  x = vconcat(x, Veccf::zeros(40e3));

  msg("Réception...");
  soit tr = rec->step(x);

  msg("Reçu : {} trames.", tr.size());

  tsd_assert(tr.size() == 1u);
  tsd_assert(tr[0].bs == bs);

  msg("ok.");
}

entier test_telecom()
{
  test_filtre_boucle_ordre_1();
  test_émetteur();
  test_filtre_adapte();
  test_discri_fm();
  test_fsk();
  test_demod();
  test_sah();
  retourne 0;
}
