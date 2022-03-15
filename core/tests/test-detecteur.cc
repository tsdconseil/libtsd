#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


struct TestMotifConfig
{
  std::string nom;
  ArrayXcf motif;
  float σ = 0.5;
  int ntests = 500;
};

struct TestMotifResult
{
  float σ_gain, σ_phase, σ_temps;
};

static TestMotifResult test_motif(const TestMotifConfig &config)
{
  int cnt_ech = 0;
  TestMotifResult res;
  ArrayXcf m =  config.motif / sqrt(config.motif.abs2().mean());

  //msg("Test motif[{}] : {} échantillons.", config.nom, config.motif.rows());


  if(tests_debug_actif)
  {
    Figure f;
    f.plot(config.motif);
    f.afficher(config.nom);
    auto [l, c] = xcorrb(config.motif, config.motif);
    Figure f2;
    f2.plot(l, c.abs(), "", "xcorrb");
    f2.afficher();
  }

  int M = m.rows();
  int N = 40 * M;

  float egain = 0, ephase = 0, etemps = 0;

  int ndet = 0;

  DetecteurConfig dconfig;
  dconfig.debug_actif    = false;
  dconfig.gere_detection = [&](const Detection &det)
  {
    ndet++;
    etemps += carré(det.position_prec + cnt_ech - 10 * M);
    ephase += carré(rad2deg(det.θ));
    egain  += carré(det.gain - 1);
  };


  int BS = 1024;
  dconfig.Ne    = BS;
  dconfig.motif = m;
  dconfig.seuil = 0.8;
  auto det = détecteur_création(dconfig);

  for(auto idt = 0; idt < config.ntests; idt++)
  {
    ndet = 0;
    cnt_ech = 0;
    ArrayXcf x = ArrayXcf::Zero(N);
    //ArrayXcf x = config.σ * randn(N);
    x.segment(10 * M, M) += m;
    //x = bruit_awgn(x, config.σ);


    ArrayXcf bruit = (config.σ / sqrt(2.0f)) * (randn(N) + cfloat(0,1) * randn(N));

    x += bruit;

    //det->step(x);
    for(auto i = 0; i < N / BS; i++)
    {
      det->step(x.segment(i * BS, BS));
      cnt_ech += BS;
    }
    if(ndet != 1)
    {
      echec("Motif non détecté (ndet = {}).", ndet);
    }
  }

  res.σ_gain  = sqrt(egain  / (config.ntests - 1));
  res.σ_phase = sqrt(ephase / (config.ntests - 1));
  res.σ_temps = sqrt(etemps / (config.ntests - 1));

  msg("Test motif[{}, {} échans] : \033[32mσ_gain={:.2e}, σ_phase={:.2e}°, σ_temps={:.2e}\033[0m",
      config.nom, config.motif.rows(), res.σ_gain, res.σ_phase, res.σ_temps);

  return res;
}



/*static BitStream pad_bs(const BitStream &bs, int nb_bits_par_symb)
{
  BitStream bs2 = bs;
  int nbits = bs.lon();
  int nzeros = nb_bits_par_symb - (nbits % nb_bits_par_symb);
  for(auto i = 0; i < nzeros; i++)
    bs2.push(0);
  return bs2;
}*/


int test_motifs()
{
  {
    ModConfig mconfig;
    mconfig.forme_onde = forme_onde_bpsk();
    mconfig.fe    = 4;
    mconfig.fsymb = 1;
    mconfig.forme_onde->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    auto mod = modulateur_création(mconfig);
    auto bs = code_mls(4);//BitStream::rand(15);

    //ArrayXcf x = mod->step(pad_bs(bs, mconfig.wf->infos.k));
    bs.pad_mult(mconfig.forme_onde->infos.k);
    ArrayXcf x = mod->step(bs);
    test_motif({"bpsk-15bits", x});


    bs = BitStream::zéros(15);
    for(auto i = 0; i < 15; i += 2)
      bs.set(i, 1);

    bs.pad_mult(mconfig.forme_onde->infos.k);
    x = mod->step(bs);
    test_motif({"bpsk-15bits-010101", x});


    bs = BitStream::zéros(7) + BitStream::uns(8);
    bs.pad_mult(mconfig.forme_onde->infos.k);
    x = mod->step(bs);
    test_motif({"bpsk-15bits-000000011111111", x});


    mconfig.forme_onde = forme_onde_qpsk();
    mconfig.forme_onde->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    mod = modulateur_création(mconfig);
    bs = code_mls(5);//BitStream::rand(31);
    bs.pad_mult(mconfig.forme_onde->infos.k);
    x = mod->step(bs);
    test_motif({"qpsk-31bits", x});


    mconfig.forme_onde = forme_onde_π4_qpsk();
    mconfig.forme_onde->filtre = SpecFiltreMiseEnForme::srrc(0.5);
    mod = modulateur_création(mconfig);
    //bs = BitStream::rand(31);
    bs.pad_mult(mconfig.forme_onde->infos.k);
    x = mod->step(bs);
    test_motif({"π4qpsk-31bits", x});
  }


  //test_motif({"essai", siggauss(100) * sigchirp2(0.001, 0.2, 100)});
  return 0;
}





void test_detecteur_unit(float σ, int BS, DetecteurConfig::Mode mode)
{
  int err = 0;
  msg_majeur("Test détecteur de motif: σ={}, BS={}, mode={}...",
              σ, BS, mode == DetecteurConfig::Mode::MODE_OLA ? "OLA" : "RIF");

  int cnt_ech = 0;
  int N = 8 * BS;

  // Dim motif
  int M = 400;
  ArrayXcf motif = siggauss(M) * sigchirp2(0.001, 0.2, M);

  //ArrayXcf motif = tsd::telecom::code_Barker(13).vers_array();
  //int M = motif.rows();

  struct Occurence
  {
    float gain, position, phase;
  };


  std::vector<Occurence> occurences = {
      {2.0f, 900.0f,   π_f/4},
      {4.0f, 2000.4f, -π_f/4},
      // Pour tester les bords...
      {1.0f, BS - 1.0f, 0},
      {1.0f, 2.0f*BS, 0},
      // Pour tester le calcul de SNR
      {0.1f, 2.2f*BS, 0},
      {0.05f, 2.5f*BS, 0},
      {0.02f, 2.7f*BS, 0},
  };

  msg("Occurences attendues :");
  for(auto &o: occurences)
    msg("  pos={:.1f}, gain={}, phase={}°", o.position, o.gain, rad2deg(o.phase));

  // On se débrouille pour que le motif soit d'énergie = 400 (juste pour aider au déboggage)
  motif *= std::sqrt(400.0f) / sqrt(motif.abs2().sum());
  msg("Energie motif = {} (devrait être 400)", motif.abs2().sum());

  // Place les deux motifs sur un vecteur de bruit
  ArrayXcf x = ArrayXcf::Zero(N);

  for(auto &o: occurences)
  {
    ArrayXcf motif2 = motif;
    float p = o.position - floor(o.position);
    if(p > 0.01)
      motif2 = tsd::fourier::délais(motif, p);
    x.segment(floor(o.position), M)  = motif2 * std::polar(o.gain, o.phase);
  }

  ArrayXcf bruit = (σ / sqrt(2.0f)) * (randn(N) + cfloat(0,1) * randn(N));

  float σ_emp = std::sqrt(bruit.abs2().mean());

  msg("σ = {}, σ_emp = {}.", σ, σ_emp);

  tsd_assert(std::abs(σ - σ_emp) / σ  < 1e-2);

  x += bruit;

  if(tests_debug_actif)
  {
    auto [lags, xc] = xcorr(motif);
    Figures f;
    f.subplot().plot(motif);
    f.subplot().plot(lags, xc.real(), "", "Auto-corrélation");
    f.afficher("Motif");
  }

  std::vector<float> dets, scores;

  std::vector<Occurence> occd;

  DetecteurConfig config;
  config.mode = mode;
  config.debug_actif = true;
  config.gere_detection = [&](const Detection &det)
  {
    msg("Detection : {}.", det);
    auto pos_abs = det.position_prec + cnt_ech;
    msg("Pos abs = {}", pos_abs);
    dets.push_back(pos_abs);
    scores.push_back(det.score);

    float theta_v = 0, gain_v = 0, pos_v = 0, SNR_v = 0;

    int n = dets.size() - 1;
    if(n < (int) occurences.size())
    {
      theta_v = occurences[n].phase;
      gain_v  = occurences[n].gain;
      pos_v   = occurences[n].position;
      SNR_v   = pow2db((motif * occurences[n].gain).abs2().mean() / (σ*σ));
    }

    msg("  détection attendue : θ={}°, g={}, p={}, σ={}, SNR={:.1f} dB.", rad2deg(theta_v), gain_v, pos_v, σ, SNR_v);
    auto err_phase = rad2deg(det.θ - theta_v);
    auto err_gain  = det.gain - gain_v;
    auto err_pos   = pos_abs - pos_v;
    auto err_σ     = det.σ_noise - σ;
    auto err_SNR   = det.SNR_dB - SNR_v;
    msg("  erreurs : \033[32mθ={:.1e}°, g={:.1e}, p={:.1e}, σ={:.1e}, SNR={:.1f} dB\033[0m.",
        err_phase, err_gain, err_pos, err_σ, err_SNR);



    if(SNR_v > 15)
    {
      if(abs(err_phase) > 1)
      {
        err = 1;
        msg_avert("Erreur de phase trop importante");
      }

      if(abs(err_gain / gain_v) > 1e-2)
      {
        err = 1;
        msg_avert("Erreur relative de gain trop importante");
      }
    }

    if(abs(err_pos) > 0.1)
    {
      err = 1;
      msg_avert("Erreur de position trop importante");
    }

    if(abs(err_σ / σ) > 0.25)
    {
      err = 1;
      msg_avert("Erreur relative σ trop importante");
    }

    if(abs(err_SNR) > 1)
    {
      err = 1;
      msg_avert("Erreur SNR trop importante");
    }
  };


  config.Ne    = BS;
  config.motif = motif;
  config.seuil = 0.8;
  auto det = détecteur_création(config);

  //det->step(x);

  for(auto i = 0; i < N / BS; i++)
  {
    det->step(x.segment(i * BS, BS));
    cnt_ech += BS;
  }

  if(tests_debug_actif)
  {
    Figure f;
    f.plot(x, "", "x");
    ArrayXf t = Eigen::Map<ArrayXf>(dets.data(), dets.size());
    ArrayXf z = Eigen::Map<ArrayXf>(scores.data(), scores.size());
    f.plot(t, z, "or", "Détections");
    f.canva().set_couleur(tsd::vue::Couleur::Rouge);
    for(auto i = 0; i < t.rows(); i++)
      f.canva().ligne(t(i), z(i), t(i) + M, z(i));
    f.afficher();
  }

  tsd_assert(dets.size() == occurences.size());
  //for(auto &s: scores)
    //tsd_assert(std::abs(s - 1) < 0.02);
  if(err)
    echec("Erreur trop importante.");
  msg("Fin.");
}

int test_detecteur()
{

  test_detecteur_unit(0.01, 4*1024, DetecteurConfig::MODE_OLA);
  test_detecteur_unit(0.01, 4*1024, DetecteurConfig::MODE_RIF);

  return 0;
}



