#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


int test_ped(const std::string &nom, Ped ped)
{
  msg("Vérification détecteur d'erreur de phase [{}]...", nom);
  //tsd::telecom::DEPQuadratique det;
  //det.configure(2.0);

  auto N = 512u;
  tsd::ArrayXcf x = tsd::ArrayXcf::Zero(N);
  float f = 0.02;
  double ϕ_moy = 0;
  for(auto i = 0u; i < N; i++)
  {
    x(i) = std::polar(1.0f, 2 * π_f * f);
    float ϕ = ped(x(i));
    ϕ_moy += ϕ / N;
  }

  float f_det = ϕ_moy / (2 * π);

  msg("  f = {}, détecté = {}, erreur = {}.", f, f_det, std::abs(f - f_det));

  return verifie_erreur_relative(f, f_det, 10.0, "det freq");
}






int test_clkrec()
{
  int res = 0;

  msg_majeur("Test recouvrement d'horloge...");

  ModConfig mc;
  mc.wf             = forme_onde_psk(2);
  mc.debug_actif    = true;
  mc.fe             = 800e3;//1e6;
  mc.fi             = 0;
  mc.fsymb          = 200e3;
  mc.sortie_reelle  = false;
  auto mod = modulateur_création(mc);

  //ArrayXf xb = randb(70);
  //ArrayXf xb = randb(150);

  BitStream xb = randstream(70);
  ArrayXcf y = mod->step(xb);


  ArrayXcf y2 = tsd::fourier::délais(y, 8.4);
  //ArrayXcf y2 = tsd::fourier::delais_entier(y, 3);

  auto f = filtre_mg<cfloat,cdouble>(mc.fe / mc.fsymb);
  ArrayXcf y3 = f->step(y2);

 // y = bruit_awgn(y, 0.1);

  ClockRecConfig config;
  config.debug_actif  = true;
  //config.itrp         = itrp_cspline<cfloat>();
  config.itrp         = itrp_lineaire<cfloat>();
  config.ted          = ted_init(TedType::GARDNER);
  config.tc           = 3;
  config.osf          = mc.fe / mc.fsymb;
  auto cr = clock_rec_init(config);

  msg("appel clkrec...");
  ArrayXcf y4 = cr->step(y3);
  msg("ok.");

  if(tests_debug_actif)
  {
    Figures f;
    f.subplot(411).plot(y.real(), "ob-", "sortie modulateur");
    f.subplot(412).plot(y2.real(), "ob-", "après délais fractionnaire");
    f.subplot(413).plot(y3.real(), "ob-", "après filtre adapté");
    f.subplot(414).plot(y4.real(), "ob-", "après clkrec");
    f.afficher("Synthèse clk rec");
  }

  return res;
}



int test_crec()
{
  int res = 0;
  msg_majeur("Test recouvrement de porteuse (CPLL)...");

  res |= test_ped("ploop - M = 2", ped_ploop(2));
  res |= test_ped("tloop - M = 2", ped_tloop(2));
  res |= test_ped("costa - M = 2", ped_costa(2));
  // res |= test_ped("costa - M = 4", ped_costa(4)); // ne marche pas

  PLLConfig config;
  config.debug = true;
  config.freq = 0;
  config.loop_filter_order = 2;
  //config.ped = ped_ploop(2);
  config.ped = ped_tloop(2);
  config.tc = 10;

  ModConfig mc;
  mc.wf = forme_onde_psk(2);
  mc.debug_actif    = true;
  mc.fe             = 1e6;
  mc.fi             = 0;
  mc.fsymb          = 100e3;
  mc.sortie_reelle  = false;
  auto mod = modulateur_création(mc);

  //ArrayXf xb = randb(1000);
  BitStream xb = randstream(1000);
  ArrayXcf y = mod->step(xb);

  y = bruit_awgn(y, 0.1);

  y *= std::polar(1.0f, π_f / 4);
  int n = y.rows();

  float df = 0.01;//0.001;

  auto ol = source_ohc(df);
  y *= ol->step(n);

  msg("Doppler : fréq = {}, pulsation = {} degrés.", df, df * 360);

  auto pll = cpll_création(config);

  ArrayXcf yc = pll->step(y);


  auto am = rad2deg(yc.square().arg().mean() / 2);

  msg("Erreur de phase moyenne en sortie : {} degrés.", am);

  tsd_assert(std::abs(am) < 10);

  /*{
    Figure f("Entrée / sortie rec porteuse");
    f.subplot(121);
    f.plot_iq(y, "ba");
    f.titre("Avant CREC");
    f.subplot(122);
    f.plot_iq(yc, "ga");
    f.titre("Après CREC");
    f.afficher();
  }*/


  return res;
}



