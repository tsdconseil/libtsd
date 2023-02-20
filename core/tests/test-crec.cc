#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


entier test_ped(const std::string &nom, Ped ped)
{
  msg("Vérification détecteur d'erreur de phase [{}]...", nom);
  //tsd::telecom::DEPQuadratique det;
  //det.configure(2.0);

  soit N = 512u;
  soit x = Veccf::zeros(N);
  soit f = 0.02f;
  soit ϕ_moy = 0.0;
  pour(auto i = 0u; i < N; i++)
  {
    x(i) = std::polar(1.0f, 2 * π_f * f);
    soit ϕ = ped(x(i));
    ϕ_moy += ϕ / N;
  }

  soit f_det = ϕ_moy / (2 * π);

  msg("  f = {}, détecté = {}, erreur = {}.", f, f_det, abs(f - f_det));

  retourne verifie_erreur_relative(f, f_det, 10.0, "det freq");
}






entier test_clkrec()
{
  entier res = 0;

  msg_majeur("Test recouvrement d'horloge...");

  ModConfig mc;
  mc.forme_onde             = forme_onde_psk(2);
  mc.debug_actif    = oui;
  mc.fe             = 800e3;//1e6;
  mc.fi             = 0;
  mc.fsymb          = 200e3;
  mc.sortie_reelle  = non;
  soit mod = modulateur_création(mc);

  //ArrayXf xb = randb(70);
  //ArrayXf xb = randb(150);

  BitStream xb = randstream(70);
  soit y = mod->step(xb);


  soit y2 = tsd::fourier::délais(y, 8.4);
  //ArrayXcf y2 = tsd::fourier::delais_entier(y, 3);

  soit f = filtre_mg<cfloat,cdouble>(mc.fe / mc.fsymb);
  soit y3 = f->step(y2);

 // y = bruit_awgn(y, 0.1);

  ClockRecConfig config;
  config.debug_actif  = oui;
  //config.itrp         = itrp_cspline<cfloat>();
  config.itrp         = itrp_lineaire<cfloat>();
  config.ted          = ted_init(TedType::GARDNER);
  config.tc           = 3;
  config.osf          = mc.fe / mc.fsymb;
  soit cr = clock_rec_init(config);

  msg("appel clkrec...");
  soit y4 = cr->step(y3);
  msg("ok.");

  si(tests_debug_actif)
  {
    Figures f;
    f.subplot(411).plot(real(y), "ob-", "sortie modulateur");
    f.subplot(412).plot(real(y2), "ob-", "après délais fractionnaire");
    f.subplot(413).plot(real(y3), "ob-", "après filtre adapté");
    f.subplot(414).plot(real(y4), "ob-", "après clkrec");
    f.afficher("Synthèse clk rec");
  }

  retourne res;
}



entier test_crec()
{
  entier res = 0;
  msg_majeur("Test recouvrement de porteuse (CPLL)...");

  res |= test_ped("ploop - M = 2", ped_ploop(2));
  res |= test_ped("tloop - M = 2", ped_tloop(2));
  res |= test_ped("costa - M = 2", ped_costa(2));
  // res |= test_ped("costa - M = 4", ped_costa(4)); // ne marche pas

  PLLConfig config;
  config.debug = oui;
  config.freq = 0;
  config.loop_filter_order = 2;
  //config.ped = ped_ploop(2);
  config.ped = ped_tloop(2);
  config.tc = 10;

  ModConfig mc;
  mc.forme_onde = forme_onde_psk(2);
  mc.debug_actif    = oui;
  mc.fe             = 1e6;
  mc.fi             = 0;
  mc.fsymb          = 100e3;
  mc.sortie_reelle  = non;
  soit mod = modulateur_création(mc);

  //ArrayXf xb = randb(1000);
  BitStream xb = randstream(1000);
  soit y = mod->step(xb);

  y = bruit_awgn(y, 0.1);

  y *= std::polar(1.0f, π_f / 4);
  entier n = y.rows();

  float df = 0.01;//0.001;

  soit ol = source_ohc(df);
  y *= ol->step(n);

  msg("Doppler : fréq = {}, pulsation = {} degrés.", df, df * 360);

  soit pll = cpll_création(config);

  soit yc = pll->step(y);

  soit am = rad2deg(arg(square(yc)).moyenne() / 2);

  msg("Erreur de phase moyenne en sortie : {} degrés.", am);

  tsd_assert(abs(am) < 10);

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


  retourne res;
}



