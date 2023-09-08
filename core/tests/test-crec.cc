#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


void test_ped(cstring nom, Ped ped)
{
  msg("Vérification détecteur d'erreur de phase [{}]...", nom);

  soit N = 512;
  Veccf x(N);
  soit f = 0.02f;
  soit ϕ_moy = 0.0;
  pour(auto i = 0; i < N; i++)
  {
    x(i) = exp(2 * π_f * ⅈ * f);
    soit ϕ = ped(x(i));
    ϕ_moy += ϕ / N;
  }

  soit f_det = ϕ_moy / (2 * π);

  msg("  f = {}, détecté = {}, erreur = {}.", f, f_det, abs(f - f_det));

  verifie_erreur_relative(f, f_det, 10.0, "det freq");
}






void test_clkrec()
{
  msg_majeur("Test recouvrement d'horloge...");

  ModConfig mc;
  mc.forme_onde     = forme_onde_psk(2);
  mc.debug_actif    = oui;
  mc.fe             = 800e3;//1e6;
  mc.fi             = 0;
  mc.fsymb          = 200e3;
  mc.sortie_reelle  = non;

  soit mod = modulateur_création(mc);

  soit f = filtre_mg<cfloat,cdouble>(mc.fe / mc.fsymb);

  BitStream xb = randstream(70);
  soit y  = mod->step(xb),
       y2 = tsd::fourier::délais(y, 8.4),
       y3 = f->step(y2);

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
    f.subplot().plot(real(y),  "ob-", "sortie modulateur");
    f.subplot().plot(real(y2), "ob-", "après délais fractionnaire");
    f.subplot().plot(real(y3), "ob-", "après filtre adapté");
    f.subplot().plot(real(y4), "ob-", "après clkrec");
    f.afficher("Synthèse clk rec");
  }
}



void test_crec()
{
  msg_majeur("Test recouvrement de porteuse (CPLL)...");

  test_ped("ploop - M = 2", ped_ploop(2));
  test_ped("tloop - M = 2", ped_tloop(2));
  test_ped("costa - M = 2", ped_costa(2));

  PLLConfig config;
  config.debug              = oui;
  config.freq               = 0;
  config.loop_filter_order  = 2;
  //config.ped = ped_ploop(2);
  config.ped                = ped_tloop(2);
  config.tc                 = 10;

  ModConfig mc;
  mc.forme_onde     = forme_onde_psk(2);
  mc.debug_actif    = oui;
  mc.fe             = 1e6;
  mc.fi             = 0;
  mc.fsymb          = 100e3;
  mc.sortie_reelle  = non;
  soit mod = modulateur_création(mc);

  soit xb = randstream(1000);
  soit y  = bruit_awgn(mod->step(xb), 0.1) * exp(π_f * ⅈ / 4.0f);

  soit n = y.rows();

  soit df = 0.01f;

  soit ol = source_ohc(df);
  y *= ol->step(n);

  msg("Doppler : fréq = {}, pulsation = {} degrés.", df, df * 360);

  soit pll = cpll_création(config);

  soit yc = pll->step(y);

  soit am = rad2deg(arg(square(yc)).moyenne() / 2);

  msg("Erreur de phase moyenne en sortie : {} degrés.", am);

  assertion(abs(am) < 10);

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
}



