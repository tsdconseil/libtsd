#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"


namespace tsd::filtrage {

 extern std::tuple<Vecf, Vecf> rifampn(const Vecf &h, entier L, bouléen symetrique);
}


static void test_plotfiltre()
{
  soit h = design_lexp(0.5);
  msg("h = {}", h);
  soit f = plot_filtre(h, oui) ;
  f.afficher() ;

  soit ri = repimp(h, 20);
  msg("ri = {}", ri);
}

static void test_rifamp()
{
  msg("Test rifamp...");
  //soit h = Vecf::valeurs({0.1, 1, 0.1});
  //soit h = design_hb(31, 0.25-0.10);

  soit h = design_rif_demi_bande(31, 0.25-0.10);


  soit [fr, amp] = rifamp(h, 128, oui);
  soit [fr2, mag] = frmag(h, 128);

  soit err = mag - abs(amp);

  Figures f;
  f.subplot().plot(fr,  amp, "", "rifamp");
  f.subplot().plot(fr2, mag, "", "frmag");
  f.subplot().plot(fr2, err, "r-", "erreur");
  f.afficher("rifamp");
}


static void test_frmag()
{
  msg("Test frmag & frphase...");
  soit n = 100, m = 233;
  soit x = sigcar(10, n);

  soit [fr, mag] = frmag(x, m);
  soit [f2, ph]  = frphase(x, m);
  tsd_assert((fr.rows() == m) && (f2.rows() == m) && (mag.rows() == m));


  Vecf vmag(m), vphase(m);
  pour(auto i = 0; i < m; i++)
  {
    soit r = (x * tsd::polar(-linspace(0, n-1, n) * (fr(i) * 2 * π_f))).somme();
    // Théoriquement :
    // mag(Omega) = somme h_k e^{i Omega k}
    vmag(i) = abs(r);
    vphase(i) = arg(r);
  }
  si(tests_debug_actif)
  {
    Figure f;
    f.plot(fr, mag,  "b-", "mag");
    f.plot(fr, vmag, "g-", "vraie mag");
    f.afficher();
    f.clear();
    f.plot(fr, mag - vmag, "r-", "erreur");
    f.afficher();
  }

  Vecf dh(m);
  for(auto k = 0; k < m; k++)
    dh(k) = modulo_pm_π(ph(k) - vphase(k));

  vphase = déplie_phase(vphase);

  {
    Figure f;
    f.plot(fr, ph,  "b-", "phase");
    f.plot(fr, vphase, "g-", "vraie phase");
    f.afficher();
    f.clear();
    f.plot(fr, dh, "r-", "erreur");
    f.afficher();
  }
  soit err = abs(vmag - mag).valeur_max();
  msg("Erreur max frmag = {}", err);
  tsd_assert(err < 1e-3);


  soit err2 = abs(dh).valeur_max();
  msg("Erreur max frphase = {}", err2);
  tsd_assert(err2 < 1e-2);
}

entier test_filtrage_analyse()
{
  test_plotfiltre();
  test_rifamp();
  test_frmag();
  retourne 0;
}

