#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"
#include "tsd/temps.hpp"
#include <tuple>


using namespace tsd::temps;


static void test_heure_composite()
{
  HeureComposite hc("23:10:02");
  assertion((hc.heure == 23) && (hc.minutes == 10) && (hc.secondes == 2) && (hc.µs == 0));

  vérifie_exception([]() {HeureComposite hc("23/10:02");});
  vérifie_exception([]() {HeureComposite hc("25:10:02");});
  vérifie_exception([]() {HeureComposite hc("ab:10:02");});
}



static void test_cal(const Calendrier &cal)
{
  DateHeure t({cal, {7, 59, 30}});
  soit cal2 = t.calendrier();
  assertion_msg(cal == cal2, "Erreur calendrier : {} VS {}", cal, cal2);
}

// A compléter !!!
void test_date_heure()
{
  test_heure_composite();

  {
    msg("Test calendrier...");
    Calendrier c(2022, 01, 01);
    assertion(c.est_valide());
    assertion(c.nb_jours_debut_année() == 0);
    c.jour = 10;
    assertion(c.nb_jours_debut_année() == 9);
    c.année = 1;
    assertion_msg(c.nb_jours_debut_ère() == 9, "nb = {}", c.nb_jours_debut_ère());
    c.année = 2;
    assertion_msg(c.nb_jours_debut_ère() == 9 + 365, "nb = {}", c.nb_jours_debut_ère());
  }

  assertion(est_bissextile(2020));
  assertion(!est_bissextile(2021));
  assertion(!est_bissextile(1900));
  assertion(est_bissextile(2000));
  assertion(!est_bissextile(2100));
  assertion(!est_bissextile(2200));
  assertion(!est_bissextile(2300));
  assertion(est_bissextile(2400));

  // Tests conversion Calendrier <-> DateHeure
  pour(auto année : {2000, 2020, 2021, 2022, 2023})
  {
    test_cal({année, 01, 01});
    test_cal({année, 01, 25});
    test_cal({année, 02, 28});
    test_cal({année, 05, 31});
    test_cal({année, 12, 30});
    test_cal({année, 12, 31});
  }



  {
    msg("Test : DateHeure...");
    DateHeure t({{2022, 01, 25}, {7, 59, 30}});

    soit cal = t.calendrier();

    assertion_msg(
        (cal.année == 2022)
     && (cal.mois  == 1)
     && (cal.jour  == 25), "échec calendrier : {}", cal);

    msg("Test : DateHeure::decomposition()");
    soit hr = t.decomposition();
    assertion(
        (hr.jour == cal)
        && (hr.heure.heure == 7)
        && (hr.heure.minutes == 59)
        && (hr.heure.secondes == 30)
        && (hr.heure.ms == 0)
        && (hr.heure.µs == 0));


    t = DateHeure::maintenant();
    msg("Maintenant = {}", t);

    // Arrondi à la seconde près
    t -= Durée::microsecondes(t.microsecondes());

    soit [sem,secs] = t.vers_GPS();
    soit t2 = DateHeure::de_GPS(sem, secs);
    soit err = t - t2;
    msg("Erreur GPS avant / après = {}", err);
    assertion(err.nb_secondes() < 1e-6);

    {
      t = DateHeure::epoque_unix();
      msg("Epoque UNIX = {}", t);
      assertion((t - DateHeure({{1970,1,1}, {0,0,0}})).tics == 0);
      t = DateHeure::epoque_J2000();
      msg("Epoque J2000 = {}", t);
      assertion((t - DateHeure({{2000,1,1}, {12,0,0}})).tics == 0);
    }
  }
}
