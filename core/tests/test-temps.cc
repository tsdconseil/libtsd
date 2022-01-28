#include "tsd/tsd.hpp"
#include "tsd/tests.hpp"
#include <tuple>
#include "../include/tsd/temps.hpp"

using namespace tsd;
using namespace tsd::temps;


static void test_heure_composite()
{
  HeureComposite hc("23:10:02");
  tsd_assert((hc.heure == 23) && (hc.minutes == 10) && (hc.secondes == 2) && (hc.µs == 0));

  vérifie_exception([]() {HeureComposite hc("23/10:02");});
  vérifie_exception([]() {HeureComposite hc("25:10:02");});
  vérifie_exception([]() {HeureComposite hc("ab:10:02");});
}



static void test_cal(const Calendrier &cal)
{
  DateHeure t({cal, {7, 59, 30}});
  auto cal2 = t.calendrier();
  tsd_assert_msg(cal == cal2, "Erreur calendrier : {} VS {}", cal, cal2);
}

// A compléter !!!
int test_date_heure()
{
  test_heure_composite();

  {
    msg("Test calendrier...");
    Calendrier c(2022, 01, 01);
    tsd_assert(c.est_valide());
    tsd_assert(c.nb_jours_debut_année() == 0);
    c.jour = 10;
    tsd_assert(c.nb_jours_debut_année() == 9);
    c.année = 1;
    tsd_assert_msg(c.nb_jours_debut_ère() == 9, "nb = {}", c.nb_jours_debut_ère());
    c.année = 2;
    tsd_assert_msg(c.nb_jours_debut_ère() == 9 + 365, "nb = {}", c.nb_jours_debut_ère());
  }

  tsd_assert(est_bissextile(2020));
  tsd_assert(!est_bissextile(2021));
  tsd_assert(!est_bissextile(1900));
  tsd_assert(est_bissextile(2000));
  tsd_assert(!est_bissextile(2100));
  tsd_assert(!est_bissextile(2200));
  tsd_assert(!est_bissextile(2300));
  tsd_assert(est_bissextile(2400));

  // Tests conversion Calendrier <-> DateHeure
  for(auto année : {2000, 2020, 2021, 2022, 2023})
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

    auto cal = t.calendrier();

    tsd_assert_msg(
        (cal.année == 2022)
     && (cal.mois  == 1)
     && (cal.jour  == 25), "echec calendrier : {}", cal);

    msg("Test : DateHeure::decomposition()");
    auto hr = t.decomposition();
    tsd_assert(
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

    auto [sem,secs] = t.vers_GPS();
    auto t2 = DateHeure::de_GPS(sem, secs);
    auto err = t - t2;
    msg("Erreur GPS avant / après = {}", err);
    tsd_assert(err.nb_secondes() < 1e-6);

    {
      t = DateHeure::epoque_unix();
      msg("Epoque UNIX = {}", t);
      tsd_assert((t - DateHeure({{1970,1,1}, {0,0,0}})).tics == 0);
      t = DateHeure::epoque_J2000();
      msg("Epoque J2000 = {}", t);
      tsd_assert((t - DateHeure({{2000,1,1}, {12,0,0}})).tics == 0);
    }
  }

  return 0;
}
