#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/temps.hpp"
#include <cstring>

namespace dsp::time {

/** @addtogroup tsd-temps
 *  @{
 */


  namespace tsdt = tsd::temps;

/** @brief %Heure de la journée, décomposée en heures, minutes, etc. */
struct HourComposite
{
  /** @brief Heure, entre 0 et 23 */
  int hour = 0,
  /** @brief Minutes, entre 0 et 59 */
      minutes = 0,
  /** @brief Secondes, entre 0 et 59 */
      seconds = 0,
  /** @brief Milli-secondes, entre 0 et 999 */
      ms = 0,
  /** @brief Micro-secondes, entre 0 et 999 */
      µs = 0;

  HourComposite(){}

  HourComposite(const tsd::temps::HeureComposite &hc)
  {
    memcpy(this, &hc, sizeof(HourComposite));
  }

  auto hc() const
  {
    tsd::temps::HeureComposite res;
    memcpy((void *) &res, this, sizeof(HourComposite));
    return res;
  }

  /** @brief Constructeur */
  HourComposite(int hour, int minutes, int seconds, int ms = 0, int µs = 0)
  {
    *this = tsdt::HeureComposite(hour, minutes, seconds, ms, µs);
  }

  /** @brief Constructeur, à partir d'une chaine de caractères de type HH::MM::SS
   *
   *  <h3>Constructeur, à partir d'une chaine de caractères de type HH::MM::SS</h3>
   *
   *  Lève une exception si la chaîne est mal formatée.
   *
   *  @par Exemple
   *  @code
   *  HeureComposite hc("23:10:02");
   *  tsd_assert((hc.heure == 23) && (hc.minutes == 10) && (hc.secondes == 2) && (hc.ms == 0) && (hc.µs == 0));
   *  @endcode
   *
   */
  HourComposite(const std::string &s)
  {
    *this = tsdt::HeureComposite(s);
  }

  /** @brief Vérifie si les heures, minutes et secondes sont dans un intervalle valide. */
  bool check_valid() const
  {
    return hc().vérifie_validité();
  }
};


/** @brief Intervalle temporel, en nombre de micro-secondes.
 *
 *  <h3>Intervalle temporel</h3>
 *
 *  Cette structure représente un intervalle de temps, exprimé en nombre
 *  de micro-secondes (champs tics).
 *
 */
struct Duration
{
  /** @brief Nombre de micro-secondes. */
  int64_t tics = 0;

  /** @brief Constructeur, d'après un nombre de micro-secondes. */
  Duration(int64_t tics_ = 0) : tics(tics_) {}

  Duration(const tsdt::Durée &d)
  {
    tics = d.tics;
  }

  tsdt::Durée dr() const
  {
    return tsdt::Durée(tics);
  }

  /** @brief Constructeur, d'après un nombre de jours, heures, etc. */
  Duration(int days, int hours, int minutes, int seconds, int microseconds = 0)
  : Duration(tsdt::Durée(days, hours, minutes, seconds, microseconds))
  {
  }

  /** @brief Constructeur, d'après une heure composite. */
  Duration(const HourComposite &hc)
  : Duration(tsdt::Durée(hc.hc()))
  {

  }

  /** @brief Constructeur statique, d'après un nombre de micro-secondes. */
  static Duration microseconds(int64_t cnt)
  {
    return tsdt::Durée::microsecondes(cnt);
  }

  /** @brief Constructeur statique, d'après un nombre de milli-secondes. */
  static Duration milliseconds(int64_t cnt)
  {
    return tsdt::Durée::millisecondes(cnt);
  }


  /** @brief Constructeur statique, d'après un nombre de secondes. */
  static Duration seconds(double cnt)
  {
    return tsdt::Durée::secondes(cnt);
  }

  /** @brief Constructeur statique, d'après un nombre de minutes. */
  static Duration minutes(double cnt)
  {
    return tsdt::Durée::minutes(cnt);
  }

  /** @brief Constructeur statique, d'après un nombre d'heures. */
  static Duration hours(double cnt)
  {
    return tsdt::Durée::heures(cnt);
  }

  /** @brief Constructeur statique, d'après un nombre de jours. */
  static Duration days(double cnt)
  {
    return tsdt::Durée::jours(cnt);
  }

  /** @brief Opérateur de comparaison. */
  std::strong_ordering operator<=>(const Duration&) const = default;


  /** @brief %Durée totale, exprimés en nombre fractionnaire de jours. */
  double nb_days() const
  {
    return dr().nb_jours();
  }

  /** @brief %Durée totale, exprimés en nombre fractionnaire d'heures. */
  double nb_hours() const
  {
    return dr().nb_heures();
  }

  /** @brief %Durée totale, exprimés en nombre fractionnaire de minutes. */
  double nb_minutes() const
  {
    return dr().nb_minutes();
  }

  /** @brief %Durée totale, exprimés en nombre fractionnaire de secondes. */
  double nb_seconds() const
  {
    return dr().nb_secondes();
  }

  /** @brief %Durée totale, exprimés en nombre fractionnaire de millisecondes. */
  double nb_milliseconds() const
  {
    return dr().nb_millisecondes();
  }

  /** @brief %Durée totale, exprimés en nombre fractionnaire de microsecondes. */
  double nb_microseconds() const
  {
    return dr().nb_microsecondes();
  }
};

/** @brief Affichage d'un intervalle de temps. */
inline std::ostream& operator<<(std::ostream& strm, const Duration& t)
{
  return strm << t.dr();
}

/** @brief Somme de 2 intervalles de temps. */
inline Duration operator+(const Duration& ts1, const Duration& ts2)
{
  return ts1.dr() + ts2.dr();
}

/** @brief Différence entre 2 intervalles de temps. */
inline Duration operator-(const Duration& ts1, const Duration& ts2)
{
  return ts1.dr() - ts2.dr();
}



/** @brief %Calendrier (date décomposée en année, mois, jour, l'heure n'est pas spécifiée). */
struct Calendar
{
  /** @brief Année. */
  int year = 0,
  /** @brief Mois (entre 1 et 12). */
      month = 0,
  /** @brief Jour (entre 1 et 31). */
      day = 0;

  Calendar(const tsdt::Calendrier &c)
  {
    memcpy(this, &c, sizeof(Calendar));
  }

  auto cl() const
  {
    tsdt::Calendrier res;
    memcpy((void *) &res, this, sizeof(Calendar));
    return res;
  }

  /** @brief Constructeur par défaut (00/00/0000) */
  Calendar(){}

  /** @brief Constructeur */
  Calendar(int year, int month, int day)
  {
    *this = tsdt::Calendrier(year, month, day);
  }

  /** @brief Constructeur, à partir d'une chaîne de caractères de type "JJ/MM/AAAA". */
  Calendar(const std::string &s)
  {
    *this = tsdt::Calendrier(s);
  }

  /** @brief Nombre de jours entiers depuis le début de l'année en cours (= 0 pour le premier Janvier). */
  int nb_days_since_beginning_of_year() const
  {
    return cl().nb_jours_debut_année();
  }

  /** @brief Vérifie si la date est valide (mois entre 1 et 12, jour entre 1 et 31, etc.) */
  bool is_valid() const
  {
    return cl().est_valide();
  }

  /** @brief Calcul du nombre de jours depuis le 1/1/1, 0h00 (= 0 pour le 1/1/1) */
  int nb_days_since_beginning_of_era() const
  {
    return cl().nb_jours_debut_ère();
  }

  /** @brief Opérateur de comparaison. */
  std::strong_ordering operator<=>(const Calendar&) const = default;
};



/** @brief Spécification composite d'un jour et d'une heure */
struct DateComposite
{
  /** @brief Jour */
  Calendar day;

  /** @brief Heure dans la journée */
  HourComposite hour;

  DateComposite(const Calendar &day, const HourComposite &hour)
  {
    this->day  = day;
    this->hour = hour;
  }

  DateComposite(const tsdt::DateComposite &dc)
  {
    day   = dc.jour;
    hour = dc.heure;
  }

  tsdt::DateComposite dc() const
  {
    return tsdt::DateComposite{day.cl(), hour.hc()};
  }

  /** @brief Opérateur de comparaison. */
  std::strong_ordering operator<=>(const DateComposite&) const = default;
};

/** @brief Date et heure, avec fonctions de conversion vers différents formats (calendriers, temps sidéral, heure UTC, etc.). */
struct DateTime
{
  /** @brief Nombre total de micro-secondes depuis le 1/1/0, 0h00. */
  int64_t ntics = 0;

  DateTime(const tsdt::DateHeure &hc)
  {
    memcpy(this, &hc, sizeof(DateTime));
  }

  auto dt() const
  {
    tsdt::DateHeure res;
    memcpy((void *) &res, this, sizeof(DateTime));
    return res;
  }

  /** @brief Constructeur (d'après un nombre total de micro-secondes depuis le 1/1/0, 0h00). */
  DateTime(int64_t tics = 0): ntics(tics){}

  /** @brief Conctructeur (d'après calendrier Grégorien).
   *
   *  <h3>Conctructeur (d'après calendrier Grégorien)</h3>
   *
   * @par Exemple
   * @code
   * // 8 Février 2021, 12h00
   * DateHeure t {{2021, 2, 8}, {12, 00, 00}}
   * @endcode
   *  */
  DateTime(const DateComposite &date)
  {
    *this = tsdt::DateHeure(date.dc());
  }

  /** @brief Heure renvoyée par le système. */
  static DateTime now()
  {
    return tsdt::DateHeure::maintenant();
  }

  /** @brief Constructeur, d'après l'année et un nombre fractionnaire de jours.
   *  @param year Année
   *  @param day  Nombre fractionnaire de jours dans l'année (ex : 1.5 = 1er Janvier, 12h00).
   */
  DateTime(int year, double day)
  {
    *this = tsdt::DateHeure(year, day);
  }

  /** @brief Calcul le temps sidéral de Greenwich.
   *
   *  <h3>Temps sidéral Greenwich</h3>
   *  Le temps sidéral de Greenwich (GST),
   *
   *  @return Temps sidéral, en radians (entre 0 et @f$2\pi@f$).
   */
  double Greenwich_sidereal_time() const
  {
    return dt().temps_sidéral_Greenwich();
  }

  /** @brief Calcul le temps sidéral local (exprimé en radians).
   *
   *  <h3>Temps sidéral local</h3>
   *  Le temps sidéral local est le temps sidéral de Greenwich, auquel on ajoute
   *  la longitude locale.
   *
   *  @param longitude Longitude locale (en radians)
   *  @return Temps sidéral local, en radians (entre 0 et @f$2\pi@f$).
   */
  double local_sidereal_time(double longitude) const
  {
    return dt().temps_sidéral_local(longitude);
  }

  /** @brief Nombre de "jours Juliens" (nombre de jours depuis le premier Janvier, 4713 BC, 12h00). */
  double nb_Julian_days() const
  {
    return dt().nb_jours_Julien();
  }

  /** @brief "Date Julienne modifiée" : nombre de jours depuis l'époque J2000 (à minuit au lieu de 12h00).
   *
   *  @sa nb_jours_Julien()
   */
  double J2000() const
  {
    return dt().J2000();
  }


  /** @brief 1er janvier 1970 00:00:00 UTC */
  static DateTime epoch_unix()
  {
    return tsdt::DateHeure::epoque_unix();
  }

  /** @brief 1er janvier 2000 00:00:00 UTC */
  static DateTime epoch_Matlab()
  {
    return tsdt::DateHeure::epoque_Matlab();
  }

  /** @brief 1er janvier 2000 12:00:00 UTC */
  static DateTime epoch_J2000()
  {
    return tsdt::DateHeure::epoque_J2000();
  }

  /** @brief 6 janvier 1980 00:00:00 UTC */
  static DateTime epoch_GPS()
  {
    return tsdt::DateHeure::epoque_GPS();
  }

  /** @brief Décompositon de la date (année, mois, jour), en UTC. */
  Calendar calendar() const
  {
    return dt().calendrier();
  }

  /** @brief Décompositon de l'heure (heure, minutes, etc.), en UTC. */
  HourComposite decomp_time() const
  {
    return dt().decomp_heure();
  }

  /** @brief Décomposition date et heure, en UTC. */
  DateComposite decomposition() const
  {
    return dt().decomposition();
  }

  /** @brief Décomposition date et heure, suivant l'heure locale. */
  DateComposite decomposition_local() const
  {
    return dt().decomposition_locale();
  }


  /** @brief Calcul de l'heure GPS, en nombre de semaines, et nombre de secondes. */
  std::tuple<int, int> to_GPS() const
  {
    return dt().vers_GPS();
  }

  /** @brief Constructeur, à partir de l'heure GPS. */
  static DateTime from_GPS(int semaine, int secs)
  {
    return tsdt::DateHeure::de_GPS(semaine, secs);
  }

  /** @brief Constructeur, à partir d'une chaîne de caractères de type "aaaa::mm::jj hh:mm:ss". */
  static DateTime parse(const std::string &s)
  {
    return tsdt::DateHeure::lis_chaine(s);
  }

  /** @brief Nombre de microsecondes depuis la dernière seconde entière, entre 0 et 1e6-1. */
  int microseconds() const
  {
    return dt().microsecondes();
  }


  std::strong_ordering operator<=>(const DateTime&) const = default;
};



/** @brief Conversion année + jour du calendrier Grégorien vers nombre de jours
 * depuis le 1° Janvier 0, 0h00.
 *  @param année : 2022, etc.
 *  @param jour_année : 1 pour le premier Janvier 0h00, 1.5 pour le premier Janvier 12h00, etc. */
inline double gregorian2days(int année, double jour_année)
{
  return tsdt::grégorien_vers_jours(année, jour_année);
}





/** @brief Nombre de jours pour un mois et une année donnée.
 *
 * <h3>Nombre de jours pour un mois et une année donnée</h3>
 *
 * @param year Numéro d'année suivant l'ère usuelle.
 * @param month Numéro de mois (entre 1 et 12).
 * @returns Nombre de jours, entre 28 et 31.
 */
inline int month_nb_days(int year, int month)
{
  return tsdt::mois_nb_jours(year, month);
}

/** @brief Vérifie si l'année est bissextile ou non.
 *
 *  <h3>Test année bissextile</h3>
 *
 *  Cette fonction vérifie si l'année passée en paramètre est bissextile (c'est-à-dire une année à 366 jours au lieu de 365) ou non.
 *
 *  Les années bissextiles (pour le calendrier grégorien) sont celles divisibles par 4 et non par 100, et ainsi que celles divisibles par 400.
 *
 *  @param year Numéro d'année suivant l'ère usuelle.
 *  @returns Vrai ou faux
 *
 */
inline bool is_bissextil(int year)
{
  return tsdt::est_bissextile(year);
}

/** @brief Vérifie si l'année est dans l'intervalle supporté (entre 1 and 9999).
 *  @param   year Numéro d'année à vérifier (comptée suivant l'ère usuelle)
 *  @returns Vrai ou faux.
 */
inline bool année_est_valide(int year)
{
  return tsdt::année_est_valide(year);
}

/** @brief Affichage d'une date / heure. */
inline std::ostream& operator<<(std::ostream& strm, const DateTime &t)
{
  return strm << t.dt();
}

/** @brief Affichage d'une date. */
inline std::ostream& operator<<(std::ostream& strm, const Calendar &date)
{
  return strm << date.cl();
}

inline std::ostream& operator<<(std::ostream& strm, const HourComposite &date)
{
  return strm << date.hc();
}

/** @brief Différence entre deux points temporels */
inline Duration operator-(const DateTime& dt1, const DateTime& dt2)
{
  return Duration(dt1.ntics - dt2.ntics);
}

/** @brief Ajout à un instant donné d'une certaine durée
 *
 *  <h3>Ajout à un instant donné d'une certaine durée</h3>
 *
 *  @par Example
 *  @code
 *  auto t0 = DateHeure::maintenant();
 *  auto t1 = t0 + IntervalleTemps::jours(3);
 *  @endcode
 */
inline DateTime operator+(const DateTime &dt1, const Duration &dt2)
{
  return DateTime(dt1.ntics + dt2.tics);
}

inline Duration operator/(const Duration &d, double r)
{
  return Duration(d.tics / r);
}

inline Duration operator*(const Duration &d, double r)
{
  return Duration(d.tics * r);
}

/** @brief Décrémente une heure d'une certaine durée. */
inline DateTime operator-(const DateTime &dt1, const Duration &dt2)
{
  return DateTime(dt1.ntics - dt2.tics);
}

/** @brief Incrémente une heure d'une certaine durée. */
inline DateTime operator+=(DateTime &dt1, const Duration &dt2)
{
  dt1.ntics += dt2.tics;
  return dt1;
}

inline DateTime operator-=(DateTime &dt1, const Duration &dt2)
{
  dt1.ntics -= dt2.tics;
  return dt1;
}

/** @brief Si vrai (par défaut), l'affichage des dates/heures est fait suivant le standard UTC, sinon, suivant l'heure locale. */
inline bool &mode_utc = tsdt::mode_utc;

/** @} */

}


ostream_formater(dsp::time::Duration)
ostream_formater(dsp::time::Calendar)
ostream_formater(dsp::time::DateTime)
ostream_formater(dsp::time::HourComposite)

