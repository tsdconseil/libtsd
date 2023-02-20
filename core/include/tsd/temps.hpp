#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"

#include <cstdint>
#include <iosfwd>
#include <compare>
#include <tuple>
#include <fmt/core.h>

namespace tsd::temps {

/** @addtogroup tsd-temps
 *  @{
 */


/** @brief %Heure de la journée, décomposée en heures, minutes, etc. */
struct HeureComposite
{
  /** @brief Heure, entre 0 et 23 */
  entier heure = 0,
  /** @brief Minutes, entre 0 et 59 */
      minutes = 0,
  /** @brief Secondes, entre 0 et 59 */
      secondes = 0,
  /** @brief Milli-secondes, entre 0 et 999 */
      ms = 0,
  /** @brief Micro-secondes, entre 0 et 999 */
      µs = 0;

  HeureComposite(){}

  /** @brief Constructeur */
  HeureComposite(entier heure, entier minutes, entier secondes, entier ms = 0, entier µs = 0);

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
  HeureComposite(const std::string &s);

  /** @brief Vérifie si les heures, minutes et secondes sont dans un intervalle valide. */
  bouléen vérifie_validité() const;
};


/** @brief Intervalle temporel, en nombre de micro-secondes.
 *
 *  <h3>Intervalle temporel</h3>
 *
 *  Cette structure représente un intervalle de temps, exprimé en nombre
 *  de micro-secondes (champs tics).
 *
 */
struct Durée
{
  /** @brief Nombre de micro-secondes. */
  int64_t tics = 0;

  /** @brief Constructeur, d'après un nombre de micro-secondes. */
  Durée(int64_t tics_ = 0) : tics(tics_) {}

  /** @brief Constructeur, d'après un nombre de jours, heures, etc. */
  Durée(entier jours, entier heures, entier minutes, entier secondes, entier microsecondes = 0);

  /** @brief Constructeur, d'après une heure composite. */
  Durée(const HeureComposite &hc);

  /** @brief Constructeur statique, d'après un nombre de micro-secondes. */
  static Durée microsecondes(int64_t cnt);

  /** @brief Constructeur statique, d'après un nombre de milli-secondes. */
  static Durée millisecondes(int64_t cnt);

  /** @brief Constructeur statique, d'après un nombre de secondes. */
  static Durée secondes(double cnt);

  /** @brief Constructeur statique, d'après un nombre de minutes. */
  static Durée minutes(double cnt);

  /** @brief Constructeur statique, d'après un nombre d'heures. */
  static Durée heures(double cnt);

  /** @brief Constructeur statique, d'après un nombre de jours. */
  static Durée jours(double cnt);

  /** @brief Opérateur de comparaison. */
  std::strong_ordering operator<=>(const Durée&) const = default;


  /** @brief %Durée totale, exprimés en nombre fractionnaire de jours. */
  double nb_jours() const;

  /** @brief %Durée totale, exprimés en nombre fractionnaire d'heures. */
  double nb_heures() const;

  /** @brief %Durée totale, exprimés en nombre fractionnaire de minutes. */
  double nb_minutes() const;

  /** @brief %Durée totale, exprimés en nombre fractionnaire de secondes. */
  double nb_secondes() const;

  /** @brief %Durée totale, exprimés en nombre fractionnaire de millisecondes. */
  double nb_millisecondes() const;

  /** @brief %Durée totale, exprimés en nombre fractionnaire de microsecondes. */
  double nb_microsecondes() const;
};

/** @brief Affichage d'un intervalle de temps. */
extern std::ostream& operator<<(std::ostream& strm, const Durée& t);

/** @brief Somme de 2 intervalles de temps. */
extern Durée operator+(const Durée& ts1, const Durée& ts2);

/** @brief Différence entre 2 intervalles de temps. */
extern Durée operator-(const Durée& ts1, const Durée& ts2);



/** @brief %Calendrier (date décomposée en année, mois, jour, l'heure n'est pas spécifiée). */
struct Calendrier
{
  /** @brief Année. */
  entier année = 0,
  /** @brief Mois (entre 1 et 12). */
      mois = 0,
  /** @brief Jour (entre 1 et 31). */
      jour = 0;

  /** @brief Constructeur par défaut (00/00/0000) */
  Calendrier(){}

  /** @brief Constructeur */
  Calendrier(entier année, entier mois, entier jour);

  /** @brief Constructeur, à partir d'une chaîne de caractères de type "JJ/MM/AAAA". */
  Calendrier(const std::string &s);

  /** @brief Nombre de jours entiers depuis le début de l'année en cours (= 0 pour le premier Janvier). */
  entier nb_jours_debut_année() const;

  /** @brief Vérifie si la date est valide (mois entre 1 et 12, jour entre 1 et 31, etc.) */
  bouléen est_valide() const;

  /** @brief Calcul du nombre de jours depuis le 1/1/1, 0h00 (= 0 pour le 1/1/1) */
  entier nb_jours_debut_ère() const;

  /** @brief Opérateur de comparaison. */
  std::strong_ordering operator<=>(const Calendrier&) const = default;
};



/** @brief Spécification composite d'un jour et d'une heure */
struct DateComposite
{
  /** @brief Jour */
  Calendrier jour;

  /** @brief Heure dans la journée */
  HeureComposite heure;

  /** @brief Opérateur de comparaison. */
  std::strong_ordering operator<=>(const DateComposite&) const = default;
};

/** @brief Date et heure, avec fonctions de conversion vers différents formats (calendriers, temps sidéral, heure UTC, etc.). */
struct DateHeure
{
  /** @brief Nombre total de micro-secondes depuis le 1/1/0, 0h00. */
  int64_t ntics = 0;

  /** @brief Constructeur (d'après un nombre total de micro-secondes depuis le 1/1/0, 0h00). */
  DateHeure(int64_t tics = 0): ntics(tics){}

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
  DateHeure(const DateComposite &date);

  /** @brief Heure renvoyée par le système. */
  static DateHeure maintenant();

  /** @brief Constructeur, d'après l'année et un nombre fractionnaire de jours.
   *  @param année Année
   *  @param jour  Nombre fractionnaire de jours dans l'année (ex : 1.5 = 1er Janvier, 12h00).
   */
  DateHeure(entier année, double jour);

  /** @brief Calcul le temps sidéral de Greenwich.
   *
   *  <h3>Temps sidéral Greenwich</h3>
   *  Le temps sidéral de Greenwich (GST),
   *
   *  @return Temps sidéral, en radians (entre 0 et @f$2\pi@f$).
   */
  double temps_sidéral_Greenwich() const;

  /** @brief Calcul le temps sidéral local (exprimé en radians).
   *
   *  <h3>Temps sidéral local</h3>
   *  Le temps sidéral local est le temps sidéral de Greenwich, auquel on ajoute
   *  la longitude locale.
   *
   *  @param longitude Longitude locale (en radians)
   *  @return Temps sidéral local, en radians (entre 0 et @f$2\pi@f$).
   */
  double temps_sidéral_local(double longitude) const;

  /** @brief Nombre de "jours Juliens" (nombre de jours depuis le premier Janvier, 4713 BC, 12h00). */
  double nb_jours_Julien() const;

  /** @brief "Date Julienne modifiée" : nombre de jours depuis l'époque J2000 (à minuit au lieu de 12h00).
   *
   *  @sa nb_jours_Julien()
   */
  double J2000() const;


  /** @brief 1er janvier 1970 00:00:00 UTC */
  static DateHeure epoque_unix();

  /** @brief 1er janvier 2000 00:00:00 UTC */
  static DateHeure epoque_Matlab();

  /** @brief 1er janvier 2000 12:00:00 UTC */
  static DateHeure epoque_J2000();

  /** @brief 6 janvier 1980 00:00:00 UTC */
  static DateHeure epoque_GPS();

  /** @brief Décompositon de la date (année, mois, jour), en UTC. */
  Calendrier calendrier() const;

  /** @brief Décompositon de l'heure (heure, minutes, etc.), en UTC. */
  HeureComposite decomp_heure() const;

  /** @brief Décomposition date et heure, en UTC. */
  DateComposite decomposition() const;

  /** @brief Décomposition date et heure, suivant l'heure locale. */
  DateComposite decomposition_locale() const;


  /** @brief Calcul de l'heure GPS, en nombre de semaines, et nombre de secondes. */
  std::tuple<entier, entier> vers_GPS() const;

  /** @brief Constructeur, à partir de l'heure GPS. */
  static DateHeure de_GPS(entier semaine, entier secs);

  /** @brief Constructeur, à partir d'une chaîne de caractères de type "aaaa::mm::jj hh:mm:ss". */
  static DateHeure lis_chaine(const std::string &s);

  /** @brief Nombre de microsecondes depuis la dernière seconde entière, entre 0 et 1e6-1. */
  entier microsecondes() const;


  std::strong_ordering operator<=>(const DateHeure&) const = default;
};



/** @brief Conversion année + jour du calendrier Grégorien vers nombre de jours
 * depuis le 1° Janvier 0, 0h00.
 *  @param année : 2022, etc.
 *  @param jour_année : 1 pour le premier Janvier 0h00, 1.5 pour le premier Janvier 12h00, etc. */
extern double grégorien_vers_jours(entier année, double jour_année);




/** @brief Vérifique que le mois et l'année sont bien dans l'intervalle supporté.
 *  @param année  Numéro d'année (suivant l'ère usuelle)
 *  @param mois   Numéro de mois (intervalle de validité : 1 - 12)
 *  @returns      Vrai ou faux
 */
extern bouléen année_mois_valide(entier année, entier mois);

/** @brief Nombre de jours pour un mois et une année donnée.
 *
 * <h3>Nombre de jours pour un mois et une année donnée</h3>
 *
 * @param année Numéro d'année suivant l'ère usuelle.
 * @param mois Numéro de mois (entre 1 et 12).
 * @returns Nombre de jours, entre 28 et 31.
 */
extern entier mois_nb_jours(entier année, entier mois);

/** @brief Vérifie si l'année est bissextile ou non.
 *
 *  <h3>Test année bissextile</h3>
 *
 *  Cette fonction vérifie si l'année passée en paramètre est bissextile (c'est-à-dire une année à 366 jours au lieu de 365) ou non.
 *
 *  Les années bissextiles (pour le calendrier grégorien) sont celles divisibles par 4 et non par 100, et ainsi que celles divisibles par 400.
 *
 *  @param année Numéro d'année suivant l'ère usuelle.
 *  @returns Vrai ou faux
 *
 */
extern bouléen est_bissextile(entier année);

/** @brief Vérifie si l'année est dans l'intervalle supporté (entre 1 and 9999).
 *  @param   année Numéro d'année à vérifier (comptée suivant l'ère usuelle)
 *  @returns Vrai ou faux.
 */
extern bouléen année_est_valide(entier année);

/** @brief Affichage d'une date / heure. */
extern std::ostream& operator<<(std::ostream& strm, const DateHeure &t);

/** @brief Affichage d'une date. */
extern std::ostream& operator<<(std::ostream& strm, const Calendrier &date);

extern std::ostream& operator<<(std::ostream& strm, const HeureComposite &date);






/** @brief Différence entre deux points temporels */
inline Durée operator-(const DateHeure& dt1, const DateHeure& dt2)
{
  return Durée(dt1.ntics - dt2.ntics);
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
inline DateHeure operator+(const DateHeure& dt1, const Durée& dt2)
{
  return DateHeure(dt1.ntics + dt2.tics);
}

inline Durée operator/(const Durée& d, double r)
{
  return Durée(d.tics / r);
}

inline Durée operator*(const Durée& d, double r)
{
  return Durée(d.tics * r);
}

/** @brief Décrémente une heure d'une certaine durée. */
inline DateHeure operator-(const DateHeure& dt1, const Durée& dt2)
{
  return DateHeure(dt1.ntics - dt2.tics);
}

/** @brief Incrémente une heure d'une certaine durée. */
inline DateHeure operator+=(DateHeure& dt1, const Durée& dt2)
{
  dt1.ntics += dt2.tics;
  return dt1;
}

inline DateHeure operator-=(DateHeure& dt1, const Durée& dt2)
{
  dt1.ntics -= dt2.tics;
  return dt1;
}

/** @brief Si vrai (par défaut), l'affichage des dates/heures est fait suivant le standard UTC, sinon, suivant l'heure locale. */
extern bouléen mode_utc;

/** @} */

}


/*template <> struct fmt::formatter<tsd::temps::DateHeure> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::temps::DateHeure& t, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    tsd::temps::DateComposite dc;
    if(tsd::temps::mode_utc)
      dc = t.decomposition();
    else
      dc = t.decomposition_locale();
    return fmt::format_to(ctx.out(), "{:0>4d}-{:0>2d}-{:0>2d} {:0>2d}:{:0>2d}:{:0>2d} {}",
      dc.jour.année, dc.jour.mois, dc.jour.jour, dc.heure.heure, dc.heure.minutes, dc.heure.secondes, tsd::temps::mode_utc ? "(UTC)" : "(local hour)");
  }
};


template <> struct fmt::formatter<tsd::temps::Calendrier> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {return ctx.begin();}
  template <typename FormatContext>
  auto format(const tsd::temps::Calendrier& date, FormatContext& ctx) const -> decltype(ctx.out()) 
  {
    return fmt::format_to(ctx.out(), "{:0>4d}-{:0>2d}-{:0>2d}", date.année, date.mois, date.jour);
  }
};*/

ostream_formater(tsd::temps::Durée)
ostream_formater(tsd::temps::DateHeure)
ostream_formater(tsd::temps::Calendrier)
ostream_formater(tsd::temps::HeureComposite)

