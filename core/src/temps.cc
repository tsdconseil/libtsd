#include "tsd/tsd.hpp"
#include <fstream>
#include <iomanip>
#include <chrono>
#include "tsd/temps.hpp"

using namespace std::chrono;
using namespace std;

namespace tsd::temps {


static const int64_t
  TICS_US       = 1,
  TICS_MS       = 1000,
  TICS_SECONDE  = 1000 * TICS_MS,
  TICS_MINUTE   =   60 * TICS_SECONDE,
  TICS_HEURE    =   60 * TICS_MINUTE,
  TICS_JOUR     =   24 * TICS_HEURE;


static const int nbsecs_par_jour = 24 * 3600;

// Nombre de jours / cycle de 400 ans
// =  400 années
//  + 100 années bissextiles, sauf 3 années multiples de 100 et pas de 400.
static const int nb_jours_cycle_400_ans = 365 * 400 + 100 - 3;


// Nombre de jours / cycle de 100 ans
// (24 années bisextiles).
static const int nb_jours_cycle_100_ans = 365 * 100 + 24;


// Nombre de jours / cycle de 4 ans
// (1 année bisextile).
static const int nb_jours_cycle_4_ans = 365 * 4 + 1;

static const int nb_jours_par_mois[2][12] = {
 //  1   2   3   4   5   6   7   8   9   10  11  12
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

static const int cumul_jours_par_mois[2][12] = {
    //  1  2   3   4   5    6    7    8    9    10   11   12
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
};

Durée operator+(const Durée& ts1, const Durée& ts2)
{
  return Durée(ts1.tics + ts2.tics);
}

Durée operator-(const Durée& ts1, const Durée& ts2)
{
  return Durée(ts1.tics - ts2.tics);
}

double Durée::nb_jours() const
{
    return ((double) tics) / TICS_JOUR;
}

double Durée::nb_heures() const
{
    return ((double) tics) / TICS_HEURE;
}

double Durée::nb_minutes() const
{
    return ((double) tics) / TICS_MINUTE;
}

double Durée::nb_secondes() const
{
    return ((double) tics) / TICS_SECONDE;
}

double Durée::nb_millisecondes() const
{
    return ((double) tics) / TICS_MS;
}

double Durée::nb_microsecondes() const
{
    return ((double) tics) / TICS_US;
}


std::ostream& operator<<(std::ostream& ss, const Durée& t)
{
  ss << t.nb_microsecondes() << " µs";
  return ss;
}

Durée Durée::microsecondes(int64_t cnt)
{
  return Durée(cnt);
}

Durée Durée::millisecondes(int64_t cnt)
{
  return Durée(cnt * TICS_MS);
}

Durée Durée::secondes(double cnt)
{
  return Durée(cnt * TICS_SECONDE);
}

Durée Durée::minutes(double cnt)
{
  return Durée(cnt * TICS_MINUTE);
}

Durée Durée::heures(double cnt)
{
  return Durée(cnt * TICS_HEURE);
}
Durée Durée::jours(double cnt)
{
  return Durée(cnt * TICS_JOUR);
}

Durée::Durée(const HeureComposite &hc)
{
  tics = (hc.heure * 3600 + hc.minutes * 60 + hc.secondes) * TICS_SECONDE +
      hc.ms * TICS_MS  + hc.µs * TICS_US;
}

Durée::Durée(
    int jours,
    int heures,
    int minutes,
    int secondes,
    int microsecondes)
{
  tics = jours * TICS_JOUR +
      (heures * 3600 + minutes * 60 + secondes) * TICS_SECONDE
      + microsecondes * TICS_US;
}



DateHeure::DateHeure(int année, double jours)
{
  ntics = Durée::jours(grégorien_vers_jours(année, jours)).tics;
}

double DateHeure::temps_sidéral_local(double longitude) const
{
  return modulo_2π(temps_sidéral_Greenwich() + longitude);
}

DateHeure::DateHeure(const DateComposite &date)
{
  ntics = (Durée(date.heure) + Durée::jours(date.jour.nb_jours_debut_ère())).tics;
}

int DateHeure::microsecondes() const
{
  return (int) (ntics % TICS_SECONDE / TICS_US);
}

bool année_mois_valide(int année, int mois)
{
  if(!année_est_valide(année))
    return false;
  return (mois >= 1) && (mois <= 12);
}

int mois_nb_jours(int année, int mois)
{
  if(!année_mois_valide(année, mois))
    echec("mois_nb_jours: année ({}) ou mois ({}) invalide.", année, mois);

  auto ptr = est_bissextile(année) ? nb_jours_par_mois[1] : nb_jours_par_mois[0];
  return ptr[mois-1];
}

bool Calendrier::est_valide() const
{
  return année_mois_valide(année, mois)
      &&  (jour >= 1) && (jour <= mois_nb_jours(année, mois));
}


static std::vector<std::string> split(const std::string& str, const std::string& delim)
{
  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do
  {
    pos = str.find(delim, prev);
    if(pos == std::string::npos)
      pos = str.length();
    std::string token = str.substr(prev, pos-prev);
    if(!token.empty())
      tokens.push_back(token);
    prev = pos + delim.length();
  }
  while (pos < str.length() && prev < str.length());
  return tokens;
}

Calendrier::Calendrier(int année, int mois, int jour)
{
  this->année = année;
  this->mois  = mois;
  this->jour  = jour;
}

Calendrier::Calendrier(const std::string &s)
{
  auto sp = split(s, "/");
  if(sp.size() != 3)
    echec("Date invalide : {}", s);
  else
  {
    jour  = stoi(sp[0]);
    mois  = stoi(sp[1]);
    année = stoi(sp[2]);
  }
  if(!est_valide())
    echec("Date invalide : {}", s);
}

HeureComposite::HeureComposite(int heure, int minutes, int secondes, int ms, int µs)
{
  this->heure     = heure;
  this->minutes   = minutes;
  this->secondes  = secondes;
  this->ms        = ms;
  this->µs        = µs;
}

HeureComposite::HeureComposite(const std::string &s)
{
  auto sp = split(s, ":");
  if(sp.size() != 3)
    echec("HeureComposite : chaine invalide : {}", s);
  else
  {
    heure    = stoi(sp[0]);
    minutes  = stoi(sp[1]);
    secondes = stoi(sp[2]);
  }
  if(!vérifie_validité())
    echec("Heure composite invalide.");
}


bool HeureComposite::vérifie_validité() const
{
  return (heure >= 0) && (heure < 24)
      && (minutes >= 0) && (minutes < 60)
      && (secondes >= 0) && (secondes < 60)
      && (ms >= 0) && (ms < 1e3)
      && (µs >= 0) && (µs < 1e3);
}


int Calendrier::nb_jours_debut_année() const
{
  if(!est_valide())
  {
    msg_avert("Date de calendrier invalide ({}/{}/{})", année, mois, jour);
    return jour;
  }

  return (jour - 1) + cumul_jours_par_mois[est_bissextile(année) ? 1 : 0][mois-1];
}

double grégorien_vers_jours(int année, double jour_année)
{
  int64_t a1 = année - 1;

  return
  // Nombre de jours de l'année précédente :
    365 * a1
  // Nombre de jours supplémentaires dues aux années bissextiles :
    + a1 / 4
  // Mais les années multiples de 100 ne sont pas bissextiles :
    - a1 / 100
  // Sauf les années multiples de 400 :
    + a1 / 400
  // Nombres de jours (moins 1, car 1° Janvier)
    + jour_année - 1;
}

int Calendrier::nb_jours_debut_ère() const
{
  return grégorien_vers_jours(année, 1 + nb_jours_debut_année());
}


// Pré-calcul
static const DateHeure cst_epoque_unix = DateHeure({{1970,1,1}, {0, 0, 0}});

DateHeure DateHeure::epoque_unix()
{
  return cst_epoque_unix;
}

DateHeure DateHeure::maintenant()
{
  return epoque_unix() + Durée::microsecondes(
      duration_cast<microseconds>(system_clock::now().time_since_epoch()).count());
}

double DateHeure::nb_jours_Julien() const
{
  // A comparer avec :     % reference: USNO http://aa.usno.navy.mil/faq/docs/JD_Formula.php
  // JD = 367.0 * yr  ...
  //         - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
  //         + floor( 275 * mon / 9.0 ) ...
  //         + jour + 1721013.5  ...
  //         + ( (sec/60.0 + mn ) / 60.0 + hr ) / 24.0;
  return Durée(ntics).nb_jours() + 1721425.5;
}

double DateHeure::J2000() const
{
  return nb_jours_Julien() - 2415020.0;
}

DateHeure DateHeure::epoque_Matlab()
{
  // 01/01/2000, 00h00
  return DateHeure{{{2000, 1, 1}, {0, 0, 0}}};
}

DateHeure DateHeure::epoque_J2000()
{
  // 01/01/2000, 12h00
  return DateHeure{{{2000, 1, 1}, {12, 0, 0}}};
}

DateHeure DateHeure::epoque_GPS()
{
  // 06/01/1980, 00h00
  return DateHeure({{1980, 1, 6}, {0, 0, 0}});
}

double DateHeure::temps_sidéral_Greenwich() const
{
  // https://lweb.cfa.harvard.edu/~jzhao/times.html

  // Nb jours Juliens du minuit précédent
  auto jd0 = floor(nb_jours_Julien() + 0.5) - 0.5;
  auto t   = (jd0 - 2451545.0) / 36525.0;
  auto jdf = nb_jours_Julien() - jd0;
  auto gt  = 24110.54841 + t * (8640184.812866 + t * (0.093104 - t * 6.2E-6));
  gt  += jdf * 1.00273790935 * 86400.0;
  return modulo_2π(deg2rad(gt * 360.0 / nbsecs_par_jour));
}


bool est_bissextile(int année)
{
  if(!année_est_valide(année))
    return false;
  return (((année % 4) == 0 && (année % 100) != 0) || (année % 400) == 0);
}

bool année_est_valide(int année)
{
  return (année >= 1) && (année <= 9999);
}


HeureComposite DateHeure::decomp_heure() const
{
  return
      {(int) ((ntics % TICS_JOUR)     / TICS_HEURE),
       (int) ((ntics % TICS_HEURE)    / TICS_MINUTE),
       (int) ((ntics % TICS_MINUTE)   / TICS_SECONDE),
       (int) ((ntics % TICS_SECONDE)  / TICS_MS),
       (int) ((ntics % TICS_MS)       / TICS_US)};
}

DateComposite DateHeure::decomposition() const
{
  return {calendrier(), decomp_heure()};
}

DateHeure DateHeure::lis_chaine(const std::string &s_)
{
  int y,M,d,h,m,s;
  sscanf(s_.c_str(), "%d-%d-%d %d:%d:%d", &y,&M,&d,&h,&m,&s);
  return DateHeure{{{y,M,d},{h,m,s}}};
}

DateHeure DateHeure::de_GPS(int semaine, int secs)
{
  // 5 Janvier 1980
  auto r = epoque_GPS();
  r += Durée::secondes(secs);
  r += Durée::jours(7 * semaine);
  r += Durée::secondes(-18);
  return r;
}

std::tuple<int, int> DateHeure::vers_GPS() const
{
  // 5 Janvier 1980
  auto e0 = epoque_GPS();

  auto nsecs = (*this - e0).nb_secondes();

  // 06/2017 : 18 leap seconds ahead of UTC
  nsecs += 18;

  int sps         = nbsecs_par_jour * 7;
  int nb_semaines = nsecs / sps;
  int nb_secs     = nsecs - nb_semaines * sps;

  return {nb_semaines, nb_secs};
}

DateComposite DateHeure::decomposition_locale() const
{
  // timezone = secondes à l'ouest de Greenwich
  // tzset();
  //msg("timezone = {}", timezone);

  /*static int decalage = -1000;

  if(decalage == -1000)
  {
    std::time_t t = std::time(nullptr);
    struct tm t1, t2;
    t1 = *(std::gmtime(&t)); // UTC
    t2 = *(std::localtime(&t));

    t2.tm_isdst = 1;

    decalage = difftime(mktime(&t1), mktime(&t2));

    msg("\n\n\n*** Décalage horaire locale vs UTC : {}\n\n\n", decalage);
  }*/

  //auto s = std::chrono::sys_info.offset;
  //int decalage = timezone;
  //decalage = 0;
  //using namespace std::chrono;
  //static int decalage = std::chrono::get_tzdb().current_zone()->get_info(system_clock::now()).offset;
  //return ajoute_secondes(-decalage).decomposition();

  // TODO!!!
  //return (*this + Durée::heures(2)).decomposition();
  return (*this + Durée::heures(1)).decomposition();
}

Calendrier DateHeure::calendrier() const
{
  int nb_jours = ntics / TICS_JOUR;

  // Décompte les cycles de 400 ans
  int ncycles_400_ans = nb_jours / nb_jours_cycle_400_ans;
  nb_jours %= nb_jours_cycle_400_ans;

  // Décompte les cycles de 100 ans
  // 1,2,3,...99,  100,101,...,199, ....
  int ncycles_100_ans = nb_jours / nb_jours_cycle_100_ans;
  if (ncycles_100_ans == 4)
      // 31 décembre de la dernière année bissextile
      ncycles_100_ans = 3;
  nb_jours -= ncycles_100_ans * nb_jours_cycle_100_ans;

  // Décompte les cycles de 4 ans
  int ncycles_4_ans = nb_jours / nb_jours_cycle_4_ans;
  nb_jours %= nb_jours_cycle_4_ans;

  // Décompte les années qui restent (dans le cycle de 4 ans en cours)
  // Un cycle de 4 ans =
  // 1, 2, 3, 4
  // ...
  // 2021, 2022, 2023, 2024
  // Soit en général :
  // [365, 365, 365, 366]
  // 0 <= nb jours < 3 * 365 + 366
  int nb_ans = 0;
  if(nb_jours < 365)
    ;
  else if(nb_jours < 2 * 365)
  {
    nb_ans = 1;
    nb_jours -= 365;
  }
  else if(nb_jours < 3 * 365)
  {
    nb_ans = 2;
    nb_jours -= 2 * 365;
  }
  else
  {
    nb_ans = 3;
    nb_jours -= 3 * 365;
  }

  int année = (ncycles_400_ans * 400) + (ncycles_100_ans * 100)
            + (ncycles_4_ans * 4) + nb_ans + 1;

  // Nombre de jours pour chaque mois
  auto ptr = est_bissextile(année) ? nb_jours_par_mois[1] : nb_jours_par_mois[0];

  int mois = 0;
  while((nb_jours >= ptr[mois]) && (mois <= 11))
    nb_jours -= ptr[mois++];

  return {année, mois + 1, nb_jours + 1};
}

bool mode_utc = false;

std::ostream& operator<<(std::ostream& ss, const Calendrier &date)
{
  ss << fmt::format("{:0>4d}-{:0>2d}-{:0>2d}", date.année, date.mois, date.jour);
  return ss;
}

std::ostream& operator<<(std::ostream& ss, const DateHeure &t)
{
  DateComposite dc;
  if(mode_utc)
    dc = t.decomposition();
  else
    dc = t.decomposition_locale();
  ss << fmt::format("{:0>4d}-{:0>2d}-{:0>2d} {:0>2d}:{:0>2d}:{:0>2d}",
      dc.jour.année, dc.jour.mois, dc.jour.jour, dc.heure.heure, dc.heure.minutes, dc.heure.secondes);
  if(mode_utc)
    ss << " (UTC)";
  else
    ss << " (local hour)";
  return ss;
}
}
