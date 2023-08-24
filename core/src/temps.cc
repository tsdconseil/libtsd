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


static const entier nbsecs_par_jour = 24 * 3600;

// Nombre de jours / cycle de 400 ans
// =  400 années
//  + 100 années bissextiles, sauf 3 années multiples de 100 et pas de 400.
static const entier nb_jours_cycle_400_ans = 365 * 400 + 100 - 3;


// Nombre de jours / cycle de 100 ans
// (24 années bisextiles).
static const entier nb_jours_cycle_100_ans = 365 * 100 + 24;


// Nombre de jours / cycle de 4 ans
// (1 année bisextile).
static const entier nb_jours_cycle_4_ans = 365 * 4 + 1;

static const entier nb_jours_par_mois[2][12] = {
 //  1   2   3   4   5   6   7   8   9   10  11  12
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

static const entier cumul_jours_par_mois[2][12] = {
    //  1  2   3   4   5    6    7    8    9    10   11   12
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
};

Durée operator+(const Durée& ts1, const Durée& ts2)
{
  retourne Durée(ts1.tics + ts2.tics);
}

Durée operator-(const Durée& ts1, const Durée& ts2)
{
  retourne Durée(ts1.tics - ts2.tics);
}

double Durée::nb_jours() const
{
    retourne ((double) tics) / TICS_JOUR;
}

double Durée::nb_heures() const
{
    retourne ((double) tics) / TICS_HEURE;
}

double Durée::nb_minutes() const
{
    retourne ((double) tics) / TICS_MINUTE;
}

double Durée::nb_secondes() const
{
    retourne ((double) tics) / TICS_SECONDE;
}

double Durée::nb_millisecondes() const
{
    retourne ((double) tics) / TICS_MS;
}

double Durée::nb_microsecondes() const
{
    retourne ((double) tics) / TICS_US;
}


std::ostream& operator<<(std::ostream& ss, const Durée& t)
{
  ss << t.nb_microsecondes() << " µs";
  retourne ss;
}

Durée Durée::microsecondes(int64_t cnt)
{
  retourne Durée(cnt);
}

Durée Durée::millisecondes(int64_t cnt)
{
  retourne Durée(cnt * TICS_MS);
}

Durée Durée::secondes(double cnt)
{
  retourne Durée(cnt * TICS_SECONDE);
}

Durée Durée::minutes(double cnt)
{
  retourne Durée(cnt * TICS_MINUTE);
}

Durée Durée::heures(double cnt)
{
  retourne Durée(cnt * TICS_HEURE);
}
Durée Durée::jours(double cnt)
{
  retourne Durée(cnt * TICS_JOUR);
}

Durée::Durée(const HeureComposite &hc)
{
  tics = (hc.heure * 3600 + hc.minutes * 60 + hc.secondes) * TICS_SECONDE +
      hc.ms * TICS_MS  + hc.µs * TICS_US;
}

Durée::Durée(
    entier jours, entier heures, entier minutes,
    entier secondes, entier microsecondes)
{
  tics = jours * TICS_JOUR +
      (heures * 3600 + minutes * 60 + secondes) * TICS_SECONDE
      + microsecondes * TICS_US;
}



DateHeure::DateHeure(entier année, double jours)
{
  ntics = Durée::jours(grégorien_vers_jours(année, jours)).tics;
}

double DateHeure::temps_sidéral_local(double longitude) const
{
  retourne modulo_2π(temps_sidéral_Greenwich() + longitude);
}

DateHeure::DateHeure(const DateComposite &date)
{
  ntics = (Durée(date.heure) + Durée::jours(date.jour.nb_jours_debut_ère())).tics;
}

entier DateHeure::microsecondes() const
{
  retourne (entier) (ntics % TICS_SECONDE / TICS_US);
}

bouléen année_mois_valide(entier année, entier mois)
{
  si(!année_est_valide(année))
    retourne non;
  retourne (mois >= 1) && (mois <= 12);
}

entier mois_nb_jours(entier année, entier mois)
{
  si(!année_mois_valide(année, mois))
    échec("mois_nb_jours: année ({}) ou mois ({}) invalide.", année, mois);

  soit ptr = est_bissextile(année) ? nb_jours_par_mois[1] : nb_jours_par_mois[0];
  retourne ptr[mois-1];
}

bouléen Calendrier::est_valide() const
{
  retourne année_mois_valide(année, mois)
      &&  (jour >= 1) && (jour <= mois_nb_jours(année, mois));
}


static vector<string> split(cstring str, cstring delim)
{
  vector<string> tokens;
  size_t prev = 0, pos = 0;
  do
  {
    pos = str.find(delim, prev);
    si(pos == string::npos)
      pos = str.length();
    soit token = str.substr(prev, pos-prev);
    si(!token.empty())
      tokens.push_back(token);
    prev = pos + delim.length();
  }
  tantque (pos < str.length() && prev < str.length());
  retourne tokens;
}

Calendrier::Calendrier(entier année, entier mois, entier jour)
{
  this->année = année;
  this->mois  = mois;
  this->jour  = jour;
}

Calendrier::Calendrier(cstring s)
{
  soit sp = split(s, "/");
  si(sp.size() != 3)
    msg_erreur("Date invalide : {}", s);
  sinon
  {
    jour  = stoi(sp[0]);
    mois  = stoi(sp[1]);
    année = stoi(sp[2]);
  }
  si(!est_valide())
    msg_erreur("Date invalide : {}", s);
}

HeureComposite::HeureComposite(entier heure, entier minutes, entier secondes,
                               entier ms, entier µs)
{
  this->heure     = heure;
  this->minutes   = minutes;
  this->secondes  = secondes;
  this->ms        = ms;
  this->µs        = µs;
}

HeureComposite::HeureComposite(cstring s)
{
  soit sp = split(s, ":");
  si(sp.size() != 3)
    échec("HeureComposite : chaine invalide : {}", s);
  sinon
  {
    heure    = stoi(sp[0]);
    minutes  = stoi(sp[1]);
    secondes = stoi(sp[2]);
  }
  si(!vérifie_validité())
    échec("Heure composite invalide.");
}


bouléen HeureComposite::vérifie_validité() const
{
  retourne (heure >= 0) && (heure < 24)
      && (minutes >= 0) && (minutes < 60)
      && (secondes >= 0) && (secondes < 60)
      && (ms >= 0) && (ms < 1e3)
      && (µs >= 0) && (µs < 1e3);
}


entier Calendrier::nb_jours_debut_année() const
{
  si(!est_valide())
  {
    msg_avert("Date de calendrier invalide ({}/{}/{})", année, mois, jour);
    retourne jour;
  }

  retourne (jour - 1) + cumul_jours_par_mois[est_bissextile(année) ? 1 : 0][mois-1];
}

double grégorien_vers_jours(entier année, double jour_année)
{
  int64_t a1 = année - 1;

  retourne
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

entier Calendrier::nb_jours_debut_ère() const
{
  retourne grégorien_vers_jours(année, 1 + nb_jours_debut_année());
}


// Pré-calcul
static const DateHeure cst_epoque_unix = DateHeure({{1970,1,1}, {0, 0, 0}});

DateHeure DateHeure::epoque_unix()
{
  retourne cst_epoque_unix;
}

DateHeure DateHeure::maintenant()
{
  retourne epoque_unix() + Durée::microsecondes(
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
  retourne Durée(ntics).nb_jours() + 1721425.5;
}

double DateHeure::J2000() const
{
  retourne nb_jours_Julien() - 2415020.0;
}

DateHeure DateHeure::epoque_Matlab()
{
  // 01/01/2000, 00h00
  retourne DateHeure{{{2000, 1, 1}, {0, 0, 0}}};
}

DateHeure DateHeure::epoque_J2000()
{
  // 01/01/2000, 12h00
  retourne DateHeure{{{2000, 1, 1}, {12, 0, 0}}};
}

DateHeure DateHeure::epoque_GPS()
{
  // 06/01/1980, 00h00
  retourne DateHeure({{1980, 1, 6}, {0, 0, 0}});
}

double DateHeure::temps_sidéral_Greenwich() const
{
  // https://lweb.cfa.harvard.edu/~jzhao/times.html

  soit nj = nb_jours_Julien();

  // Nb jours Juliens du minuit précédent
  soit jd0 = floor(nj + 0.5) - 0.5;
  soit t   = (jd0 - 2451545.0) / 36525.0;
  soit jdf = nj - jd0;
  soit gt  = 24110.54841 + t * (8640184.812866 + t * (0.093104 - t * 6.2E-6));
  gt  += jdf * 1.00273790935 * 86400.0;
  retourne modulo_2π(deg2rad(gt * 360.0 / nbsecs_par_jour));
}


bouléen est_bissextile(entier année)
{
  si(!année_est_valide(année))
    retourne non;
  retourne (((année % 4) == 0 && (année % 100) != 0) || (année % 400) == 0);
}

bouléen année_est_valide(entier année)
{
  retourne (année >= 1) && (année <= 9999);
}


HeureComposite DateHeure::decomp_heure() const
{
  retourne
      {(entier) ((ntics % TICS_JOUR)     / TICS_HEURE),
       (entier) ((ntics % TICS_HEURE)    / TICS_MINUTE),
       (entier) ((ntics % TICS_MINUTE)   / TICS_SECONDE),
       (entier) ((ntics % TICS_SECONDE)  / TICS_MS),
       (entier) ((ntics % TICS_MS)       / TICS_US)};
}

DateComposite DateHeure::decomposition() const
{
  retourne {calendrier(), decomp_heure()};
}

DateHeure DateHeure::lis_chaine(cstring s_)
{
  entier y,M,d,h,m,s;
  sscanf(s_.c_str(), "%d-%d-%d %d:%d:%d", &y,&M,&d,&h,&m,&s);
  retourne DateHeure{{{y,M,d},{h,m,s}}};
}

DateHeure DateHeure::de_GPS(entier semaine, entier secs)
{
  // 5 Janvier 1980
  soit r = epoque_GPS();
  r += Durée::secondes(secs);
  r += Durée::jours(7 * semaine);
  r += Durée::secondes(-18);
  retourne r;
}

tuple<entier, entier> DateHeure::vers_GPS() const
{
  // 5 Janvier 1980
  soit e0 = epoque_GPS();

  soit nsecs = (*this - e0).nb_secondes();

  // 06/2017 : 18 leap seconds ahead of UTC
  nsecs += 18;

  entier sps         = nbsecs_par_jour * 7,
         nb_semaines = nsecs / sps,
         nb_secs     = nsecs - nb_semaines * sps;

  retourne {nb_semaines, nb_secs};
}

DateComposite DateHeure::decomposition_locale() const
{
  // timezone = secondes à l'ouest de Greenwich
  // tzset();
  //msg("timezone = {}", timezone);

  /*static entier decalage = -1000;

  si(decalage == -1000)
  {
    std::time_t t = std::time(nullptr);
    struct tm t1, t2;
    t1 = *(std::gmtime(&t)); // UTC
    t2 = *(std::localtime(&t));

    t2.tm_isdst = 1;

    decalage = difftime(mktime(&t1), mktime(&t2));

    msg("\n\n\n*** Décalage horaire locale vs UTC : {}\n\n\n", decalage);
  }*/

  //soit s = std::chrono::sys_info.offset;
  //entier decalage = timezone;
  //decalage = 0;
  //using namespace std::chrono;
  //static entier decalage = std::chrono::get_tzdb().current_zone()->get_info(system_clock::now()).offset;
  //retourne ajoute_secondes(-decalage).decomposition();

  // TODO!!!
  //retourne (*this + Durée::heures(2)).decomposition();
  retourne (*this + Durée::heures(1)).decomposition();
}

Calendrier DateHeure::calendrier() const
{
  entier nb_jours = ntics / TICS_JOUR;

  // Décompte les cycles de 400 ans
  entier ncycles_400_ans = nb_jours / nb_jours_cycle_400_ans;
  nb_jours %= nb_jours_cycle_400_ans;

  // Décompte les cycles de 100 ans
  // 1,2,3,...99,  100,101,...,199, ....
  entier ncycles_100_ans = nb_jours / nb_jours_cycle_100_ans;
  si (ncycles_100_ans == 4)
      // 31 décembre de la dernière année bissextile
      ncycles_100_ans = 3;
  nb_jours -= ncycles_100_ans * nb_jours_cycle_100_ans;

  // Décompte les cycles de 4 ans
  entier ncycles_4_ans = nb_jours / nb_jours_cycle_4_ans;
  nb_jours %= nb_jours_cycle_4_ans;

  // Décompte les années qui restent (dans le cycle de 4 ans en cours)
  // Un cycle de 4 ans =
  // 1, 2, 3, 4
  // ...
  // 2021, 2022, 2023, 2024
  // soit en général :
  // [365, 365, 365, 366]
  // 0 <= nb jours < 3 * 365 + 366
  entier nb_ans = 0;
  si(nb_jours < 365)
    ;
  sinon si(nb_jours < 2 * 365)
  {
    nb_ans = 1;
    nb_jours -= 365;
  }
  sinon si(nb_jours < 3 * 365)
  {
    nb_ans = 2;
    nb_jours -= 2 * 365;
  }
  sinon
  {
    nb_ans = 3;
    nb_jours -= 3 * 365;
  }

  entier année = (ncycles_400_ans * 400) + (ncycles_100_ans * 100)
            + (ncycles_4_ans * 4) + nb_ans + 1;

  // Nombre de jours pour chaque mois
  soit ptr = est_bissextile(année) ? nb_jours_par_mois[1] : nb_jours_par_mois[0];

  entier mois = 0;
  tantque((nb_jours >= ptr[mois]) && (mois <= 11))
    nb_jours -= ptr[mois++];

  retourne {année, mois + 1, nb_jours + 1};
}

bouléen mode_utc = non;

std::ostream& operator<<(std::ostream& ss, const Calendrier &date)
{
  ss << sformat("{:0>4d}-{:0>2d}-{:0>2d}", date.année, date.mois, date.jour);
  retourne ss;
}


std::ostream& operator<<(std::ostream& ss, const HeureComposite &hc)
{
  ss << sformat("{:0>2d}:{:0>2d}:{:0>2d},{:0>3d}:{:0>3d}",
      hc.heure, hc.minutes, hc.secondes, hc.ms, hc.µs);
  retourne ss;
}

std::ostream& operator<<(std::ostream& ss, const DateHeure &t)
{
  DateComposite dc;
  si(mode_utc)
    dc = t.decomposition();
  sinon
    dc = t.decomposition_locale();

  ss << sformat("{:0>4d}-{:0>2d}-{:0>2d} {:0>2d}:{:0>2d}:{:0>2d}",
      dc.jour.année, dc.jour.mois, dc.jour.jour,
      dc.heure.heure, dc.heure.minutes, dc.heure.secondes);

  ss << sformat(":{:0>3d}", dc.heure.ms);

  si(mode_utc)
    ss << " (UTC)";
  sinon
    ss << " (local hour)";
  retourne ss;
}
}
