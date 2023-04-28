#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"

namespace tsd::vue::unites {



string valeur_vers_chaine(double t, cstring unite, entier nb_chiffres)
{
  soit [expo, nc] = calc_expo_nb_chiffres(t, unite);
  retourne valeur_vers_chaine(t, unite, expo, nb_chiffres);
}

string valeur_vers_chaine(double t, cstring unite)
{
  soit [expo, nc] = calc_expo_nb_chiffres(t, unite);
  retourne valeur_vers_chaine(t, unite, expo, nc);
}

string valeur_vers_chaine(double t, cstring unite, entier expo, entier nb_chiffres)
{
  string un = unite;

  si(!unite.empty())
  {
    si(expo == 3)
      un = "K" + un;
    sinon si(expo == 6)
      un = "M" + un;
    sinon si(expo == 9)
      un = "G" + un;
    sinon si(expo == -3)
      un = "m" + un;
    sinon si(expo == -6)
      un = "u" + un;
    sinon si(expo == 0)
      ;
    sinon
      un = "?" + un;
  }
  sinon
  {
    si(expo != 0)
      un = sformat("e{}", expo);
  }

  t *= pow(10, -expo);

  si(!un.empty())
    un = " " + un;

  si(nb_chiffres == 0)
    retourne sformat("{}{}", (entier) round(t), un);
  sinon
  {
    char buf[50];
    sprintf(buf, "{:.%df}{}", nb_chiffres);
    retourne sformat(FMT_RUNTIME(string(buf)), t, un);
  }
}


static entier ndigits(double a)
{
  pour(auto i = 0; i < 8; i++)
  {
    float at = a * pow(10, i);
    si(abs(at - round(at)) < 2 * pow(10, i) * std::numeric_limits<float>::epsilon())
    {
      retourne i;
    }
  }
  retourne 8;
}


entier ndigits(double t, entier expo, cstring unite)
{
  soit at = abs(t);
  retourne ndigits(at * pow(10.0, (double) -expo));
}

// unité : unité de base,
// t : valeur
// retourne : {exposant, nb chiffres significatifs}
tuple<entier,entier> calc_expo_nb_chiffres(double t, cstring unite)
{
  // TODO : ici ndigits calculé en doublon
  float at = abs(t);

  si(unite.empty())
  {
    si((abs(t) < 70000) && (abs(t) >= 1))
      retourne {0, ndigits(at)};
    sinon si(t == 0)
      retourne {0, 0};
    sinon si(abs(t) >= 0.1)
      retourne {0, ndigits(at)};
    sinon
    {
      entier p = floor(log10(at));
      retourne {p, ndigits(at * pow(10.0, (double) -p))}; // A VERIFIER
    }
  }

  si((at < 1000) && (at >= 1))
    retourne {0, ndigits(at)};
  sinon si((at >= 1000) && (at < 1e6))
    retourne {3, ndigits(at * 1e-3)};
  sinon si((at >= 1e6) && (at < 1e9))
    retourne {6, ndigits(at * 1e-6)};
  sinon si(at >= 1e9)
    retourne {9, ndigits(at * 1e-9)};
  sinon si(t == 0)
    retourne {0, 0};
  sinon si(at < 1e-3)
    retourne {-6, ndigits(at * 1e6)};
  sinon si(at < 1)
    retourne {-3, ndigits(at * 1e3)};
  sinon
    retourne {0, ndigits(at)};
}


tuple<entier, entier>
  calc_expo_nb_chiffres_commun(const vector<double> &tics, cstring unite)
{
  si(tics.empty())
    retourne {0, 0};

  bouléen has_max = non;

  entier expo_max = 0;
  //soit [expo_max, nbchiffres_max] = calc_expo_nb_chiffres(tics[0], unite);
  pour(auto i = 0u; i < tics.size(); i++)
  {
    // TODO : ici nbchiffres inutiles
    soit [expo, nbchiffres] = calc_expo_nb_chiffres(tics[i], unite);
    si(tics[i] != 0)
    {
      si(has_max)
        expo_max  = min(expo_max, expo);
      sinon
        expo_max = expo;
    }
  }

  soit nbchiffres_max = ndigits(tics[0], expo_max, unite);
  pour(auto i = 1u; i < tics.size(); i++)
  {
    soit nbchiffres = ndigits(tics[i], expo_max, unite);
    nbchiffres_max  = max(nbchiffres_max, nbchiffres);
  }

  retourne {expo_max, nbchiffres_max};
}



}
