#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"

namespace tsd::vue::unites {



std::string valeur_vers_chaine(double t, const std::string &unite, int nb_chiffres)
{
  auto [expo, nc] = calc_expo_nb_chiffres(t, unite);
    return valeur_vers_chaine(t, unite, expo, nb_chiffres);
}

std::string valeur_vers_chaine(double t, const std::string &unite)
{
  auto [expo, nc] = calc_expo_nb_chiffres(t, unite);
  return valeur_vers_chaine(t, unite, expo, nc);
}

std::string valeur_vers_chaine(double t, const std::string &unite, int expo, int nb_chiffres)
{
  std::string un = unite;

  if(!unite.empty())
  {
    if(expo == 3)
      un = "K" + un;
    else if(expo == 6)
      un = "M" + un;
    else if(expo == 9)
      un = "G" + un;
    else if(expo == -3)
      un = "m" + un;
    else if(expo == -6)
      un = "u" + un;
    else if(expo == 0)
      ;
    else
      un = "?" + un;
  }
  else
  {
    if(expo != 0)
      un = fmt::format("e{}", expo);
  }

  t *= std::pow(10, -expo);

  if(!un.empty())
    un = " " + un;

  if(nb_chiffres == 0)
    return fmt::format("{}{}", (int) t, un);
  else
  {
    char buf[50];
    sprintf(buf, "{:.%df}{}", nb_chiffres);
    return fmt::format(FMT_RUNTIME(std::string(buf)), t, un);
  }
}


static int ndigits(double a)
{
  for(auto i = 0; i < 8; i++)
  {
    float at = a * std::pow(10, i);
    //if(at == std::floor(at))
    if(std::abs(at - std::round(at)) < 2 * std::pow(10, i) * std::numeric_limits<float>::epsilon())
    {
      //msg("ndigits({}) = {} (at = {}, at-round(at) = {})", a, i, at, std::abs(at - std::round(at)));
      return i;
    }
  }
  //msg("ndigits({}) = 8 !", a);
  return 8;
}


int ndigits(double t, int expo, const std::string &unite)
{
  float at = std::abs(t);
  return ndigits(at * std::pow((double) 10, (double) -expo));
}

// unité : unité de base,
// t : valeur
// retourne : {exposant, nb chiffres significatifs}
std::tuple<int,int> calc_expo_nb_chiffres(double t, const std::string &unite)
{
  // TODO : ici ndigits calculé en doublon


  float at = std::abs(t);

  if(unite.empty())
  {
    if((std::abs(t) < 70000) && (std::abs(t) >= 1))
      return {0, ndigits(at)};
    else if(t == 0)
      return {0, 0};
    else if(std::abs(t) >= 0.1)
      return {0, ndigits(at)};
    else
    {
      int pow = floor(std::log10(at));
      return {pow, ndigits(at * std::pow((double) 10, (double) -pow))}; // A VERIFIER
    }
  }

  if((at < 1000) && (at >= 1))
    return {0, ndigits(at)};
  else if((at >= 1000) && (at < 1e6))
    return {3, ndigits(at * 1e-3)};
  else if((at >= 1e6) && (at < 1e9))
    return {6, ndigits(at * 1e-6)};
  else if(at >= 1e9)
    return {9, ndigits(at * 1e-9)};
  else if(t == 0)
    return {0, 0};
  else if(at < 1e-3)
    return {-6, ndigits(at * 1e6)};
  else if(at < 1)
    return {-3, ndigits(at * 1e3)};
  else
    return {0, ndigits(at)};
}


std::tuple<int, int> calc_expo_nb_chiffres_commun(const std::vector<double> &tics, const std::string &unite)
{
  if(tics.empty())
    return {0, 0};

  bool has_max = false;

  int expo_max = 0;
  //auto [expo_max, nbchiffres_max] = calc_expo_nb_chiffres(tics[0], unite);
  for(auto i = 0u; i < tics.size(); i++)
  {
    // TODO : ici nbchiffres inutiles
    auto [expo, nbchiffres] = calc_expo_nb_chiffres(tics[i], unite);
    if(tics[i] != 0)
    {
      if(has_max)
        expo_max        = std::min(expo_max, expo);
      else
        expo_max = expo;
    }
  }

  int nbchiffres_max = ndigits(tics[0], expo_max, unite);
  //msg("ndigits({}, {}) = {}", tics[0], expo_max, nbchiffres_max);
  for(auto i = 1u; i < tics.size(); i++)
  {
    auto nbchiffres = ndigits(tics[i], expo_max, unite);
    //msg("ndigits({}, {}) = {}", tics[i], expo_max, nbchiffres);
    nbchiffres_max  = std::max(nbchiffres_max, nbchiffres);
  }

  return {expo_max, nbchiffres_max};
}



}
