#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"


namespace tsd {

class MoniteurCpu
{
public:
  struct Stats
  {
    string nom;
    float conso_cpu_pourcents = 0;
    entier   nb_appels = 0;
  };

  MoniteurCpu(cstring nom = "");
  string &nom();
  void commence_op();
  void fin_op();
  void reset();
  Stats stats() const;

  _PIMPL_
};

struct MoniteursStats
{
  vector<MoniteurCpu::Stats> lst;
  void ajoute(MoniteurCpu &m);
  MoniteurCpu::Stats get(cstring nom) const;
};


extern MoniteurCpu moniteur_spectrum;

}
