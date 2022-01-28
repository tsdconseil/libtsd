#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"


namespace tsd {

class MoniteurCpu
{
public:
  struct Stats
  {
    std::string nom;
    float conso_cpu_pourcents = 0;
    int   nb_appels = 0;
  };

  MoniteurCpu(const std::string &nom = "");
  std::string &nom();
  void commence_op();
  void fin_op();
  void reset();
  Stats stats() const;
private:
  struct Impl;
  sptr<Impl> impl;
};

struct MoniteursStats
{
  std::vector<MoniteurCpu::Stats> lst;
  void ajoute(MoniteurCpu &m);
  MoniteurCpu::Stats get(const std::string &nom) const;
};


extern MoniteurCpu moniteur_spectrum;

}
