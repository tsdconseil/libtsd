#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"
#include "tsd/tests.hpp"

extern "C"
{
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
}

#include <string>
#include <vector>
#include <chrono>
#include <iostream>

namespace tsd {



bool tests_debug_actif = false;

extern bool erreur_attendue;

void vérifie_exception(std::function<void()> func)
{
  erreur_attendue = true;
  bool erreur_détectée = false;
  try
  {
    func();
  }
  catch(...)
  {
    erreur_détectée = true;
    msg("Erreur bien détectée.");
  }
  erreur_attendue = false;
  tsd_assert_msg(erreur_détectée, "Une exception était attendue, elle n'a pas été détectée.");
}


int verifie_erreur_relative(float v, float ref, float precision, const std::string &refname)
{
  // TODO: not true
  if((ref < 0.000001) && (v < 0.000001))
    return 0;

  if(((ref == 0.0) || (ref == -0.0)) && (v == 0.0))
    return 0;

  float err = 100.0 * std::abs((v - ref) / ref);

  if(err > precision)
  {
    msg_erreur("{}: erreur trop grande. Valeur = {}, référence = {}, erreur relative = {} %, erreur relative max = {} %.",
        refname.c_str(), v, ref, err, precision);
    return -1;
  }

  return 0;
}


static uint64_t get_tick_count_us()
{
  struct timespec ts;

  if(clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
  {
    perror("clock_gettime().");
    return 0;
  }
  return (uint64_t) (ts.tv_nsec / 1000) + (((uint64_t) ts.tv_sec) * 1000 * 1000);
}


static int teste(int i, std::vector<Test> &tests)
{
  auto &t = tests[i];
  fmt::print("\033[34;1mTest");
  fmt::print(" [{}/{}]", (i+1), tests.size());

  tsd::vue::stdo.def_dossier_sortie("./build/test-log/" + t.nom);

  fmt::print(" \033[34;1m{}\033[0m...", t.nom);
  std::cout << std::endl;

  float t0 = get_tick_count_us();

  int res = -1;
  try
  {
    res = t.fonction();
  }
  catch(std::runtime_error &e)
  {
    msg("Exception : {}", e.what());
  }

  float t1 = get_tick_count_us();

  float dms = (t1 - t0) / 1000.0;

  if(res == 0)
  {
    fmt::print("  \033[32msuccès,\033[0m durée = {:.2f} ms.\n", dms);
  }
  else
  {
    fmt::print("  \033[31mEchec test {},\033[0m durée = {:.2f} ms.\n", t.nom, dms);
    fmt::print("Abandon des tests.\n");
  }

  tsd::vue::stdo.flush();

  return res;
}




int fait_tests(int argc, const char *argv[], std::vector<Test> &tests)
{
# ifdef WIN
  // Pour éviter la fenêtre Windows en cas d'échec d'assertion
  // (celle-ci empêche le point d'arrêt gdb)
  _set_error_mode(_OUT_TO_STDERR);
# endif

  int opt;
  msg("Tests automatiques libtsd...");

  tsd::vue::stdo.def_dossier_sortie("./build/test-log");

  while ((opt = getopt(argc, (char * const *) argv, "t:hald")) != -1)
  {
    std::string nom;

    if(optarg != nullptr)
      nom = optarg;

    switch (opt)
    {
      case 'd':
        tests_debug_actif = true;
        break;
      case 't':
      {
        for(auto i = 0u; i < tests.size(); i++)
        {
          if(tests[i].nom == nom)
          {
            int res = teste(i, tests);
            fmt::print("Fin des tests.\n");
            return res;
          }
        }
        msg_erreur("Test non trouvé : [{}]", nom);
        return -1;
      }
      case 'h':
      case 'a':
      case 'l':
      {
        fmt::print("Liste des tests:\n");
        for(auto &t: tests)
          fmt::print(" - {}\n", t.nom);
        return 0;
      }
    }
  }

  for(auto i = 0u; i < tests.size(); i++)
  {
    int res = teste(i, tests);
    if(res)
      return res;
  }


  fmt::print("Fin des tests.\n");

  tsd::vue::stdo.fin();

  return 0;
}
}


