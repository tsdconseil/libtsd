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


extern void tsd_log_msg(const char *fn, entier ligne, entier niveau, cstring str);

namespace tsd {



bouléen tests_debug_actif = non;

extern bouléen erreur_attendue;

void vérifie_exception(fonction<void()> func)
{
  erreur_attendue = oui;
  bouléen erreur_détectée = non;
  try
  {
    func();
  }
  catch(...)
  {
    erreur_détectée = oui;
    msg("Erreur bien détectée.");
  }
  erreur_attendue = non;
  assertion_msg(erreur_détectée, "Une exception était attendue, elle n'a pas été détectée.");
}


void verifie_erreur_relative(float v, float ref, float precision, cstring refname)
{
  // TODO: not oui
  si((ref < 0.000001) && (v < 0.000001))
    retourne;

  si(((ref == 0.0) || (ref == -0.0)) && (v == 0.0))
    retourne;

  float err = 100.0 * abs((v - ref) / ref);

  si(err > precision)
    échec("{}: erreur trop grande. Valeur = {}, référence = {}, erreur relative = {} %, erreur relative max = {} %.",
        refname.c_str(), v, ref, err, precision);
}


static uint64_t get_tick_count_us()
{
  struct timespec ts;

  si(clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
  {
    perror("clock_gettime().");
    retourne 0;
  }
  retourne (uint64_t) (ts.tv_nsec / 1000) + (((uint64_t) ts.tv_sec) * 1000 * 1000);
}


static entier teste(entier i, vector<Test> &tests)
{
  soit &t = tests[i];
  fmt::print("\033[34;1mTest");
  fmt::print(" [{}/{}]", (i+1), tests.size());


  tsd::vue::stdo.def_dossier_sortie("./build/test-log/" + t.nom);

  fmt::print(" \033[34;1m{}\033[0m...", t.nom);
  std::cout << std::endl;

  float t0 = get_tick_count_us();

  entier res = -1;
  try
  {
    t.fonction();
    res = 0;
  }
  catch(std::runtime_error &e)
  {
    msg("Exception : {}", e.what());
    res = -1;
  }

  float t1 = get_tick_count_us();

  float dms = (t1 - t0) / 1000.0;

  si(res == 0)
  {
    fmt::print("  \033[32msuccès,\033[0m durée = {:.2f} ms.\n", dms);
  }
  sinon
  {
    fmt::print("  \033[31mEchec test {},\033[0m durée = {:.2f} ms.\n", t.nom, dms);
    fmt::print("Abandon des tests.\n");
  }

  tsd::vue::stdo.flush();

  retourne res;
}




entier fait_tests(entier argc, const char *argv[], vector<Test> &tests)
{
# ifdef WIN
  // pour éviter la fenêtre Windows en cas d'échec d'assertion
  // (celle-ci empêche le point d'arrêt gdb)
  _set_error_mode(_OUT_TO_STDERR);
# endif

  entier opt;
  msg("Tests automatiques libtsd...");

  get_logger() = tsd_log_msg;

  tsd::vue::stdo.def_dossier_sortie("./build/test-log");

  tantque ((opt = getopt(argc, (char * const *) argv, "t:hald")) != -1)
  {
    string nom;

    si(optarg)
      nom = optarg;

    switch (opt)
    {
      case 'd':
        tests_debug_actif = oui;
        break;
      case 't':
      {
        pour(auto i = 0u; i < tests.size(); i++)
        {
          si(tests[i].nom == nom)
          {
            entier res = teste(i, tests);
            fmt::print("Fin des tests.\n");
            retourne res;
          }
        }
        msg_erreur("Test non trouvé : [{}]", nom);
        retourne -1;
      }
      case 'h':
      case 'a':
      case 'l':
      {
        fmt::print("Liste des tests:\n");
        pour(auto &t: tests)
          fmt::print(" - {}\n", t.nom);
        retourne 0;
      }
    }
  }

  pour(auto i = 0u; i < tests.size(); i++)
  {
    entier res = teste(i, tests);
    si(res)
      retourne res;
  }

  fmt::print("Fin des tests.\n");

  tsd::vue::stdo.fin();

  retourne 0;
}
}


