#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include <string>
#include <vector>
#include <functional>

#include "tsd/tsd.hpp"

namespace tsd {

typedef int (test_routine_t)();

struct Test
{
  std::string nom;
  test_routine_t *fonction;
};


extern bouléen tests_debug_actif;

extern entier fait_tests(entier argc, const char *argv[], std::vector<Test> &tests);


extern entier verifie_erreur_relative(float v, float ref, float precision, const std::string &refname);

extern void vérifie_exception(std::function<void()> func);

}

