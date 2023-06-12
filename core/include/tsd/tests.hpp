#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include <functional>

#include "tsd/tsd.hpp"

namespace tsd {

typedef void (test_routine_t)();

struct Test
{
  string nom;
  test_routine_t *fonction;
};


extern bouléen tests_debug_actif;

extern entier fait_tests(entier argc, const char *argv[], vector<Test> &tests);


extern void verifie_erreur_relative(float v, float ref, float precision, cstring refname);

extern void vérifie_exception(fonction<void()> func);

}

