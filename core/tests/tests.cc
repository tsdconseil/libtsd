#include "tsd/tests.hpp"

using namespace tsd;
using namespace std;


extern test_routine_t test_detecteur, test_unites, test_filtre_fft, test_crec, test_clkrec,
  test_tsd, test_wav, test_demod, test_filtres, test_figure, test_fourier, test_ra,
  test_fenetres, test_tod, test_cqt, test_image, test_poly, test_stats, test_telecom, test_kalman,
  test_geometrie, test_date_heure, test_itrp, test_delais_filtres, test_filtre_adapte, test_formes_ondes,
  test_itrp_irreg, test_recepteur, bench_recepteur, test_motifs, test_bitstream,
  test_filtrage_analyse, test_dsp;


static std::vector<Test> tests =
{
  {"dsp",           &test_dsp},
  {"tsd",           &test_tsd},
  {"bitstream",     &test_bitstream},
  {"motifs",        &test_motifs},
  {"date-heure",    &test_date_heure},
  {"geo",           &test_geometrie},
  {"kalman",        &test_kalman},
  {"stats",         &test_stats},
  {"poly",          &test_poly},
  {"cqt",           &test_cqt},
  {"tod",           &test_tod},
  {"fenetres",      &test_fenetres},
  {"ra",            &test_ra},
  {"fourier",       &test_fourier},
  {"figure",        &test_figure},
  {"unites",        &test_unites},
  {"image",         &test_image},
  {"demod",         &test_demod},
  {"wav",           &test_wav},
  {"fft-filtre",    &test_filtre_fft},
  {"detecteur",     &test_detecteur},
  {"filtres",          &test_filtres},
  {"filtrage-analyse", &test_filtrage_analyse},
  {"crec",          &test_crec},
  {"telecom-fo",    &test_formes_ondes},
  {"telecom-recepteur", &test_recepteur},
  {"telecom-recepteur-bench", &bench_recepteur},
  {"clkrec",        &test_clkrec},
  {"itrp",          &test_itrp},
  {"itrpi",         &test_itrp_irreg},
  {"dfiltre",       &test_delais_filtres},
  {"telecom",       &test_telecom},
  {"fa",            &test_filtre_adapte},
};


int main(int argc, const char **argv)
{
  return fait_tests(argc, argv, tests);
}

