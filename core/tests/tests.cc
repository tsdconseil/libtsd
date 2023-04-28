#include "tsd/tests.hpp"
#include "tsd/tsd.hpp"

using namespace tsd;
using namespace std;


extern test_routine_t test_detecteur, test_unites, test_filtre_fft, test_crec, test_clkrec,
  test_tsd, test_wav, test_demod, test_filtres, test_figure, test_fourier, test_ra,
  test_fenetres, test_tod, test_cqt, test_image, test_poly, test_stats, test_telecom, test_kalman,
  test_geometrie, test_date_heure, test_itrp, test_delais_filtres, test_filtre_adapte, test_formes_ondes,
  test_itrp_irreg, test_recepteur, bench_recepteur, test_motifs, test_bitstream,
  test_filtrage_analyse, test_dsp, test_tab;


static vector<Test> tests =
{
  {"tab",               &test_tab},
  {"tsd",               &test_tsd},
  {"dsp",               &test_dsp},
  {"bitstream",         &test_bitstream},
  {"fourier",           &test_fourier},
  {"date-heure",        &test_date_heure},
  {"geo",               &test_geometrie},
  {"fenetres",          &test_fenetres},
  {"ra",                &test_ra},
  {"figure",            &test_figure},
  {"unites",            &test_unites},
  {"image",             &test_image},
  {"filtres",           &test_filtres},
  {"filtrage-analyse",  &test_filtrage_analyse},
  {"fft-filtre",        &test_filtre_fft},
  {"itrp",              &test_itrp},
  {"itrpi",             &test_itrp_irreg},
  {"motifs",            &test_motifs},
  {"stats",             &test_stats},
  {"poly",              &test_poly},
  {"cqt",               &test_cqt},
  {"tod",               &test_tod},
  {"wav",               &test_wav},
  {"crec",              &test_crec},
  {"telecom-fo",        &test_formes_ondes},
  {"clkrec",            &test_clkrec},
  {"dfiltre",           &test_delais_filtres},
  {"fa",                &test_filtre_adapte},
  {"detecteur",         &test_detecteur},
  {"demod",             &test_demod},
  {"telecom",           &test_telecom},
  {"telecom-recepteur", &test_recepteur},
  {"telecom-recepteur-bench", &bench_recepteur},
  //  {"kalman",        &test_kalman},
};


entier main(entier argc, const char **argv)
{
  retourne fait_tests(argc, argv, tests);
}

