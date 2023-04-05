#include "tsd/moniteur-cpu.hpp"

#include <ctime>
#include <deque>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <mutex>
#include <map>
#include <thread>
#ifdef WIN
# include <windows.h>
#endif

using namespace std;

namespace tsd {


MoniteurCpu moniteur_spectrum{"spectrum"}, moniteur_fft{"fft"};

int64_t tic_μs(bouléen monotonique = oui)
{
  // Avec MSYS2, sous Windows, les compteurs Posix ne sont pas suffisament prècis
# ifdef WIN
  static LARGE_INTEGER base_tick = {0}, frequency = {0};
  static bouléen tick_init_done = non;
  si(!tick_init_done)
  {
    si(!QueryPerformanceFrequency(&frequency))
    {
      echec("Failed to initialize 64 bits counter.");
      frequency.QuadPart = 1000 * 1000;
    }
    QueryPerformanceCounter(&base_tick);
    tick_init_done = oui;
  }
  LARGE_INTEGER tick;
  QueryPerformanceCounter(&tick);
  retourne (int64_t) ((tick.QuadPart-base_tick.QuadPart)*1000.0*1000.0 / frequency.QuadPart);
# elif USE_STD_CLOCK
  retourne clock() * (1e6 / CLOCKS_PER_SEC);
# else
  struct timespec ts;
  clock_gettime(monotonique ? CLOCK_MONOTONIC : CLOCK_THREAD_CPUTIME_ID, &ts);
  retourne ts.tv_sec * 1e6 + ts.tv_nsec * 1e-3;
# endif
}




struct MoniteurCpu::Impl
{
  string nom;
  mutable mutex mut;

  struct PerThread
  {
    bouléen en_cours     = non;

    // cl0  : par thread
    // cl00 : monotonique
    int64_t cl0, cl00;

    bouléen cl00_init = non;

    int64_t total_μs  = 0;
    int64_t μs_en_cours = 0;
    entier cnt     = 0, total_cnt = 0;
    float pourcent_cpu   = 0;

    float last_add = 0;
  };


  map<thread::id, PerThread> pt;

  void commence_op()
  {
    soit &p = pt[this_thread::get_id()];
    si(p.en_cours)
      msg_avert("Moniteur [{}] : deja en cours.", nom);
    p.en_cours  = oui;

    p.cl0 = tic_μs(non);

    si(!p.cl00_init)
    {
      p.cl00 = tic_μs(oui);
      p.cl00_init = oui;
    }
  }
  void fin_op()
  {
    soit &p = pt[this_thread::get_id()];
    si(!p.en_cours)
      msg_avert("Moniteur [{}] : pas en cours.", nom);
    p.en_cours = non;

    // Attention, comptage à la ms près seulement !!!

    soit now_thread = tic_μs(non);
    soit now_mono   = tic_μs(oui);

    soit df = now_thread - p.cl0;


    //msg("df = {}", df);

    p.total_cnt++;
    p.total_μs    += df;
    p.μs_en_cours += df;
    p.cnt++;

    /*si(df != 0)
    {
      msg("DF = {}, us en cours = {}", df, p.us_en_cours);
    }*/

    soit df0 = (now_mono - p.cl00) * 1e-6f;

    // Mise à jour toute les 500 ms
    si(df0 > 0.5)
    {
      p.pourcent_cpu = (1e2f * p.μs_en_cours * 1e-6) / df0;
      //msg("MAJ {}: PCPU = {} (EN COURS = {}, TOT = {}, cl00={}, df0={})", nom, p.pourcent_cpu, p.us_en_cours, p.total_us, (float) now_mono, (float) df0);
      p.cnt          = 0;
      p.μs_en_cours  = 0;
      p.cl00         = now_mono;
    }
  }

  // PB : si stats pas mis à jour depuis longtemps par un thread...
  Stats lis_stats() const
  {
    Stats res;
    res.nom = nom;

    soit now = tic_μs(oui);
    entier nta = 0;

    mut.lock();
    pour(auto &p: pt)
    {
      soit df = (now - p.second.cl00) * 1e-6f;
      si(df > 1)
      {
        // Thread plus actif.
        //msg("Thread plus actif : {}, df = {}", nom, df);
      }
      sinon
      {
        nta++;
        res.conso_cpu_pourcents += p.second.pourcent_cpu;
        res.nb_appels           += p.second.total_cnt;
      }
    }
    mut.unlock();

    //msg("Stats[{}]: cpu = {:e} %, nthreads actifs : {}, nb appels : {}", nom, res.conso_cpu_pourcents, nta, res.nb_appels);

    retourne res;
  }

  /*void affiche_stats()
  {
    printf("%s: %d appels, temps cpu total = %.3f ms\n",
        nom.c_str(), total_cnt,
        ((float) total_us) * 1e-3f);
  }*/

  void reset()
  {
    //msg("Reset mon cpu.");
    pour(auto &p: pt)
      p.second = PerThread();
  }
};



string &MoniteurCpu::nom()
{
  retourne impl->nom;
}

MoniteurCpu::MoniteurCpu(const string &nom)
{
  impl = make_shared<Impl>();
  impl->nom = nom;
}

void MoniteurCpu::commence_op()
{
  impl->commence_op();
}

void MoniteurCpu::fin_op()
{
  impl->fin_op();
}

/*void MoniteurCpu::affiche_stats()
{
  impl->affiche_stats();
}*/

MoniteurCpu::Stats MoniteurCpu::stats() const
{
  retourne impl->lis_stats();
}

void MoniteurCpu::reset()
{
  impl->reset();
}

void MoniteursStats::ajoute(MoniteurCpu &m)
{
  lst.push_back(m.stats());
}

MoniteurCpu::Stats MoniteursStats::get(const string &nom) const
{
  pour(auto &s: lst)
  {
    si(s.nom == nom)
      retourne s;
  }
  msg_avert("Moniteur : stat non trouvée [{}]", nom);
  retourne MoniteurCpu::Stats();
}

}

