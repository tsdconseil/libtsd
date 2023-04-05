#include "tsd/tsd.hpp"
#include "tsd/filtrage/frat.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include <fmt/color.h>
#include <random>
#include <cstdarg>
#include <iostream>
#include <fmt/chrono.h>
#include <chrono>
#include <stdexcept>

using namespace std;




namespace tsd {


default_random_engine generateur_aleatoire;





Vecf sigscie(entier p, entier n)
{
  retourne Vecf::int_expr(n, IMAP(
      ((i % p) - ((p-1) * 0.5)) / (0.5 * (p-1))));
}

Vecf sigtri(entier p, entier n)
{
  Vecf x(n);
  pour(auto i = 0; i < n; i++)
  {
    entier j = i % p;
    si(j < p/2)
      x(i) = j;
    sinon
      x(i) = p - j;
    x(i) -= 0.5 * (p/2);
    x(i) /= p;
  }
  retourne 4*x; // entre -1 et 1
}

Vecf sigcar(entier p, entier n)
{
  retourne Vecf::int_expr(n, IMAP(
      2 * (((i / (p/2)) % 2) - 0.5)));
}

Vecf sigimp(entier n, entier p)
{
  Vecf x = Vecf::zeros(n);
  tsd_assert_msg((p >= 0) && (p < n), "sigimp(n={},p={}) : p devrait être compris entre 0 et n-1.", n, p);
  x(p) = 1;
  retourne x;
}



Veccf sigexp(float f, entier n)
{
  Veccf x(n);
  cdouble r = 1,
          w0 = std::polar<double>(1.0, f*2*π);

  entier i = 0;
  tantque(i < n)
  {
    pour(auto j = 0; (j < 1000) && (i < n); j++)
    {
      x(i++) = r;
      r *= w0;
    }
    // Tout les 1000 échantillons, corrige la phase
    r /= abs(r);
  }

  retourne x;
}

Vecf sigsin(float f, entier n)
{
  retourne imag(sigexp(f, n));
}

Vecf sigcos(float f, entier n)
{
  retourne real(sigexp(f, n));
}

Vecf signyquist(entier n)
{
  entier i;
  Vecf x(n);
  pour(i = 0; i + 1 < n; i += 2)
  {
    x(i) = -1;
    x(i+1) = 1;
  }
  si(i == n - 1)
    x(i) = -1;
  retourne x;
}

Vecf sigchirp(float f0, float f1, entier n, char mode)
{
  Vecf freq;
  si(mode == 'l')
   freq = linspace(f0, f1, n);
  sinon si(mode == 'q')
    freq  = f0 + (f1 - f0) * square(linspace(0, 1, n));
  sinon
  {
    msg_erreur("sigchirp: mode invalide '{}' (devrait être 'l' ou 'q').", mode);
    retourne {};
  }
  soit phase = (2 * π) * cumsum(freq);
  retourne cos(phase);
}

Vecf siggauss(entier n, float a)
{
  soit t = (linspace(0, n-1, n) - n/2.0f) / (n/2.0f);
  retourne exp(- a * t * t);
}

Vecf siggsin(float f, entier n, float a)
{
  soit x = sigsin(f, n);
  soit t = (linspace(0, n-1, n) - n/2.0f) / (n/2.0f);
  retourne x * exp(- a * t * t);
}

entier prochaine_puissance_de_2(unsigned int i)
{
  entier lg2 = (entier) ceil(log((float) i) / log(2.0f));
  retourne 1l << lg2;
}

uint64_t get_tick_count_µs()
{

# if WIN32
  echec("TODO: VSTUDIO / get_tick_count_µs");
  retourne 0;
# else
  struct timespec ts;
  si(clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
  {
    perror("clock_gettime().");
    retourne 0;
  }
  retourne (ts.tv_nsec / 1000) + (((uint64_t)ts.tv_sec) * 1000 * 1000);
# endif
}




// En mode test
bouléen erreur_attendue = non;


void log_default(const char *fn, const entier ligne, entier niveau, const char *fonction, const string &str)
{
  fmt::text_style drapeaux;

  si(!erreur_attendue)
  {
    si(niveau == 2)
      drapeaux = fmt::emphasis::bold;
    sinon si(niveau == 3)
      drapeaux = fg(fmt::color::crimson);
    sinon si(niveau >= 4)
      drapeaux = fg(fmt::color::crimson) | fmt::emphasis::bold;
  }

  soit tot = get_tick_count_µs();

  static int64_t first = -1, last = -1;
  si(first == -1)
    first = last = tot;

  soit aff = tot - first;
  float diff = tot - last;

  //soit aff = tot - last;


  soit µs  = aff % 1000;
  soit ms  = (aff / 1000) % 1000;
  soit sec = (aff / 1000000);

  string unit = "µs";
  si(diff > 1000)
  {
    diff /= 1000;
    unit = "ms";
  }
  si(diff > 1000)
  {
    diff /= 1000;
    unit = " s";
  }

  soit s1 = fmt::format("[+{: >5.1f} {}]", diff, unit);

  si(unit == " s")
    s1 = "\033[1;32m" + s1 + "\033[0m";
  sinon si(unit == "ms")
    s1 = "\033[1m" + s1 + "\033[0m";



  soit s_dt  = fmt::format("[{:03d},{:03d},{:03d}] {}", sec, ms, µs, s1);
  soit s_msg = fmt::format(" {}\n", str);

  last = tot;

  si(niveau >= 4)
  {
    s_msg += fmt::format("(fichier : {}, ligne : {})\n", fn, ligne);
  }

  tsd::vue::stdo.printf(s_msg.c_str());
  fmt::print(FMT_RUNTIME(s_dt));
  fmt::print(drapeaux, "{}", s_msg);

  fflush(0);
  flush(cout);

  si(niveau >= 4)
  {
    printf("Throw...\n"); fflush(0);
    throw runtime_error(str);
  }
}





static logger_t the_log = log_default;


void reset_logger()
{
  the_log = log_default;
}

void set_logger(logger_t log_)
{
  the_log = log_;
}

void msg_impl2(const char *fn, const entier ligne, entier niveau, const char *fonction, const string &str)
{
  si(the_log)
    the_log(fn, ligne, niveau, fonction, str);
  sinon
    // Car l'affectation log=log_default n'est pas forcément déjà faite ici.
    // (par exemple, si appel depuis l'initialisation d'une variable statique)
    log_default(fn, ligne, niveau, fonction, str);
}

/** @brief Cette classe tamponne les données d'entrée
 *  pour les sortir sous la forme de paquets de taille constante
 *  et spécifiée en paramètre.
 *  Dans cette nouvelle version, il n'y a pas de gestion de mutex. */
template<typename T>
struct TamponNv2: Sink<T, entier>
{
  entier N = 0, windex = 0;
  function<void (const Vecteur<T> &)> callback;
  Vecteur<T> tampon;

  TamponNv2(entier N, function<void (const Vecteur<T> &)> callback)
  {
    this->callback = callback;
    Configurable<entier>::configure(N);
  }

  entier configure_impl(const entier &N)
  {
    this->N = N;
    windex  = 0;
    // Allocation faite plus tard, afin d'éviter d'allouer de la mémoire si elle n'est pas utilisée...
    retourne 0;
  }

  inline entier dim()
  {
    retourne windex;
  }

  void step(const Vecteur<T> &x)
  {
    si((windex == 0) && (x.rows() == N))
    {
      si(callback)
        callback(x);
      retourne;
    }

    si((tampon.rows() == 0) && (N > 0))
      tampon.resize(N);

    entier nb_ech_entree = x.rows();
    entier i = 0;
    tantque(i < nb_ech_entree)
    {
      /* Write index + number of bytes remaining in
       * the input datablock. */
      entier new_windex = windex + nb_ech_entree - i;

      /* If too much data pour the FIFO */
      si(new_windex >= N)
        new_windex = N;

      /* Copy input data to the FIFO */
      entier nech = new_windex - windex;

      tampon.segment(windex, nech) = x.segment(i, nech);

      i      += nech;
      windex = new_windex;

      si(windex == N)
      {
        si(callback)
          callback(tampon);
        windex = 0;
      }
    }
  }
};



template<typename T>
  sptr<Sink<T,entier>> tampon_création(entier N,
      function<void (const Vecteur<T> &)> callback)
{
  retourne make_shared<TamponNv2<T>>(N, callback);
}


template sptr<Sink<float,entier>> tampon_création<float>(entier N,
    function<void (const Vecteur<float> &)> callback);

template sptr<Sink<cfloat,entier>> tampon_création<cfloat>(entier N,
    function<void (const Vecteur<cfloat> &)> callback);

string vstrprintf( const char* format, va_list argv)
{
  va_list argv2;
  va_copy(argv2, argv);
  entier length = vsnprintf(nullptr, 0, format, argv);
  tsd_assert(length >= 0);

  soit buf = new char[length + 1];
  vsnprintf(buf, length + 1, format, argv2);

  va_end(argv2);

  string str(buf);
  delete[] buf;
  retourne str;
  //retourne move(str);
}



Vecf randu(entier n, float a, float b)
{
  uniform_real_distribution<float> distribution(a, b);

  Vecf X(n);

  soit ptr = X.data();
  pour(auto i = 0; i < n; i++)
    *ptr++ = distribution(generateur_aleatoire);

  retourne X;
}

float randu()
{
  uniform_real_distribution<float> distribution(-1, 1);
  retourne distribution(generateur_aleatoire);
}

Vecb randb(entier n)
{
  Vecb y(n);
  default_random_engine generator;
  uniform_int_distribution<> dis(0, 1);
  pour(auto j = 0; j < n; j++)
    y(j) = dis(generator);
  retourne y;
}


entier randi(entier M)
{
  uniform_int_distribution<entier> distribution(0, M-1);
  retourne distribution(generateur_aleatoire);
}

Veci randi(entier M, entier n)
{
  uniform_int_distribution<entier> distribution(0, M-1);
  Veci X(n);
  soit ptr = X.raw_data;
  pour(auto i = 0; i < n; i++)
    *ptr++ = distribution(generateur_aleatoire);
  retourne X;
}

float randn()
{
  normal_distribution<float> distribution(0.0f, 1.0f);
  retourne distribution(generateur_aleatoire);
}

Vecf randn(entier n)
{
  //static /*default_random_engine*//*mt19937_64*/minstd_rand0 generateur;
  //generateur.seed()

  normal_distribution<float> distribution(0.0f, 1.0f);
  //msg("      randn({})...", n);
  Vecf X(n);
  soit ptr = X.raw_data;
  pour(auto i = 0; i < n; i++)
    *ptr++ = distribution(generateur_aleatoire);
  //msg("      ok.");
  retourne X;
}


Veccf randcn(entier n)
{
  normal_distribution<float> distribution(0.0f, 1.0f);
  Veccf X(n);
  soit ptr = X.raw_data;
  pour(auto i = 0; i < n; i++)
    *ptr++ = {distribution(generateur_aleatoire), distribution(generateur_aleatoire)};
  retourne X;
}


Tabf randn_2d(unsigned int n, unsigned int m)
{
  //static default_random_engine generateur;
  normal_distribution<float> distribution(0.0f, 1.0f);

  Tabf X(n, m);
  soit ptr = X.raw_data;
  pour(auto i = 0u; i < m*n; i++)
    *ptr++ = distribution(generateur_aleatoire);
  retourne X;
}

template<typename T>
Vecteur<T> déplie_phase(const Vecteur<T> &x, float r)
{
  soit n = x.rows();
  Vecteur<T> y(n);

  y(0) = x(0);
  pour(auto i = 1; i < n; i++)
  {
    T d = fmod(x(i)-y(i-1), r);

    // d est entre -2pi et 2pi
    si(d > r/2)
      d -= r;
    si(d < -π)
      d += r;
    y(i) = y(i-1) + d;
  }

  retourne y;
}


vector<entier> trouve(const Vecb &x)
{
  vector<entier> res;
  pour(auto i = 0; i < x.rows(); i++)
    si(x(i))
      res.push_back(i);
  retourne res;
}

/** @brief retourne l'index du premier élément vrai, ou -1 si aucun n'est trouvé */
entier trouve_premier(const Vecb &x)
{
  pour(auto i = 0; i < x.rows(); i++)
    si(x(i))
      retourne i;
  retourne -1;
}

/** @brief retourne l'index du premier élément vrai, ou -1 si aucun n'est trouvé */
entier trouve_dernier(const Vecb &x)
{
  pour(auto i = x.rows() - 1; i >= 0; i--)
    si(!x(i))
      retourne i + 1;
  retourne -1;
}


struct OLUT::Impl
{
  entier N;
  Veccf lut;
  Impl(entier N)
  {
    this->N = N;
    soit θ = linspace(0, 2*π_f*(1 - 1.0f/N), N);
    lut = polar(θ);
  }
  cfloat step(float θ0)
  {
    soit θ = modulo_2π(θ0);

    // Interpolation linéaire
    tsd_assert((θ >= 0) && (θ < 2 * π_f));
    float idf = θ * N / (2 * π_f);
    entier idi = floor(idf);
    float f = idf - idi;
    tsd_assert((idi >= 0) && (idi < N));
    retourne lut(idi) * (1-f) + lut((idi+1)%N) * f;
  }
};


OLUT::OLUT(entier resolution)
{
  impl = make_shared<Impl>(resolution);
}
cfloat OLUT::step(float θ)
{
  retourne impl->step(θ);
}



struct OHC: Source<cfloat, OHConfig>
{
  float freq = 1;
  cdouble acc = 1,
          rotation = 1;

  OHC(){}

  OHC(const float &f)
  {
    configure({f});
  }

  entier configure_impl(const OHConfig &cfg)
  {
    freq      = cfg.freq;
    rotation  = std::polar(1.0, 2.0 * π * freq);
    retourne 0;
  }

  Veccf step(entier n)
  {
    Veccf y(n);
    // A chaque première itération,
    // corrige le module
    acc /= abs(acc);
    pour(auto i = 0; i < n; i++)
    {
      y(i) = acc;
      acc *= rotation;
    }
    retourne y;
  }

};

struct OHR: Source<float, OHConfig>
{
  OHC ohc;

  OHR(const float &nu)
  {
    configure({nu});
  }

  entier configure_impl(const OHConfig &config)
  {
    retourne ohc.configure(config);
  }

  Vecf step(entier n)
  {
    retourne real(ohc.step(n));
  }
};



sptr<Source<float, OHConfig>> source_ohr(float freq)
{
  retourne make_shared<OHR>(freq);
}

sptr<Source<cfloat, OHConfig>> source_ohc(float freq)
{
  retourne make_shared<OHC>(freq);
}


namespace hidden{
soit déplie_phase1 = déplie_phase<float>;
soit déplie_phase2 = déplie_phase<double>;
}


}
