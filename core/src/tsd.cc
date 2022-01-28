#include "tsd/tsd.hpp"
#include "tsd/filtrage/frat.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"
#include <random>
#include <cstdarg>
#include <iostream>
#include <fmt/chrono.h>
#include <chrono>
#include <stdexcept>

using namespace std;

namespace tsd {



/** @brief Calcul de la matrice de covariance
 *  @param m Dimension de la matrice (nombre de délais examinés) */
/*template<typename T>
Eigen::MatrixXf cov(const Vecteur<T> &x, int m)
{
  int n = x.rows();
  Eigen::MatrixXf R = Eigen::MatrixXf::Zero(m, m);
  for(auto i = m; i < n; i++)
  {
    auto t = x.segment(i - m, m).reverse().matrix();
    R += t * t.adjoint();
  }
  return R / n;
}*/

ArrayXcf sigexp(float f, int n)
{
  ArrayXcf x(n);
  cdouble r = 1, w0 = std::polar<double>(1.0, f*2*π);

  for(auto i = 0; i < n; i++)
  {
    x(i) = r;
    r *= w0;
  }
  return x;
}

ArrayXf sigscie(int p, int n)
{
  ArrayXf x(n);
  for(auto i = 0; i < n; i++)
    x(i) = ((i % p) - ((p-1) * 0.5)) / (0.5 * (p-1));
  return x;
}

ArrayXf sigtri(int p, int n)
{
  ArrayXf x(n);
  for(auto i = 0; i < n; i++)
  {
    int j = i % (2*p);
    if(j < p)
      x(i) = j;
    else
      x(i) = 2 * p - j;
    x(i) -= p/2;
    x(i) /= p;
  }
 return 2*x; // entre -1 et 1
}

ArrayXf sigimp(int n, int p)
{
  ArrayXf x = ArrayXf::Zero(n);
  x(p) = 1;
  return x;
}

ArrayXf sigcar(int p, int n)
{
  ArrayXf x(n);
  for(auto i = 0; i < n; i++)
    x(i) = 2 * (((i / p) % 2) - 0.5);
  return x;
}

ArrayXf sigsin(float f, int n)
{
  ArrayXf x(n);
  cdouble r = 1, w0 = std::polar<double>(1.0, f*2*π);

  for(auto i = 0; i < n; i++)
  {
    x(i) = r.imag();
    r *= w0;
  }
  return x;
}

ArrayXf sigcos(float f, int n)
{
  ArrayXf x(n);
  cdouble r = 1, w0 = std::polar<double>(1.0, f*2*π);

  for(auto i = 0; i < n; i++)
  {
    x(i) = r.real();
    r *= w0;
  }
  return x;
}

ArrayXf sigchirp(float f0, float f1, int n)
{
  ArrayXf freq  = linspace(f0, f1, n);
  ArrayXf phase = (2 * π) * cumsum(freq);
  return phase.cos();
}

ArrayXf sigchirp2(float f0, float f1, int n)
{
  ArrayXf freq  = f0 + (f1 - f0) * (linspace(0, 1, n)).square();
  ArrayXf phase = (2 * π) * cumsum(freq);
  return phase.cos();
}

ArrayXf siggauss(int n, float a)
{
  ArrayXf t = (linspace(0, n-1, n) - n/2.0f) / (n/2.0f);
  return (- a * t * t).exp();
}

ArrayXf siggsin(float f, int n, float a)
{
  ArrayXf x = sigsin(f, n);
  ArrayXf t = (linspace(0, n-1, n) - n/2.0f) / (n/2.0f);
  return x * (- a * t * t).exp();
}

int prochaine_puissance_de_2(unsigned int i)
{
  int lg2 = (int) ceil(log((float) i) / log(2.0f));
  return 1l << lg2;
}

uint64_t get_tick_count_µs()
{

# if WIN32
  echec("TODO: VSTUDIO / get_tick_count_µs");
  return 0;
# else
  struct timespec ts;
  if(clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
  {
    perror("clock_gettime().");
    return 0;
  }
  return (ts.tv_nsec / 1000) + (((uint64_t)ts.tv_sec) * 1000 * 1000);
# endif
}




// En mode test
bool erreur_attendue = false;


void log_default(const char *fn, const int ligne, int niveau, const char *fonction, const string &str)
{
  fmt::text_style drapeaux;

  if(!erreur_attendue)
  {
    if(niveau == 2)
      drapeaux = fmt::emphasis::bold;
    else if(niveau == 3)
      drapeaux = fg(fmt::color::crimson);
    else if(niveau >= 4)
      drapeaux = fg(fmt::color::crimson) | fmt::emphasis::bold;
  }

  auto tot = get_tick_count_µs();

  static int64_t first = -1, last = -1;
  if(first == -1)
    first = last = tot;

  auto aff = tot - first;
  float diff = tot - last;

  //auto aff = tot - last;


  auto µs  = aff % 1000;
  auto ms  = (aff / 1000) % 1000;
  auto sec = (aff / 1000000);

  string unit = "µs";
  if(diff > 1000)
  {
    diff /= 1000;
    unit = "ms";
  }
  if(diff > 1000)
  {
    diff /= 1000;
    unit = " s";
  }

  auto s1 = fmt::format("[+{: >5.1f} {}]", diff, unit);

  if(unit == " s")
    s1 = "\033[1;32m" + s1 + "\033[0m";
  else if(unit == "ms")
    s1 = "\033[1m" + s1 + "\033[0m";



  auto s_dt  = fmt::format("[{:03d},{:03d},{:03d}] {}", sec, ms, µs, s1);
  auto s_msg = fmt::format(" {}\n", str);

  last = tot;

  if(niveau >= 4)
  {
    s_msg += fmt::format("(fichier : {}, ligne : {})\n", fn, ligne);
  }

  tsd::vue::stdo.printf(s_msg.c_str());
  fmt::print(FMT_RUNTIME(s_dt));
  fmt::print(drapeaux, "{}", s_msg);

  fflush(0);
  flush(cout);

  if(niveau >= 4)
  {
    printf("Throw...\n"); fflush(0);
    throw runtime_error(str);
  }
}





static logger_t log = log_default;


void reset_logger()
{
  log = log_default;
}

void set_logger(logger_t log_)
{
  log = log_;
}

void msg_impl2(const char *fn, const int ligne, int niveau, const char *fonction, const string &str)
{
  if(log)
    log(fn, ligne, niveau, fonction, str);
  else
    // Car l'affectation log=log_default n'est pas forcément déjà faite ici.
    // (par exemple, si appel depuis l'initialisation d'une variable statique)
    log_default(fn, ligne, niveau, fonction, str);
}

/** @brief Cette classe tamponne les données d'entrée
 *  pour les sortir sous la forme de paquets de taille constante
 *  et spécifiée en paramètre.
 *  Dans cette nouvelle version, il n'y a pas de gestion de mutex. */
template<typename T>
struct TamponNv2: Sink<T, int>
{
  int N = 0, windex = 0;
  function<void (const Vecteur<T> &)> callback;
  Vecteur<T> tampon;

  TamponNv2(int N, function<void (const Vecteur<T> &)> callback)
  {
    this->callback = callback;
    Configurable<int>::configure(N);
  }

  int configure_impl(const int &N)
  {
    this->N = N;
    windex  = 0;
    // Allocation faite plus tard, afin d'éviter d'allouer de la mémoire si elle n'est pas utilisée...
    return 0;
  }

  inline int dim()
  {
    return windex;
  }

  void step(const Eigen::Ref<const Vecteur<T>> x)
  {
    if((windex == 0) && (x.rows() == N))
    {
      if(callback)
        callback(x);
      return;
    }

    if((tampon.rows() == 0) && (N > 0))
      tampon.resize(N);

    int nb_ech_entree = x.rows();
    int i = 0;
    while(i < nb_ech_entree)
    {
      /* Write index + number of bytes remaining in
       * the input datablock. */
      int new_windex = windex + nb_ech_entree - i;

      /* If too much data for the FIFO */
      if(new_windex >= N)
        new_windex = N;

      /* Copy input data to the FIFO */
      int nech = new_windex - windex;

      tampon.segment(windex, nech) = x.segment(i, nech);

      i      += nech;
      windex = new_windex;

      if(windex == N)
      {
        if(callback)
          callback(tampon);
        windex = 0;
      }
    }
  }
};



template<typename T>
  sptr<Sink<T,int>> tampon_création(int N,
      function<void (const Vecteur<T> &)> callback)
{
  return make_shared<TamponNv2<T>>(N, callback);
}


template sptr<Sink<float,int>> tampon_création<float>(int N,
    function<void (const Vecteur<float> &)> callback);

template sptr<Sink<cfloat,int>> tampon_création<cfloat>(int N,
    function<void (const Vecteur<cfloat> &)> callback);

string vstrprintf( const char* format, va_list argv)
{
  va_list argv2;
  va_copy(argv2, argv);
  int length = vsnprintf(nullptr, 0, format, argv);
  tsd_assert(length >= 0);

  auto buf = new char[length + 1];
  vsnprintf(buf, length + 1, format, argv2);

  va_end(argv2);

  string str(buf);
  delete[] buf;
  return str;
  //return move(str);
}

ArrayXXf randu_2d(unsigned int n, unsigned int m)
{
  return 0.5 * ArrayXXf::Random(n, m) + 0.5;
}

ArrayXf randu(int n)
{
  return 0.5 * ArrayXf::Random(n) + 0.5;
}

ArrayXf randb(int n)
{
  ArrayXf y(n);
  default_random_engine generator;
  uniform_int_distribution<> dis(0, 1);
  for(auto j = 0; j < n; j++)
    y(j) = dis(generator);
  return y;
}

ArrayXi randi(int M, int n)
{
  static default_random_engine generateur;
  uniform_int_distribution<int> distribution(0, M-1);
  ArrayXi X(n);
  auto ptr = X.data();
  for(auto i = 0; i < n; i++)
    *ptr++ = distribution(generateur);
  return X;
}


ArrayXf randn(int n)
{
  //static /*default_random_engine*//*mt19937_64*/minstd_rand0 generateur;
  static default_random_engine generateur;
  normal_distribution<float> distribution(0.0f, 1.0f);
  //msg("      randn({})...", n);
  ArrayXf X(n);
  auto ptr = X.data();
  for(auto i = 0; i < n; i++)
    *ptr++ = distribution(generateur);
  //msg("      ok.");
  return X;
}



ArrayXXf randn_2d(unsigned int n, unsigned int m)
{
  static default_random_engine generateur;
  normal_distribution<float> distribution(0.0f, 1.0f);

  ArrayXXf X(n, m);
  auto ptr = X.data();
  for(auto i = 0u; i < m*n; i++)
    *ptr++ = distribution(generateur);
  return X;
}

/*Eigen::ArrayXf diff(const Eigen::ArrayXf &x)
{
  auto n = x.rows();
  return x.tail(n-1) - x.head(n-1);
}

Eigen::ArrayXi diff(const Eigen::ArrayXi &x)
{
  auto n = x.rows();
  return x.tail(n-1) - x.head(n-1);
}*/

ArrayXf unwrap(const ArrayXf &x, float r)
{
  auto n = x.rows();
  ArrayXf y(n);

  y(0) = x(0);
  for(auto i = 1; i < n; i++)
  {
    float d = fmod(x(i)-y(i-1), r);

    // d est entre -2pi et 2pi
    if(d > r/2)
      d -= r;
    if(d < -π)
      d += r;
    y(i) = y(i-1) + d;
  }

  return y;
}

/*ArrayXf cumsum(const ArrayXf &x)
{
  auto n = x.rows();
  ArrayXf y(n);
  double cs = 0;
  for(auto i = 0; i < n; i++)
  {
    cs += x(i);
    y(i) = cs;
  }
  return y;
}

Eigen::ArrayXi cumsum(const Eigen::ArrayXi &x)
{
  auto n = x.rows();
  Eigen::ArrayXi y(n);
  int32_t cs = 0;
  for(auto i = 0; i < n; i++)
  {
    cs = cs + x(i);
    y(i) = cs;
  }
  return y;
}*/


vector<int> trouve(IArrayXb x)
{
  vector<int> res;
  for(auto i = 0; i < x.rows(); i++)
    if(x(i))
      res.push_back(i);
  return res;
}

/** @brief Retourne l'index du premier élément vrai, ou -1 si aucun n'est trouvé */
int trouve_premier(IArrayXb x)
{
  for(auto i = 0; i < x.rows(); i++)
    if(x(i))
      return i;
  return -1;
}

ArrayXcf polar(const ArrayXf &ρ, const ArrayXf &θ)
{
  auto n = θ.rows();
  tsd_assert(n == ρ.rows());
  ArrayXcf res(n);
  for(auto i = 0; i < n; i++)
    res(i) = std::polar(ρ(i), θ(i));
  return res;
}

ArrayXcf polar(const ArrayXf &θ)
{
  auto n = θ.rows();
  ArrayXcf res(n);
  for(auto i = 0; i < n; i++)
    res(i) = std::polar(1.0f, θ(i));
  return res;
}




struct OLUT::Impl
{
  int N;
  ArrayXcf lut;
  Impl(int res)
  {
    N = res;
    ArrayXf θ = linspace(0, 2*π_f*(1 - 1.0f/N), N);
    lut = polar(θ);
  }
  cfloat step(float θ0)
  {
    auto θ = wrap_2pi(θ0);

    // Interpolation linéaire
    tsd_assert((θ >= 0) && (θ < 2 * π_f));
    float idf = θ * N / (2 * π_f);
    int idi = floor(idf);
    float f = idf - idi;
    tsd_assert((idi >= 0) && (idi < N));
    return lut(idi) * (1-f) + lut((idi+1)%N) * f;
  }
};


OLUT::OLUT(int resolution)
{
  impl = make_shared<Impl>(resolution);
}
cfloat OLUT::step(float θ)
{
  return impl->step(θ);
}



struct OHC: Source<cfloat, OHConfig>
{
  float freq = 1;
  cfloat phaseur = 1, rotation = 1;

  OHC(){}

  OHC(const float &f)
  {
    configure({f});
  }

  int configure_impl(const OHConfig &cfg)
  {
    freq      = cfg.freq;

    if(cfg.shift)
    {
      rotation *= cfloat(1.0f, 2 * π_f * cfg.df);
    }
    else
      rotation  = std::polar(1.0f, 2.0f * π_f * freq);
    return 0;
  }

  ArrayXcf step(int n)
  {
    ArrayXcf out(n);
    phaseur = phaseur / abs(phaseur); // Periodic update
    for(auto i = 0; i < n; i++)
    {
      out(i) = phaseur;
      phaseur *= rotation;
    }
    return out;
  }

};

struct OHR: Source<float, OHConfig>
{
  OHC ohc;

  OHR(const float &nu)
  {
    configure({nu});
  }

  int configure(const float &nu)
  {
    ohc.configure({nu});
    return 0;
  }

  int configure_impl(const OHConfig &config)
  {
    return configure(config.df);
  }

  ArrayXf step(int n)
  {
    return ohc.step(n).real();
  }
};



sptr<Source<float, OHConfig>> source_ohr(float freq)
{
  return make_shared<OHR>(freq);
}

sptr<Source<cfloat, OHConfig>> source_ohc(float freq)
{
  return make_shared<OHC>(freq);
}




/*ArrayXf resample(IArrayXf x, float lom)
{
  auto f = tsd::filtrage::filtre_reechan<float>(lom);
  return f->step(x);
}
ArrayXcf resample(IArrayXcf x, float lom)
{
  auto f = tsd::filtrage::filtre_reechan<cfloat>(lom);
  return f->step(x);
}*/

/*ArrayXf fenetre(const string &code, unsigned int n)
{
  Fenetre f(code.c_str());
  return f.construire(n);
}

Fenetre::Fenetre()
{

}*/











/*ArrayXf fenetre(Fenetre type, unsigned int n)
{
  Fenetre f(type);
  return f.construire(n);
}*/



}
