#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/vue.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/tests.hpp"
#include "tsd/divers.hpp"

#include <Eigen/Eigenvalues>


using namespace tsd::vue;

namespace tsd::filtrage
{

// Renvoie un intervalle entre -1/2 et 1/2 de n points
// avec symétrie ou périodicité
// Si intervalle périodique, va entre
// -1/2 et 1/2 - (n-1)/n² ???
ArrayXf fen_inter(int n, bool sym)
{
  float tmin, tmax;

  if((n & 1) == 0)
  {
    if(sym)
    {
      // -2, . , ., 2
      tmin = -n/2;
      tmax = n/2;
    }
    else
    {
      // Périodique :
      // [-2, -1, 0, 1] 2
      tmin = -n/2;
      tmax = (n-1)/2;
    }
  }
  // IMPAIR
  else
  {
    if(sym)
    {
      // -1, 0, 1
      tmin = -n/2;
      tmax = n/2;
    }
    else
    {
      // Ex : n = 3, n/2 = 1
      // Périodique :
      // [-1 . .] 1
      // [-2 . . . .] 2
      tmin = -n/2;
      tmax = n/2 - ((float) n - 1) / n;
      // -> 1/3, ok
    }
  }

  tmin *= 1.0f / n;
  tmax *= 1.0f / n;

  return linspace(tmin, tmax, n);
}

ArrayXf fenêtre_kaiser1(int n, float β, bool symetrique)
{
  β *= π;

  auto d = 1.0 / std::cyl_bessel_i(0, β);
  ArrayXf x(n);

  ArrayXf t = fen_inter(n, symetrique);

  for(auto i = 0; i < n; i++)
  {
    // -1 ... 1
    //auto K = (2.0 * i / (n-1)) - 1.0;

    auto K = 2 * t(i);

    // -1  : bessel(0, 0)
    // 1   : idem
    // 1-ε : sqrt(2ε)

    x(i) = std::cyl_bessel_i(0, β * sqrt(1.0 - (K * K))) * d;

    // auto k = (2 * beta / (n-1)) * std::sqrt(i * (n-1-i));
    //x(i) = std::cyl_bessel_i(0, k) * d;

  }

  //msg("kaiser : {}", x);

  return x;
}

std::tuple<float, int> kaiser_param(float atten_db, float df)
{
  float domega = 2*π*df;
  int n = (int) ceil((atten_db - 7.95) / (2.285 * domega));
  if(n <= 0)
    n = 1;
  if((n & 1) == 0)
    n++;
  float β = 0;
  if(atten_db > 50)
    β = 0.1102 * (atten_db - 8.7);
  else if(atten_db >= 21)
    β = 0.5842 * std::pow(atten_db - 21, 0.4) + 0.07886 * (atten_db - 21);
  msg("Fen kaiser spec : atten={:.1f} dB, df={} => n={}, β={}", atten_db, df, n, β);
  return {β / π, n};
}

ArrayXf fenêtre_kaiser(float atten_db, float δf, bool symetrique)
{
  auto [β, n] = kaiser_param(atten_db, δf);
  return fenêtre_kaiser1(n, β, symetrique);
}

// https://ccrma.stanford.edu/%7Ejos/sasp/Hann_Hanning_Raised_Cosine.html
// - Si n est pair :
//  - Fenêtre périodique :
//     -n/2 -> (n-1)/2
//  - Fenêtre symétrique :
//
// (plus facile de faire une fenêtre périodique dans ce cas là !!!)


// Si n est impair :
// -n/2 -> n/2 = (n-1)/2
// -pi -> pi
//  - Symétrie ? oui




static ArrayXf Hamming_generalise(float a, int n, bool sym)
{
  //ArrayXf t = 2 * π * fen_inter(n, sym);
  //ArrayXf r2 = a + (1-a) * t.cos();
  //return r2;
  return a + (1-a) * (2 * π * fen_inter(n, sym)).cos();
}

std::string Fenetre2string(Fenetre f)
{
  switch(f)
  {
  case Fenetre::AUCUNE:
    return "re";
  case Fenetre::HANN:
    return "hn";
  case Fenetre::TRIANGLE:
    return "tr";
  case Fenetre::HAMMING:
    return "hm";
  case Fenetre::BLACKMAN:
    return "bm";
  case Fenetre::CHEBYCHEV:
    return "ch";
  }
  return "?";
}

std::ostream& operator<<(std::ostream& ss, const Fenetre &t)
{
  switch(t)
  {
  case Fenetre::AUCUNE:
    ss << "rectangulaire";
    break;
  case Fenetre::HANN:
    ss << "Hann";
    break;
  case Fenetre::TRIANGLE:
    ss << "Triangulaire";
    break;
  case Fenetre::HAMMING:
    ss << "Hamming";
    break;
  case Fenetre::BLACKMAN:
    ss << "Blackman";
    break;
  case Fenetre::CHEBYCHEV:
    ss << "Chebychev";
    break;
  default:
    ss << "?";
  }
  return ss;
}

static Fenetre parse_fenêtre(const std::string &nom)
{
  struct FenId
  {
    std::string nom;
    Fenetre type;
  };
  FenId lst[] =
  {
    {"",        Fenetre::AUCUNE},
    {"aucune",  Fenetre::AUCUNE},
    {"none",    Fenetre::AUCUNE},
    {"re",      Fenetre::AUCUNE},
    {"hn",      Fenetre::HANN},
    {"hann",    Fenetre::HANN},
    {"tr",      Fenetre::TRIANGLE},
    {"hamming", Fenetre::HAMMING},
    {"hm",      Fenetre::HAMMING},
    {"cheb",    Fenetre::CHEBYCHEV},
    {"ch",      Fenetre::CHEBYCHEV},
    {"bm",      Fenetre::BLACKMAN}
  };
  for(auto &fi: lst)
    if(fi.nom == nom)
      return fi.type;
  msg_erreur("Type de fenêtre non reconnu : {}", nom);
  return Fenetre::AUCUNE;
}

ArrayXf fenetre(const std::string &type, int n, bool symetrique)
{
  return fenetre(parse_fenêtre(type), n, symetrique);
}

ArrayXf fenetre(Fenetre type, int n, bool symétrique)
{
  ArrayXf x = ArrayXf::Zero(n);

  //if(((n & 1) == 0) && symétrique)
    //msg_avert("Demande de création d'une fenêtre symétrique avec un nombre pair de coefficients (n = {})", n);

  if(type == Fenetre::AUCUNE)
    x.setOnes();
  else if(type == Fenetre::HANN)
  {
    x = Hamming_generalise(0.5, n, symétrique);
  }
  else if(type == Fenetre::HAMMING)
  {
    float a = 0.54;
    //auto a = 0.53836f; // TO CHECK
    x = Hamming_generalise(a, n, symétrique);
  }
  else if(type == Fenetre::TRIANGLE)
  {
    ArrayXf t = fen_inter(n, symétrique);

    // -0.5 -> 0
    // 0    -> 1
    // 0.5  -> 0

    for(auto i = 0; i < n; i++)
    {
      if(t(i) < 0)
        x(i) = 2 * (0.5 + t(i));
      else
        x(i) = 2 * (0.5 - t(i));
    }

    // Laisse un zéro à la fin
    /*if((n & 1) == 0)
    {
      x.head(n/2) = linspace(0, 1, n/2);
      x.segment(n/2-1,n/2) = linspace(1, 0, n/2);
    }
    else
    {
      x.head((n+1)/2) = linspace(0, 1, (n+1)/2);
      x.segment((n+1)/2 - 1,(n+1)/2) = linspace(1, 0, (n+1)/2);
    }*/
  }
  else if(type == Fenetre::BLACKMAN)
  {
    // https://en.wikipedia.org/wiki/Window_function#Blackman_window
    ArrayXf t = 2 * π * fen_inter(n, symétrique) + π;
    auto α = 0.16f;
    auto a0 = (1 - α)/2;
    auto a1 = 0.5f;
    auto a2 = α/2;
    x = a0 - a1 * t.cos() + a2 * (2 * t).cos();
  }
  else if(type == Fenetre::CHEBYCHEV)
  {
    // TODO: remplacer l'énum fenêtre par une vrai classe
    x = fenêtre_chebychev(n, 80);
  }
  else
  {
    echec("fenêtre non gérée ({}).", (int) type);
  }

  return x;
}

auto cheb_eval(auto x, int n)
{
  if((x >= -1) && (x <= 1))
    return std::cos(n * std::acos(x));
  else if(x >= 1)
    return std::cosh(n * std::acosh(x));

  return /*std::pow(-1.0, n)*/ ((n & 1) ? -1 : 1) * std::cosh(n * std::acosh(-x));
}


// D'après https://www.dsprelated.com/showarticle/42.php
// TODO : gestion symétrie, parité
ArrayXf fenêtre_chebychev(int n, float atten_db, bool symetrique)
{
  auto Alin  = std::pow(10.0, atten_db/20.0);
  auto β = real(std::cosh(std::acosh((cdouble) Alin) * (1.0/(n-1))));

  //msg("atten_db={}, beta={}", atten_db, beta);

  // Commençons par faire une fenêtre symétrique
  // Avec 9 points

  ArrayXf A(n-1);
  ArrayXf X(n-1);

  // TODO : remplacer -1^n par un simple fftshift

  for(auto i = 0; i < n - 1; i++)
  {
    double x = β * std::cos((π * i)/(n-1));
    A(i) = x;
    X(i) = /*std::pow(-1, i) **/ ((i & 1) ? -1 : 1) * cheb_eval(x, n-1);
  }


  ArrayXf w(n);
  w.head(n-1) = real(tsd::fourier::ifft(X));

  w(0) /= 2;
  w(n-1) = w(0);

  w /= w.maxCoeff();

  msg("Fen cheby({}): w(0) = {}, w(n-1) = {}", n, w(0), w(n-1));

  return w;
}


MatrixXf slepian_matrice(int N, float B)
{
  MatrixXf X(N, N);

  // TODO : la moitié des coefs pourrait être implicites (symétrie)
  for(auto i = 0; i < N; i++)
  {
    for(auto j = 0; j < N; j++)
    {
      // sin (π * B * (i-j)) / (π * B * (i-j)
      X(i,j) = π * B * sinc((i-j) * B);
    }
  }
  return X;
}

MatrixXf slepian_evec(int N, float B)
{
  MatrixXf X = slepian_matrice(N, B);
  Eigen::SelfAdjointEigenSolver<MatrixXf> eig;//(X);
  eig.compute(X);
  return eig.eigenvectors();
}

ArrayXf fenêtre_slepian(int N, float B)
{
  return slepian_evec(N,B).col(N-1).array();
}


static int premier_max_local(const ArrayXf &x)
{
  for(auto i = 1; i + 1 < x.rows(); i++)
  {
    if((x(i) > x(i-1)) && (x(i) > x(i+1)))
      return i;
  }
  return -1;
}

/*
static int premier_min_local(const ArrayXf &x)
{
  return premier_max_local(-x);
}*/


FenInfos filtre_pb_analyse(const ArrayXf &h)
{
  auto [fr, mag] = frmag(h);
  return filtre_pb_analyse(h.rows(), fr, mag);
}

FenInfos filtre_pb_analyse(int ncoefs, const ArrayXf &fr, const ArrayXf &mag, bool do_plot)
{
  std::string nom = "";
  FenInfos res;

  res.symetrique = res.periodique = false;

  {
    ArrayXf Hl = 20*log10(mag);
    // Recherche deuxième lobe
    int i = premier_max_local(Hl);

    auto lst = trouve(Hl < -3);

    int j = -1;//mag.rows();

    if(!lst.empty())
      j = lst[0];


    Figure f;

    if(do_plot)
    {
      f = res.fig.subplot(211);
      f.plot(fr, mag, "b-");
      if(i >= 0)
        f.plot(fr(i), mag(i), "rs");
      f.titre("Vue linéaire");
    }

    res.largeur_lp = j / (1.0f * fr.rows());

    int i2 = i;

    if(i >= 0)
    {
      res.atten_ls  = -Hl.tail(Hl.rows() - i - 1).maxCoeff(&i2);
      //res.atten_ls  = -Hl.segment(i, Hl.rows() - i - 1).maxCoeff(&i2);
      i2 += i;
      res.atten_pls = -Hl(i);
    }

    if(do_plot)
    {
      if(j >= 0)
        f.plot(fr(j), mag(j), "gs");

      f.canva().set_couleur(tsd::vue::Couleur::Bleu);
      f.canva().set_dim_fonte(0.6);
      f.canva().texte({0.1f, 0.9f},
              fmt::format("Largeur lobe principal : {:.5f} (={:.2f}/N)\nAttén. premier lobe sec. : {:.1f} dB\nAttén. pire lob sec. : {:.1f} dB",
              res.largeur_lp, res.largeur_lp * ncoefs, res.atten_pls, res.atten_ls),
              {0.4f, 0.2f});

      f = res.fig.subplot(212);

      f.plot(fr, Hl, "b-");
      f.titre("Vue logarithmique");
      if(i >= 0)
        f.plot(fr(i), Hl(i), "rs");
      if(i2 >= 0)
        f.plot(fr(i2), Hl(i2), "ms");
    }

    msg("  Atténuation lobe principal : {:.1f} dB.", -Hl(0));
    msg("  Largeur lobe principal : {:.3f} (j={}, = {}/N)", res.largeur_lp, j, res.largeur_lp * ncoefs);
    msg("  Atténuation premier lobe secondaire : {:.1f} dB.", res.atten_pls);
    msg("  Atténuation pire lobe secondaire :    {:.1f} dB.", res.atten_ls);
  }
  return res;
}



FenInfos fenetre_analyse(const std::string &nom, const ArrayXf &x, bool do_plot)
{
  int n = x.rows();
  FRat<float> frat(x / x.sum());
  auto [fr,xm] = frmag(frat, 4096);


  auto res = filtre_pb_analyse(x.rows(), fr, xm, do_plot);

  {
    res.symetrique = true;

    int imax = (n-1)/2;

    float emax = 0, etot = 0;
    for(auto i = 0; i <= imax; i++)
    {
      auto err = std::abs(x(i)-x(n-i-1));
      if(err > 1e-6)
        res.symetrique = false;
      emax = std::max(err, emax);
      etot += err;
    }
    msg("  Delta fenêtre symétrique : {} (max {})", etot, emax);
  }

  return res;
}



}
