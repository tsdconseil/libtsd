#include "tsd/tsd-all.hpp"
#include "tsd/tests.hpp"
#include "tsd/divers.hpp"

#include "tsd/eig-util.hpp"

#include <Eigen/Eigenvalues>


namespace tsd::filtrage
{

// Renvoie un intervalle entre -1/2 et 1/2 de n points
// avec symétrie ou périodicité
// si intervalle périodique, va entre
// -1/2 et 1/2 - (n-1)/n² ???
Vecf fen_inter(entier n, bouléen sym)
{
  float tmin, tmax;

  si((n & 1) == 0)
  {
    si(sym)
    {
      // -2, . , ., 2
      tmin = -n/2;
      tmax = n/2;
    }
    sinon
    {
      // Périodique :
      // [-2, -1, 0, 1] 2
      tmin = -n/2;
      tmax = (n-1)/2;
    }
  }
  // IMPAIR
  sinon
  {
    si(sym)
    {
      // -1, 0, 1
      tmin = -n/2;
      tmax = n/2;
    }
    sinon
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

  retourne linspace(tmin / n, tmax / n, n);
}

Vecf fenêtre_kaiser1(entier n, float β, bouléen symetrique)
{
  β *= π;

  soit d = 1.0 / std::cyl_bessel_i(0, β);
  Vecf x(n);
  soit t = fen_inter(n, symetrique);

  pour(auto i = 0; i < n; i++)
  {
    // -1 ... 1
    //soit K = (2.0 * i / (n-1)) - 1.0;

    soit K = 2 * t(i);

    // -1  : bessel(0, 0)
    // 1   : idem
    // 1-ε : sqrt(2ε)

    x(i) = std::cyl_bessel_i(0, β * sqrt(1.0 - (K * K))) * d;
  }
  retourne x;
}

tuple<float, entier> kaiser_param(float atten_db, float df)
{
  soit dω = 2*π_f*df;
  soit n = (entier) ceil((atten_db - 7.95) / (2.285 * dω));
  si(n <= 0)
    n = 1;
  si((n & 1) == 0)
    n++;
  soit β = 0.0f;
  si(atten_db > 50)
    β = 0.1102 * (atten_db - 8.7);
  sinon si(atten_db >= 21)
    β = 0.5842 * pow(atten_db - 21, 0.4) + 0.07886 * (atten_db - 21);
  msg("Fen kaiser spec : atten={:.1f} dB, df={} => n={}, β={}", atten_db, df, n, β);
  retourne {β / π, n};
}

Vecf fenêtre_kaiser(float atten_db, float δf, bouléen symetrique)
{
  soit [β, n] = kaiser_param(atten_db, δf);
  retourne fenêtre_kaiser1(n, β, symetrique);
}

// https://ccrma.stanford.edu/%7Ejos/sasp/Hann_Hanning_Raised_Cosine.html
// - si n est pair :
//  - Fenêtre périodique :
//     -n/2 -> (n-1)/2
//  - Fenêtre symétrique :
//
// (plus facile de faire une fenêtre périodique dans ce cas là !!!)


// si n est impair :
// -n/2 -> n/2 = (n-1)/2
// -pi -> pi
//  - Symétrie ? oui




static Vecf Hamming_generalise(float a, entier n, bouléen sym)
{
  retourne a + (1-a) * cos(2 * π * fen_inter(n, sym));
}

string Fenetre2string(Fenetre f)
{
  switch(f)
  {
    case Fenetre::AUCUNE:
      retourne "re";
    case Fenetre::HANN:
      retourne "hn";
    case Fenetre::TRIANGLE:
      retourne "tr";
    case Fenetre::HAMMING:
      retourne "hm";
    case Fenetre::BLACKMAN:
      retourne "bm";
    case Fenetre::CHEBYCHEV:
      retourne "ch";
  }
  retourne "?";
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
  retourne ss;
}

static Fenetre parse_fenêtre(cstring nom)
{
  struct FenId
  {
    string nom;
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
  pour(auto &fi: lst)
    si(fi.nom == nom)
      retourne fi.type;
  msg_erreur("Type de fenêtre non reconnu : {}", nom);
  retourne Fenetre::AUCUNE;
}

Vecf fenetre(cstring type, entier n, bouléen symetrique)
{
  retourne fenetre(parse_fenêtre(type), n, symetrique);
}

Vecf fenetre(Fenetre type, entier n, bouléen symétrique)
{
  soit x = Vecf::zeros(n);

  //si(((n & 1) == 0) && symétrique)
    //msg_avert("Demande de création d'une fenêtre symétrique avec un nombre pair de coefficients (n = {})", n);

  si(type == Fenetre::AUCUNE)
    x.setConstant(1);
  sinon si(type == Fenetre::HANN)
    x = Hamming_generalise(0.5, n, symétrique);
  sinon si(type == Fenetre::HAMMING)
  {
    soit a = 0.54f;
    //soit a = 0.53836f; // TO CHECK
    x = Hamming_generalise(a, n, symétrique);
  }
  sinon si(type == Fenetre::TRIANGLE)
  {
    soit t = fen_inter(n, symétrique);

    // -0.5 -> 0
    // 0    -> 1
    // 0.5  -> 0

    pour(auto i = 0; i < n; i++)
    {
      si(t(i) < 0)
        x(i) = 2 * (0.5 + t(i));
      sinon
        x(i) = 2 * (0.5 - t(i));
    }
  }
  sinon si(type == Fenetre::BLACKMAN)
  {
    // https://en.wikipedia.org/wiki/Window_function#Blackman_window
    soit t  = 2 * π * fen_inter(n, symétrique) + π;
    soit α  = 0.16f,
         a0 = (1 - α)/2,
         a1 = 0.5f,
         a2 = α/2;
    x = - a1 * cos(t) + a2 * cos(2 * t) + a0;
  }
  sinon si(type == Fenetre::CHEBYCHEV)
  {
    // TODO: remplacer l'énum fenêtre par une vrai classe
    x = fenêtre_chebychev(n, 80);
  }
  sinon
  {
    échec("fenêtre non gérée ({}).", (entier) type);
  }

  retourne x;
}

auto cheb_eval(auto x, entier n)
{
  si((x >= -1) && (x <= 1))
    retourne cos(n * std::acos(x));
  sinon si(x >= 1)
    retourne std::cosh(n * std::acosh(x));

  retourne /*pow(-1.0, n)*/ ((n & 1) ? -1 : 1) * std::cosh(n * std::acosh(-x));
}


// D'après https://www.dsprelated.com/showarticle/42.php
// TODO : gestion symétrie, parité
Vecf fenêtre_chebychev(entier n, float atten_db, bouléen symetrique)
{
  si(n <= 1)
    retourne Vecf::ones(n);

  soit Alin = pow(10.0, atten_db/20.0),
          β = real(std::cosh(std::acosh((cdouble) Alin) * (1.0/(n-1))));

  Vecf A(n-1), X(n-1);

  // TODO : remplacer -1^n par un simple fftshift

  pour(auto i = 0; i < n - 1; i++)
  {
    soit x = β * std::cos((π * i)/(n-1));
    A(i) = x;
    X(i) = ((i & 1) ? -1 : 1) * cheb_eval(x, n-1);
  }


  Vecf w(n);
  w.head(n-1) = real(ifft(X));

  w(0) /= 2;
  w(n-1) = w(0);
  w /= w.valeur_max();

  //msg("Fen cheby({}): w(0) = {}, w(n-1) = {}", n, w(0), w(n-1));

  retourne w;
}


Eigen::MatrixXf slepian_matrice(entier N, float B)
{
  Eigen::MatrixXf X(N, N);

  // TODO : la moitié des coefs pourrait être implicites (symétrie)
  pour(auto i = 0; i < N; i++)
  {
    pour(auto j = 0; j < N; j++)
    {
      // sin (π * B * (i-j)) / (π * B * (i-j)
      X(i,j) = π * B * sinc((i-j) * B);
    }
  }
  retourne X;
}

Tabf slepian_evec(entier N, float B)
{
  Eigen::MatrixXf X = slepian_matrice(N, B);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig;//(X);
  eig.compute(X);
  retourne etab2tab(eig.eigenvectors());
}

Vecf fenêtre_slepian(entier N, float B)
{
  retourne slepian_evec(N,B).col(N-1);
}


FenInfos filtre_pb_analyse(const Vecf &h)
{
  soit [fr, mag] = frmag(h);
  retourne filtre_pb_analyse(h.rows(), fr, mag);
}

FenInfos filtre_pb_analyse(entier ncoefs, const Vecf &fr, const Vecf &mag, bouléen do_plot)
{
  string nom;
  FenInfos res;

  res.symetrique = res.periodique = non;

  {
    soit Hl = mag2db(mag);
    // Recherche deuxième lobe
    soit i = trouve_premier_max_local(Hl);
    soit lst = trouve(Hl < -3);
    soit j = -1;

    si(!lst.empty())
      j = lst[0];


    Figure f;

    si(do_plot)
    {
      f = res.fig.subplot(211);
      f.plot(fr, mag, "b-");
      si(i >= 0)
        f.plot(fr(i), mag(i), "rs");
      f.titre("Vue linéaire");
    }

    res.largeur_lp = j / (1.0f * fr.rows());

    soit i2 = i;

    si(i >= 0)
    {
      res.atten_ls  = -Hl.tail(Hl.rows() - i - 1).maxCoeff(&i2);
      i2 += i;
      res.atten_pls = -Hl(i);
    }

    si(do_plot)
    {
      si(j >= 0)
        f.plot(fr(j), mag(j), "gs");

      f.canva().set_couleur(Couleur::Bleu);
      f.canva().set_dim_fonte(0.6);
      f.canva().texte({0.1f, 0.9f},
              sformat("Largeur lobe principal : {:.5f} (={:.2f}/N)\nAttén. premier lobe sec. : {:.1f} dB\nAttén. pire lob sec. : {:.1f} dB",
              res.largeur_lp, res.largeur_lp * ncoefs, res.atten_pls, res.atten_ls),
              {0.4f, 0.2f});

      f = res.fig.subplot(212);

      f.plot(fr, Hl, "b-");
      f.titre("Vue logarithmique");
      si(i >= 0)
        f.plot(fr(i), Hl(i), "rs");
      si(i2 >= 0)
        f.plot(fr(i2), Hl(i2), "ms");
    }

    msg("  Atténuation lobe principal : {:.1f} dB.", -Hl(0));
    msg("  Largeur lobe principal : {:.3f} (j={}, = {}/N)", res.largeur_lp, j, res.largeur_lp * ncoefs);
    msg("  Atténuation premier lobe secondaire : {:.1f} dB.", res.atten_pls);
    msg("  Atténuation pire lobe secondaire :    {:.1f} dB.", res.atten_ls);
  }
  retourne res;
}



FenInfos fenetre_analyse(cstring nom, const Vecf &x, bouléen do_plot)
{
  soit n = x.rows();
  FRat<float> frat(x / x.somme());
  soit [fr,xm] = frmag(frat, 4096);
  soit res = filtre_pb_analyse(n, fr, xm, do_plot);

  {
    res.symetrique = oui;
    soit imax = (n-1)/2;
    soit emax = 0.0f, etot = 0.0f;
    pour(auto i = 0; i <= imax; i++)
    {
      soit err = abs(x(i)-x(n-i-1));
      si(err > 1e-6)
        res.symetrique = non;
      emax = max(err, emax);
      etot += err;
    }
    msg("  Delta fenêtre symétrique : {} (max {})", etot, emax);
  }

  retourne res;
}



}
