#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"

using namespace std;

// Interpolation quadratique,
// d'après https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
static float qint_loc(float ym1, float y0, float yp1)
{
  retourne (yp1 - ym1) / (2 * (2*y0 - yp1 - ym1));
}


namespace tsd::fourier {


/** @brief Estimation du délais le
 *  plus probable entre deux signaux (via une corrélation linéaire) */
tuple<float,float> estimation_délais(const Veccf &x, const Veccf &y)
{
  soit [xp, yp] = pad_zeros(x, y);

  soit N = xp.rows();

  assertion((xp.rows() == N) && (yp.rows() == N));

  /*si(fenêtrage)
  {
    soit f = tsd::filtrage::fenetre("hann", N, non);
    xp *= f;
    yp *= f;
  }*/


  // Utilise la version biaisée car sinon amplification des délais importants
  soit [lags, corr] = xcorrb(xp, yp);

  // Score normalisé
  soit e1 = sqrt(abs2(xp).moyenne()),
       e2 = sqrt(abs2(yp).moyenne());
  soit cn = abs(corr) / (e1 * e2 + 1e-50);

  soit [score, index] = cn.max();

  /*{
    tsd::vue::Figure f;
    f.plot(xp, "b-", "xp");
    f.plot(yp, "r-", "yp");
    f.enregistrer("./essai0.png");
  }

  {
    tsd::vue::Figure f;
    f.plot(lags.segment(index-10, 20), cn.segment(index - 10, 20), "o-g");
    f.plot(lags(index), cn(index), "or");
    f.enregistrer("./essai0-zoom.png");
  }

  {
    tsd::vue::Figure f;
    f.plot(lags, cn);
    f.plot(lags(index), cn(index), "or");
    f.enregistrer("./essai.png");

    msg("cn(index - 1) = {}", cn(index-1));
    msg("cn(index) = {}", cn(index));
    msg("cn(index + 1) = {}", cn(index+1));
  }*/

  soit δ = 0.0f;

  si((index > 0) && (index + 1 < cn.rows()))
  {
    soit cn0 = cn(index-1),
         cn1 = cn(index+1);

    // sortie du corrélateur = corrélation non normalisée
    δ = qint_loc(cn0, score, cn1);

    msg("cn0 = {}, cn = {}, cn1 = {} -> δ = {}", cn0, score, cn1, δ);

    si((δ < -0.5) || (δ > 0.5))
    {
      msg_avert("Echec interpolation quadratique: pics = {} {} {}, δ={}, idx={}, N={}",
          cn0, score, cn1, δ, index, N);
      δ = clamp(δ, -0.5f, 0.5f);
    }
  }

  //msg("index = {}, lag = {}, dim = {}, lagf = {}", index, lags(index), xp.rows(), lags(index) + δ);
  retourne {lags(index) + δ, score};

  /*soit [xp, yp] = pad_zeros(x, y);
  EstimateurDelais est;
  est.configure(xp.rows());
  float délais, score;
  est.step(xp, yp, délais, score);
  retourne {délais, score};*/
}

entier estimation_délais_entier(const Veccf &x, const Veccf &y, float &score)
{
  soit [d, s] = estimation_délais(x, y);
  score = s;
  retourne (entier) round(d);
}


template<typename T>
  tuple<Vecteur<T>, Vecteur<T>, entier, float> aligne_entier(const Vecteur<T> &x1, const Vecteur<T> &y1)
{
  float score;

  soit x = x1.as_complex(),
       y = y1.as_complex();


  // Zéro-padding pour avoir la même dimension et qu'elle soit une puissance de 2
  soit [xp, yp] = pad_zeros(x, y, oui);

  //infos("aligne entier : #x=%d, #y=%d, #xp=%d, #yp=%d.", (entier) x.rows(), (entier) y.rows(), (entier) xp.rows(), (entier) yp.rows());

  soit d = estimation_délais_entier(xp, yp, score);

  Veccf xa, ya;
  si(d == 0)
  {
    xa = x;
    ya = y;
  }
  sinon si(d < 0)
  {
    xa = x.tail(x.rows() + d);
    ya = y;
  }
  sinon
  {
    xa = x;
    ya = y.tail(y.rows() - d);
  }

  si(xa.rows() > ya.rows())
  {
    // Attention: xa peut être une référence ici !

    //msg("Changement dim...");
    xa = xa.head(ya.rows());
  }
  sinon si(xa.rows() < ya.rows())
  {
    ya = ya.head(xa.rows());
  }

  assertion(xa.rows() == ya.rows());

  si constexpr(est_complexe<T>())
    retourne {xa, ya, d, score};
  sinon
    retourne {real(xa), real(ya), d, score};
}

template
  tuple<Vecteur<float>, Vecteur<float>, entier, float>
    aligne_entier<float>(const Vecteur<float> &x1, const Vecteur<float> &y1);

template
  tuple<Vecteur<cfloat>, Vecteur<cfloat>, entier, float>
    aligne_entier<cfloat>(const Vecteur<cfloat> &x1, const Vecteur<cfloat> &y1);


}

