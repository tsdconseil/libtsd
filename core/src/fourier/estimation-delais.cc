#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"

using namespace std;

// Interpolation quadratique,
// d'après https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
static float qint_loc(float ym1, float y0, float yp1)
{
  return (yp1 - ym1) / (2 * (2*y0 - yp1 - ym1));
}


namespace tsd::fourier {


/** @brief Estimation du délais le
 *  plus probable entre deux signaux (via une corrélation linéaire) */
tuple<float,float> estimation_délais(IArrayXcf x, IArrayXcf y)//, bool fenêtrage)
{
  auto [xp, yp] = pad_zeros(x, y);

  int N = xp.rows();

  tsd_assert(xp.rows() == N);
  tsd_assert(yp.rows() == N);

  /*if(fenêtrage)
  {
    ArrayXf f = tsd::filtrage::fenetre("hann", N, false);
    xp *= f;
    yp *= f;
  }*/


  // Utilise la version biaisée car sinon amplification des délais importants
  auto [lags, corr] = xcorrb(xp, yp);

  // Score normalisé
  auto e1 = sqrt(xp.abs2().mean());
  auto e2 = sqrt(yp.abs2().mean());
  ArrayXf cn = corr.abs() / (e1 * e2 + 1e-50);

  int index;
  float score = cn.maxCoeff(&index);

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

  float δ = 0;

  if((index > 0) && (index + 1 < cn.rows()))
  {
    float cn0 = cn(index-1);
    float cn1 = cn(index+1);

    // sortie du corrélateur = corrélation non normalisée
    δ = qint_loc(cn0, score, cn1);

    msg("cn0 = {}, cn = {}, cn1 = {} -> δ = {}", cn0, score, cn1, δ);

    if((δ < -0.5) || (δ > 0.5))
    {
      msg_avert("Echec interpolation quadratique: pics = {} {} {}, δ={}, idx={}, N={}",
          cn0, score, cn1, δ, index, N);
      δ = clamp(δ, -0.5f, 0.5f);
    }
  }

  //msg("index = {}, lag = {}, dim = {}, lagf = {}", index, lags(index), xp.rows(), lags(index) + δ);

  return {lags(index) + δ, score};



  /*auto [xp, yp] = pad_zeros(x, y);
  EstimateurDelais est;
  est.configure(xp.rows());
  float délais, score;
  est.step(xp, yp, délais, score);
  return {délais, score};*/
}

int estimation_délais_entier(IArrayXcf x, IArrayXcf y, float &score)
{
  auto [d, s] = estimation_délais(x, y);
  score = s;
  return (int) round(d);
}


template<typename T>
  tuple<Vecteur<T>, Vecteur<T>, int, float> aligne_entier(const Vecteur<T> &x1, const Vecteur<T> &y1)
{
  float score;

  ArrayXcf x = x1, y = y1;

  // Zéro-padding pour avoir la même dimension et qu'elle soit une puissance de 2
  auto [xp, yp] = pad_zeros(x, y, true);


  //infos("aligne entier : #x=%d, #y=%d, #xp=%d, #yp=%d.", (int) x.rows(), (int) y.rows(), (int) xp.rows(), (int) yp.rows());

  int d = estimation_délais_entier(xp, yp, score);

  ArrayXcf xa, ya;
  if(d == 0)
  {
    xa = x1;
    ya = y1;
    //xa = xp;
    //ya = yp;
  }
  else if(d < 0)
  {
    xa = x1.tail(x1.rows() + d);
    ya = y1;

    //xa = x1.tail(x1.rows()+d);
    //ya = y1;//.head(x1.rows());
    //xa = xp.tail(xp.rows()+d);
    //ya = yp.head(xa.rows());
  }
  else
  {
    xa = x1;
    ya = y1.tail(y1.rows() - d);
    //xa = y1.tail(y1.rows()-d);
    //ya = x1.head(ya.rows());
    //ya = yp.tail(yp.rows()-d);
    //xa = xp.head(ya.rows());
  }

  if(xa.rows() > ya.rows())
    xa = xa.head(ya.rows()).eval();
  else if(xa.rows() < ya.rows())
    ya = ya.head(xa.rows()).eval();

  //infos("aligne entier : d = %d, xa.rows() = %d, ya.rows() = %d", d, (int) xa.rows(), (int) ya.rows());

  tsd_assert(xa.rows() == ya.rows());

  if constexpr(est_complexe<T>())
    return {xa, ya, d, score};
  else
    return {xa.real(), ya.real(), d, score};
}

template
  tuple<Vecteur<float>, Vecteur<float>, int, float>
    aligne_entier<float>(const Vecteur<float> &x1, const Vecteur<float> &y1);

template
  tuple<Vecteur<cfloat>, Vecteur<cfloat>, int, float>
    aligne_entier<cfloat>(const Vecteur<cfloat> &x1, const Vecteur<cfloat> &y1);


}

