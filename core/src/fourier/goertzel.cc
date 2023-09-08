#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"


namespace tsd::fourier {

float goertzel(const Vecf &x, float f)
{
  soit n   = x.rows();
  soit c   = cos(2 * π * f);
  float w0 = 0, w1 = 0,
        en = square(x).somme();

  // Dénominateur
  pour(auto i = 0; i < n; i++)
    std::tie(w0, w1) = std::pair(2 * c * w0 - w1  + x(i), w0);

  // i = 0   : w0 = x0, w1 = 0
  // i = 1   : w0 = ..., w1 = x0
  // ....
  // i = N-1 : ...

  // Dernière étape, avec xn = 0 :
  // s_N = ...
  std::tie(w0, w1) = std::pair(2 * c * w0 - w1, w0);

  // Normalisation par rapport à l'énergie du signal
  retourne 2 * (w0 * w0 - 2 * c * w0 * w1 + w1 * w1) / (en * n);
}


struct Goertzel: FiltreGen<float>
{
  float w0 = 0, w1 = 0, c = 0;
  entier  cnt = 0,
          R   = 0; // Nombre de points de la FFT équivalente

  sptr<FiltreGen<float, float>> filtre_energie;

  /** Initialisation contexte de goertzel */
  Goertzel(float frequence, entier R)
  {
    c               = cos(2.0f * π * frequence);
    this->R         = R;
    filtre_energie  = tsd::filtrage::filtre_mg<float,double>(R);
  }


  void step(const Vecteur<float> &x, Vecteur<float> &y)
  {
    soit n = x.rows(), idx = 0;
    y.resize((n + cnt) / R);

    soit en = filtre_energie->step(square(x));

    pour(auto i = 0; i < n; i++)
    {
      std::tie(w0, w1) = std::pair(2 * c * w0 - w1  + x(i), w0);
      cnt++;

      si(cnt >= R)
      {
        // Dernière étape, avec x_N = 0
        soit [w0p, w1p] = std::pair(2 * c * w0 - w1, w0);

        assertion(idx < y.rows());
        // Normalisation par rapport à l'énergie du signal
        y(idx++) = 2 * (w0p * w0p - 2 * c * w0p * w1p + w1p * w1p) / (R * en(i) * R);

        // Redémmarage
        cnt = 0;
        w0 = w1 = 0;
      }
    }
    assertion(idx == y.rows());
  }
};

sptr<FiltreGen<float>> filtre_goertzel(float frequence, entier N)
{
  retourne make_shared<Goertzel>(frequence, N);
}

}

