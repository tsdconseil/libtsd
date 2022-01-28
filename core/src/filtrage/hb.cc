///////////////////////////////////////
// Conception de filtres demi-bande
///////////////////////////////////////

#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"

using namespace std;
using namespace tsd::vue;


namespace tsd::filtrage {


# if 0
  template<typename T, typename Tc>
  struct FiltreRIFNv: FiltreGen<T>
  {
    Vecteur<T> fenêtre;
    Vecteur<Tc> coefs;
    int index, K;
    int odd = 0, R = 0;

    Vecteur<T> futur;

    FiltreRIFNv(const Vecteur<Tc> &c, unsigned int R)
    {
      odd = 0;
      this->R = R;
      this->coefs = c;
      index = 0;
      K = coefs.rows();
      fenêtre.setZero(K);
    }

    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      tsd_assert(K > 0);
      int n = x.rows();
      auto iptr = x.data();
      auto optr = y.data();
      if(iptr != optr)
      {
        y.resize((n + odd) / R);
        optr = y.data();
      }

      for(auto j = 0; j < n; j++)
      {
        auto cptr = coefs.data();
        auto xi = *iptr++;

        for(auto i = 0; i < n; i++)
          futur((index + i) % K) += xi * *cptr++;

        *optr++ = futur(0);
        index++;
      }
    }
  };
# endif



  // Implémentation optimisée pour filtre demi-bande
  template<typename T, typename Tc>
  struct FiltreRIFDemiBande: FiltreGen<T>
  {
    Vecteur<T> fenêtre;
    Vecteur<Tc> coefs;
    int index, K;
    int odd = 0, R = 2;

    FiltreRIFDemiBande(const Vecteur<Tc> &c)
    {
      odd = 0;
      this->coefs = c;
      index = 0;
      K = coefs.rows();
      fenêtre.setZero(K);
    }

    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      tsd_assert(K > 0);
      int n = x.rows();
      auto iptr = x.data();
      auto optr = y.data();
      if(iptr != optr)
      {
        y.resize((n + odd) / R);
        optr = y.data();
      }

      for(auto j = 0; j < n; j++)
      {
        auto cptr = coefs.data();

        T somme = 0;

        // La moitié des coefficients sont nuls.
        // .... x3 x2 x1 x0
        // -> y0 = x0 * h0 + x2 * h1 + ...
        // -> y1 = x2 * h0 + ?

        // Donc on n'utilise que le signal pair ?
        // xpair   -> H0
        // ximpair -> H1 (nul)

        // Attention, 3 coefficients centrals non nuls puis zéros réguliérement.
        // Disons filtre de type h0 0 h1 0 h2
        //

        // Remplace les deux plus ancien éléments
        fenêtre(index) = *iptr++;
        index = (index + 1) % K;

        // Ne calcule qu'un échantillon sur R
        if(odd < R - 1)
        {
          odd++;
          continue;
        }
        odd = 0;


        T *wptr = fenêtre.data() + index;

        // Nombre d'échantillons à la fin de la ligne à retard
        int K1 = K - index;
        // Nombre d'échantillons au début de la ligne à retard
        int K2 = K - K1;

        int i;
        for(i = 0; i < K1; i += 2) // 1 coefficient sur 2 du filtre est nul
        {
          somme += *wptr * *cptr;
          wptr  += 2;
          cptr  += 2;
        }

        i = i - K1;

        // Redémarre au début de la ligne à retard
        wptr = fenêtre.data() + i;
        for(; i < K2; i += 2)
        {
          somme += *wptr * *cptr;
          wptr  += 2;
          cptr  += 2;
        }

        // Coefficient central non nul = 0.5 (forcément)
        somme += /*coefs(K/2)*/ 0.5f * fenêtre((index + K/2) % K);


        /*
         *
         * Pour i = 0; i < N; i += R
         * (fenêtre(index + i [N]) + fenêtre(index - 1 - i [N])) * coefs(i)
         *
         *
         *
         */

        *optr++ = somme;
      }
    }
  };



///////////////////////////////////////////////////////////////////////////////
// Filtre RIF avec décimation 1:R (suppose qu'un échantillon sur R est nul) ///
///////////////////////////////////////////////////////////////////////////////
template<typename T, typename Tc>
struct FiltreRIFDecim: FiltreGen<T>
{
  Vecteur<T> fenêtre;
  Vecteur<Tc> coefs;
  int index, K;
  int odd = 0, R = 0;

  FiltreRIFDecim(const Vecteur<Tc> &c, unsigned int R)
  {
    odd = 0;
    this->R = R;
    this->coefs = c;
    index = 0;
    K = coefs.rows();
    fenêtre.setZero(K);
  }

  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    tsd_assert(K > 0);
    int n = x.rows();
    auto iptr = x.data();
    auto optr = y.data();
    if(iptr != optr)
    {
      y.resize((n + odd) / R);
      optr = y.data();
    }

    for(auto j = 0; j < n; j++)
    {
      auto cptr = coefs.data();

      T somme = 0;

      // La moitié des coefficients sont nuls.
      // .... x3 x2 x1 x0
      // -> y0 = x0 * h0 + x2 * h1 + ...
      // -> y1 = x2 * h0 + ?

      // Donc on n'utilise que le signal pair ?
      // xpair   -> H0
      // ximpair -> H1 (nul)

      // Attention, 3 coefficients centrals non nuls puis zéros réguliérement.
      // Disons filtre de type h0 0 h1 0 h2
      //

      // Remplace les deux plus ancien éléments
      fenêtre(index) = *iptr++;
      index = (index + 1) % K;

      // Ne calcule qu'un échantillon sur R
      if(odd < R - 1)
      {
        odd++;
        continue;
      }
      odd = 0;


      T *wptr = fenêtre.data() + index;

      // Nombre d'échantillons à la fin de la ligne à retard
      int K1 = K - index;
      // Nombre d'échantillons au début de la ligne à retard
      int K2 = K - K1;

      for(auto i = 0; i < K1; i++)
        somme += *wptr++ * *cptr++;

      // Redémarre au début de la ligne à retard
      wptr = fenêtre.data();
      for(auto i = 0; i < K2; i++)
        somme += *wptr++ * *cptr++;



      /*
       *
       * Pour i = 0; i < N; i += R
       * (fenêtre(index + i [N]) + fenêtre(index - 1 - i [N])) * coefs(i)
       *
       *
       *
       */

      *optr++ = somme;
    }
  }
};


// Principe : filtrage polyphase
// au lieu d'insérer R-1 zéros et filtrer
// on partage mle filtre en R branches qui produisent chacun une partie de la sortie.
//
template<typename T, typename Tc>
struct FiltreRIFUps: FiltreGen<T>
{
  Vecteur<T> fenêtre;
  int index, K;
  Vecteur<Tc> coefs;
  int odd = 0, R = 0;


  FiltreRIFUps(const Vecteur<Tc> &c, int R)
  {
    this->R     = R;
    // Afin de préserver l'amplitude du signal
    this->coefs = c * R;
    odd   = 0;
    index = 0;
    K     = coefs.rows();

    // Complète avec des zéros
    if((K % R) != 0)
    {
      coefs = vconcat(coefs, Vecteur<Tc>::Zero(R - (K % R)));
      K = coefs.rows();
    }

    tsd_assert((K % R) == 0);

    fenêtre.setZero(K/R);
  }

  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    tsd_assert(K > 0);
    int n = x.rows();
    auto iptr = x.data();
    auto optr = y.data();
    if(iptr != optr)
    {
      y.resize(n*R);
      optr = y.data();
    }

    for(auto j = 0; j < n; j++)
    {

      // La moitié des coefficients sont nuls.
      // .... x3 x2 x1 x0
      // -> y0 = x0 * h0 + x2 * h1 + ...
      // -> y1 = x2 * h0 + ?

      // Donc on n'utilise que le signal pair ?
      // xpair   -> H0
      // ximpair -> H1 (nul)

      // Attention, 3 coefficients centrals non nuls puis zéros réguliérement.
      // Disons filtre de type h0 0 h1 0 h2
      //

      // Remplace les deux plus ancien éléments
      fenêtre(index) = *iptr++;
      index = (index + 1) % (K/R);

      // Pour chaque échantillon d'entrée, calcule R sorties différentes

      for(auto i = 0; i < R; i++)
      {
        T sum = 0;
        auto cptr = coefs.data() + (R - 1) - i;
        T *wptr = fenêtre.data() + index;

        int K1, K2;

        // Number of samples at end of delay line
        K1 = K/R - index;
        // Number at begin of delay line
        K2 = K/R - K1;

        for(auto i = 0; i < K1; i++)
        {
          sum += *wptr++ * *cptr;
          cptr += R;
        }

        wptr = fenêtre.data(); // Restart at beginning of delay line
        for(auto i = 0; i < K2; i++)
        {
          sum += *wptr++ * *cptr;
          cptr += R;
        }

        *optr++ = sum;
      }
    }
  }
};


template<typename Tc, typename T = Tc>
sptr<FiltreGen<T>> filtre_rif_decim(const Eigen::Ref<const Vecteur<Tc>> c, int R)
{
  return make_shared<FiltreRIFDecim<T,Tc>>(c, R);
}



template<typename Tc, typename T = Tc>
sptr<FiltreGen<T>> filtre_rif_demi_bande(const Eigen::Ref<const Vecteur<Tc>> c)
{
  return make_shared<FiltreRIFDemiBande<T,Tc>>(c);
}

template<typename Tc, typename T = Tc>
sptr<FiltreGen<T>> filtre_rif_ups(const Eigen::Ref<const Vecteur<Tc>> c, int R)
{
  return make_shared<FiltreRIFUps<T,Tc>>(c, R);
}


float filtre_rif_ups_délais(int nc, int R)
{
  int pad = 0;
  if((nc % R) != 0)
  {
    pad = R - (nc % R);
  }
  return (nc - 1) / 2.0 + pad;
}

// R = 2, nc = 15 -> pad = 1

float rif_delais(int nc)
{
  return (nc - 1) / 2.0f;
}

namespace hidden {
auto filtre_rif_decim1 = filtre_rif_decim<float, float>;
auto filtre_rif_decim2 = filtre_rif_decim<float, cfloat>;
auto filtre_rif_demi_bande1 = filtre_rif_demi_bande<float, float>;
auto filtre_rif_demi_bande2 = filtre_rif_demi_bande<float, cfloat>;
auto filtre_rif_ups1 = filtre_rif_ups<float, float>;
auto filtre_rif_ups2 = filtre_rif_ups<float, cfloat>;
}


}
















