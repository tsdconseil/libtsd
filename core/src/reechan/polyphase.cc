///////////////////////////////////////
// Conception de filtres demi-bande
///////////////////////////////////////

#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"

using namespace std;
using namespace tsd::vue;


namespace tsd::filtrage {


  template<typename T>
  TabT<T,2> forme_polyphase(const Vecteur<T> &x, unsigned int M)
  {
    soit n = x.rows();

    si(n == 0)
      retourne {};

    soit r = n % M;
    Vecteur<T> x2;
    si(r != 0)
      x2 = x | Vecteur<T>::zeros(M-r);
    sinon
      x2 = x.clone();
    n = x2.rows();

    // n / M = nb échantillons par canal
    soit X = x2.reshape(M, n / M);

    retourne X;
    // Inversion des lignes
    //retourne X.reverse_rows();
  }


  template<typename T>
  Vecteur<T> iforme_polyphase(const TabT<T,2> &x)
  {
    retourne x.reshape(x.cols() * x.rows());
  }






  // Implémentation optimisée pour filtre demi-bande
  template<typename T, typename Tc>
  struct FiltreRIFDemiBande: FiltreGen<T>
  {
    Vecteur<T> fenêtre;
    Vecteur<Tc> coefs;
    entier index = 0, K;
    entier odd = 0, R = 2;

    FiltreRIFDemiBande(const Vecteur<Tc> &c)
    {
      coefs = c;
      K     = coefs.rows();
      fenêtre = Vecteur<T>::zeros(K);
    }

    void step(const Vecteur<T> &x, Vecteur<T> &y)
    {
      tsd_assert(K > 0);
      soit n    = x.rows();
      soit iptr = x.data();
      soit optr = y.data();
      si(iptr != optr)
      {
        y.resize((n + odd) / R);
        optr = y.data();
      }

      pour(auto j = 0; j < n; j++)
      {
        soit cptr = coefs.data();

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
        si(odd < R - 1)
        {
          odd++;
          continue;
        }
        odd = 0;

        soit wptr = fenêtre.data() + index;

        // Nombre d'échantillons à la fin de la ligne à retard
        soit K1 = K - index;
        // Nombre d'échantillons au début de la ligne à retard
        soit K2 = K - K1;

        entier i;
        pour(i = 0; i < K1; i += 2) // 1 coefficient sur 2 du filtre est nul
        {
          somme += *wptr * *cptr;
          wptr  += 2;
          cptr  += 2;
        }

        i = i - K1;

        // Redémarre au début de la ligne à retard
        wptr = fenêtre.data() + i;
        pour(; i < K2; i += 2)
        {
          somme += *wptr * *cptr;
          wptr  += 2;
          cptr  += 2;
        }

        // Coefficient central non nul = 0.5 (forcément)
        somme += /*coefs(K/2)*/ 0.5f * fenêtre((index + K/2) % K);


        /*  pour i = 0; i < N; i += R
         * (fenêtre(index + i [N]) + fenêtre(index - 1 - i [N])) * coefs(i)
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
  entier index = 0, K;
  entier cnt = 0, R;

  FiltreRIFDecim(const Vecteur<Tc> &c, entier R)
  {
    this->R = R;
    coefs = c;
    K = coefs.rows();
    fenêtre = Vecteur<T>::zeros(K);
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    tsd_assert(K > 0);
    soit n = x.rows();
    soit iptr = x.data();
    soit optr = y.data();
    si(iptr != optr)
    {
      y.resize((n + cnt) / R);
      optr = y.data();
    }

    pour(auto j = 0; j < n; j++)
    {
      soit cptr = coefs.data();

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
      si(cnt < R - 1)
      {
        cnt++;
        continue;
      }
      cnt = 0;


      soit wptr = fenêtre.data() + index;

      // Nombre d'échantillons à la fin de la ligne à retard
      soit K1 = K - index;
      // Nombre d'échantillons au début de la ligne à retard
      soit K2 = K - K1;

      pour(auto i = 0; i < K1; i++)
        somme += *wptr++ * *cptr++;

      // Redémarre au début de la ligne à retard
      wptr = fenêtre.data();
      pour(auto i = 0; i < K2; i++)
        somme += *wptr++ * *cptr++;



      /* pour i = 0; i < N; i += R
       * (fenêtre(index + i [N]) + fenêtre(index - 1 - i [N])) * coefs(i) */

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
  entier index, K;
  Vecteur<Tc> coefs;
  entier odd = 0, R = 0;


  FiltreRIFUps(const Vecteur<Tc> &c, entier R)
  {
    this->R     = R;
    // Afin de préserver l'amplitude du signal
    this->coefs = c * R;
    odd   = 0;
    index = 0;
    K     = coefs.rows();

    // Complète avec des zéros
    si((K % R) != 0)
    {
      coefs = vconcat(coefs, Vecteur<Tc>::zeros(R - (K % R)));
      K = coefs.rows();
    }

    tsd_assert((K % R) == 0);

    fenêtre = Vecteur<T>::zeros(K/R);

  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    tsd_assert(K > 0);
    soit n    = x.rows();
    soit iptr = x.data();
    soit optr = y.data();
    si(iptr != optr)
    {
      y.resize(n*R);
      optr = y.data();
    }

    pour(auto j = 0; j < n; j++)
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

      // pour chaque échantillon d'entrée, calcule R sorties différentes

      pour(auto i = 0; i < R; i++)
      {
        T sum = 0;
        soit cptr = coefs.data() + (R - 1) - i;
        soit wptr = fenêtre.data() + index;

        // Number of samples at end of delay line
        soit K1 = K/R - index;
        // Number at begin of delay line
        soit K2 = K/R - K1;

        pour(auto i = 0; i < K1; i++)
        {
          sum += *wptr++ * *cptr;
          cptr += R;
        }

        wptr = fenêtre.data(); // Restart at beginning of delay line
        pour(auto i = 0; i < K2; i++)
        {
          sum += *wptr++ * *cptr;
          cptr += R;
        }

        *optr++ = sum;
      }
    }
  }
};


template<typename Tc, typename T>
sptr<FiltreGen<T>> filtre_rif_decim(const Vecteur<Tc> &c, entier R)
{
  retourne make_shared<FiltreRIFDecim<T,Tc>>(c, R);
}

template<typename Tc, typename T>
sptr<FiltreGen<T>> filtre_rif_demi_bande(const Vecteur<Tc> &c)
{
  retourne make_shared<FiltreRIFDemiBande<T,Tc>>(c);
}

template<typename Tc, typename T>
sptr<FiltreGen<T>> filtre_rif_ups(const Vecteur<Tc> &c, entier R)
{
  retourne make_shared<FiltreRIFUps<T,Tc>>(c, R);
}


float filtre_rif_ups_délais(entier nc, entier R)
{
  soit pad = 0;
  si((nc % R) != 0)
    pad = R - (nc % R);
  retourne (nc - 1) / 2.0 + pad;
}

// R = 2, nc = 15 -> pad = 1

float rif_delais(entier nc)
{
  retourne (nc - 1) / 2.0f;
}

namespace hidden {
soit filtre_rif_decim1 = filtre_rif_decim<float, float>;
soit filtre_rif_decim2 = filtre_rif_decim<float, cfloat>;
soit filtre_rif_demi_bande1 = filtre_rif_demi_bande<float, float>;
soit filtre_rif_demi_bande2 = filtre_rif_demi_bande<float, cfloat>;
soit filtre_rif_ups1 = filtre_rif_ups<float, float>;
soit filtre_rif_ups2 = filtre_rif_ups<float, cfloat>;
soit forme_polyphase1 = forme_polyphase<float>;
soit forme_polyphase2 = forme_polyphase<cfloat>;
soit iforme_polyphase1 = iforme_polyphase<float>;
soit iforme_polyphase2 = iforme_polyphase<cfloat>;

}

/* TODEL template
  TabT<float,2> forme_polyphase<float>(const Vecteur<float> &x, unsigned int M);

template
  TabT<cfloat,2> forme_polyphase<cfloat>(const Vecteur<cfloat> &x, unsigned int M);

template
  Vecteur<float> iforme_polyphase<float>(const TabT<float,2> &x);

template
  Vecteur<cfloat> iforme_polyphase<cfloat>(const TabT<cfloat, 2> &x);*/


}











# if 0
  template<typename T, typename Tc>
  struct FiltreRIFNv: FiltreGen<T>
  {
    Vecteur<T> fenêtre;
    Vecteur<Tc> coefs;
    entier index, K;
    entier cnt = 0, R = 0;

    Vecteur<T> futur;

    FiltreRIFNv(const Vecteur<Tc> &c, unsigned int R)
    {
      cnt = 0;
      this->R = R;
      this->coefs = c;
      index = 0;
      K = coefs.rows();
      fenêtre.setZero(K);
    }

    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      tsd_assert(K > 0);
      entier n = x.rows();
      auto iptr = x.data();
      auto optr = y.data();
      si(iptr != optr)
      {
        y.resize((n + cnt) / R);
        optr = y.data();
      }

      pour(auto j = 0; j < n; j++)
      {
        auto cptr = coefs.data();
        auto xi = *iptr++;

        pour(auto i = 0; i < n; i++)
          futur((index + i) % K) += xi * *cptr++;

        *optr++ = futur(0);
        index++;
      }
    }
  };
# endif





