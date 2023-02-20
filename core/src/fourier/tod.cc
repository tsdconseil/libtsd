#include "tsd/fourier/tod.hpp"
#include "tsd/tsd.hpp"
#include "tsd/filtrage/frat.hpp"

namespace tsd::tf::tod {

  std::ostream& operator<<(std::ostream& os, const Laurent &p)
  {
    os << format("({}) * z^-{}", p.polynome, p.n0);
    retourne os;
  }

  std::ostream& operator<<(std::ostream& os, const FormePolyphase &p)
  {
    os << "Matrice polyphase:\n"
       << "H00 = " << p.H00 << "\n"
       << "H01 = " << p.H01 << "\n"
       << "H10 = " << p.H10 << "\n"
       << "H11 = " << p.H11 << "\n";
    retourne os;
  }

  std::ostream& operator<<(std::ostream& os, const QMF &p)
  {
    os << "Forme QMF:\n"
       << "H0 = " << p.H0 << "\n"
       << "H1 = " << p.H1 << "\n"
       << "G0 = " << p.G0 << "\n"
       << "G1 = " << p.G1 << "\n";
    retourne os;

  }

  // a + b * c
  Laurent laurent_mac(const Laurent &a, const Laurent &b, const Laurent &c)
  {
    Laurent bc;
    bc.polynome = b.polynome * c.polynome;
    bc.n0       = b.n0 + c.n0;

    Laurent r;

    r.n0 = min(a.n0, bc.n0);

    Poly<float> ashift  = a.polynome << (a.n0 - r.n0);
    Poly<float> bcshift = bc.polynome << (bc.n0 - r.n0);

    r.polynome = ashift + bcshift;

    retourne r;
  }

  FormePolyphase::FormePolyphase(const Lift &lift)
  {
    H00.polynome = Poly<float>::one();
    H11.polynome = Poly<float>::one();
    H01.polynome = Poly<float>(0.0f);
    H10.polynome = Poly<float>(0.0f);

    pour(const LiftElem &etape: lift.etapes)
    {
      si(etape.predict)
      {
        H10 = laurent_mac(H10, etape.polynome, H00);
        H11 = laurent_mac(H11, etape.polynome, H01);
      }
      sinon
      {
        H00 = laurent_mac(H00, etape.polynome, H10);
        H01 = laurent_mac(H01, etape.polynome, H11);
      }
    }

    H00.polynome = H00.polynome * lift.K;
    H01.polynome = H01.polynome * lift.K;
    H10.polynome = H10.polynome * (1.0f / lift.K);
    H11.polynome = H11.polynome * (1.0f / lift.K);
  }

  QMF::QMF(const FormePolyphase &fp)
  {
    // (1) rend tous les filtres causaux
    Poly<float> H00, H11, H01, H10;

    entier md = -min({fp.H00.n0, fp.H01.n0, fp.H10.n0, fp.H11.n0});

    tsd_assert(md >= 0);

    H00 = fp.H00.polynome << md;
    H01 = fp.H01.polynome << md;
    H10 = fp.H10.polynome << md;
    H11 = fp.H11.polynome << md;

    soit z  = Poly<float>::z;
    soit z2 = z * z;
    H0 = H00.horner(z2) + z * H01.horner(z2);
    H1 = H10.horner(z2) + z * H11.horner(z2);

    // TODO : G

  }


  // n=8: Split(0,1,2,3|4,5,6,7)   = (0,2,4,6|1,3,5,7)
  // half = 4
  // n=9: Split(0,1,2,3|4,5,6,7,8) = (0,2,4,6,8|1,3,5,7)
  // half = 5
  /** @brief Sépare les éléments pairs et impairs d'un signal */
  template<typename T>
  void split(Vecteur<T> &x, entier n)
  {
    //entier n = x.rows();
    entier half = (n + 1) >> 1;

    // Stockage des coefs pairs


    soit tmp = Vecteur<T>::int_expr(half,
        IMAP(x(2*i)));

    // Ecriture des coefs impairs
    // x[n-1] := x[n-1] (to remove)
    // x[n-2] := x[n-3]
    // ....
    // x[n-half] := x[n-(n+1)+1] = x[0] si n impair
    // x[n-half] := x[n-n+1] = x[1] si n pair

    entier max = ((n & 1) == 0) ? half : half - 1;
    entier dec = ((n & 1) == 0) ? 0 : 1;
    pour(auto i = 1; i <= max; i++)
      x(n-i) = x(n-(2*i)+1-dec);

    x.head(half) = tmp;
  }


  // n=8: Split(0,1,2,3|4,5,6,7)   = (0,2,4,6|1,3,5,7)
  // half = 4
  // n=9: Split(0,1,2,3|4,5,6,7,8) = (0,2,4,6,8|1,3,5,7)
  // half = 5
  template<typename T>
  void merge(Vecteur<T> &x, entier n)
  {
    //entier n = x.rows();
    entier half = (n + 1) >> 1;

    // Stockage coefs pairs
    Vecteur<T> tmp = x.head(half);

    // Ecriture coefs impairs
    entier max = ((n & 1) == 0) ? half : half - 1;
    pour(auto i = 0; i < max; i++)
      x(2*i+1) = x(i+half);

    pour(auto i = 0; i < half; i++)
      x(2*i) = tmp(i);
  }

/*static Laurent lift_element(bouléen pred, entier n0, const tsd::Poly<float> &p)
{
  retourne Laurent{p, n0, pred};
}*/

static const auto z0 = tsd::Poly<float>::one();
static const auto z  = tsd::Poly<float>::z;

Lift lift_haar()
{
  Lift res;

  res.nom = "haar";
  res.K = sqrt(2);


  Laurent p{-z0, 0};
  Laurent m{0.5f * z0, 0};

  res.etapes = {LiftElem{p, oui}, LiftElem{m, non}};

  retourne res;
}

Lift lift_db2()
{
  Lift res;

  float s3 = sqrt(3);

  res.nom = "db2";
  res.K = (s3-1)/sqrt(2);

  LiftElem p1{{sqrt(3.0f) * z0, 0}, non};
  LiftElem p2{{-(s3-2)/4 - (s3/4) * z, -1}, oui};
  LiftElem p3{{-z, 0}, non};

  /*p2.predict  = oui;
  p2.n0 = -1;
  ArrayXf c(2);
  c << -(s3-2)/4, -s3/4;
  c(0) = -(s3-2)/4;
  c(1) = -s3/4;
  p2.polynome = {c};

  msg("p2 : {}", p2.polynome);*/

  //p3.predict  = non;
  //p3.polynome = -tsd::Poly<float>::z;

  res.etapes = {p1, p2, p3};
  retourne res;
}


template<typename T = float, typename Tcoefs = float>
struct OndeletteGen: Ondelette<T>
{
  Lift lift;
  OndeletteGen(const Lift &lift)
  {
    this->nom  = lift.nom;
    this->lift = lift;
  }
  void lift_step(Vecteur<T> &x, entier n)
  {
    entier half = n >> 1;

    split(x, n);

    pour(auto &etape: lift.etapes)
    {
      entier idxi = etape.predict ? 0 : half;
      entier idxo = etape.predict ? half : 0;

      pour(auto j = 0; j < half; j++)
      {
        float sm = 0;
        pour(auto l = 0; l < etape.polynome.polynome.coefs.rows(); l++)
        {
          si((j + etape.polynome.n0 + l >= 0) && (j + etape.polynome.n0 + l < half))
            sm += etape.polynome.polynome.coefs(l) * x(idxi+j+etape.polynome.n0+l);
        }
        x(idxo+j) += sm;
      }
    }
  }
  void ilift_step(Vecteur<T> &x, entier n)
  {
    entier half = n >> 1;

    pour(auto i = lift.etapes.rbegin(); i != lift.etapes.rend(); i++)
    {
      soit &etape = *i;
      entier idxi = etape.predict ? 0 : half;
      entier idxo = etape.predict ? half : 0;

      pour(auto j = 0; j < half; j++)
      {
        float sm = 0;
        pour(auto l = 0; l < etape.polynome.polynome.coefs.rows(); l++)
        {
          si((j + etape.polynome.n0 + l >= 0) && (j + etape.polynome.n0 + l < half))
            sm += etape.polynome.polynome.coefs(l) * x(idxi+j+etape.polynome.n0+l);
        }
        x(idxo+j) -= sm;
      }
    }

    merge(x, n);
  }
};




#if 0
template<typename T = float, typename Tcoefs = float>
struct OndeletteBiortho35: Ondelette<T>
{
  OndeletteBiortho35()
  {
    this->nom = "bior-3-5";
  }
  void lift_step(Vecteur<T> &x, entier n)
  {
    //entier n = x.rows();
    entier half = n >> 1;

    split(x, n);

    // predict
    pour(auto j = 0; j < half - 1; j++)
      x(half+j) -= (x(j) + x(j+1)) / 2;
    x(n-1) -= x(half-1);

    // update
    x(0)   += (x(half) + 1) / 2;
    pour(auto j = 1; j < half; j++)
      x(j)   += (x(half+j) + x(half+j-1) + 2) / 4;
  }
  void ilift_step(Vecteur<T> &x, entier n)
  {
    //entier n = x.rows();
    entier half = n >> 1;

    // update
    x(0)   -= (x(half) + 1) / 2;
    pour(auto j = 1; j < half; j++)
      x(j)   -= (x(half+j) + x(half+j-1) + 2) / 4;

    // predict
    pour(auto j = 0; j < half - 1; j++)
      x(half+j) += (x(j) + x(j+1)) / 2;
    x(n-1) += x(half-1);

    merge(x, n);
  }
};


template<typename T = float, typename Tcoefs = float>
struct OndeletteHaar: Ondelette<T>
{
  OndeletteHaar()
  {
    this->nom = "haar";
  }
  void lift_step(Vecteur<T> &x, entier n)
  {
    //entier n = x.rows();
    entier half = (n + 1) >> 1;

    Vecteur<T> y(n);

    split(x, n);

    pour(auto i = 0; i < n/2; i++)
      x(half+i) -= x(i);
    pour(auto i = 0; i < n/2; i++)
      x(i) += 0.5 * x(half+i);

    /*pour(auto i = 0; i < n/2; i++)
      x(2*i+1) = x(2*i+1) - x(2*i);
    pour(auto i = 0; i < n/2; i++)
      x(2*i) = x(2*i) - 0.5 * x(2*i+1);*/

    // PB : pas splitté à la fin
  }
  void ilift_step(Vecteur<T> &x, entier n)
  {
    //entier n = x.rows();
    entier half = (n + 1) >> 1;

    pour(auto i = 0; i < n/2; i++)
      x(i) -= 0.5 * x(half+i);
    pour(auto i = 0; i < n/2; i++)
      x(half+i) += x(i);

    merge(x, n);
  }
};

template<typename T = float, typename Tcoefs = float>
struct OndeletteDb4: Ondelette<T>
{
  Tcoefs sqr3, cf1, nrm1, nrm2;
  OndeletteDb4()
  {
    sqr3 = (Tcoefs) std::sqrt(3);
    cf1  = (Tcoefs) ((2.0 - std::sqrt(3.0)) / 4.0);
    nrm1 = (Tcoefs) ((std::sqrt(3)+1) / std::sqrt(2));
    nrm2 = (Tcoefs) ((std::sqrt(3)-1) / std::sqrt(2));
    this->nom  = "db4";
  }

  // http://fr.wikipedia.org/wiki/Lifting_en_ondelettes
  // http://ronan.lepage1.free.fr/repository/ondelettes/lifting13/node13.html
  void lift_step(Vecteur<T> &x, entier n)
  {
    //entier n = x.rows();
    entier half = (n + 1) >> 1;

    split(x, n);

    si((n & 1) == 0)
    {
      //pour(j = 0; j < half; j++)
      //  x(j) += sqr3 * x(j+half);
      x.head(half) += sqr3 * x.segment(half, half);

      //(sqrt3/4.0)*S[0] + (((sqrt3-2)/4.0)*S[0]);
      x(half) -= (sqr3 * x(0)) / 4 - cf1 * x(0);
      pour(auto j = 1; j < half; j++)
      {
        // sqrt3/4.0)*S[n] + (((sqrt3-2)/4.0)*S[n-1]);
        x(half+j) -= (sqr3 * x(j)) / 4 - cf1 * x(j-1);
      }

      //x.head(half-1) -= x.segment(half+j, half-1);

      pour(auto j = 0; j + 1 < half; j++)
        x(j) -= x(half+j+1);

      x(half-1) -= x(n-1);
    }
    sinon
    {
      // n=8: Split(0,1,2,3|4,5,6,7)   = (0,2,4,6|1,3,5,7)
      // half = 4
      // n=9: Split(0,1,2,3|4,5,6,7,8) = (0,2,4,6,8|1,3,5,7)
      // half = 5
      //pour(j = 0; j + 1 < half; j++)
      //  x[j] += sqr3 * x[j+half];

      x.head(half-1) += sqr3 * x.segment(half, half-1);

      x(half-1) += sqr3 * x(n-1);

      //(sqrt3/4.0)*S[0] + (((sqrt3-2)/4.0)*S[0]);
      x(half) -= (sqr3 * x(0)) / 4 - cf1 * x(0);
      pour(auto j = 1; j + 1 < half; j++)
      {
        // sqrt3/4.0)*S[n] + (((sqrt3-2)/4.0)*S[n-1]);
        x(half+j) -= (sqr3 * x(j)) / 4 - cf1 * x(j-1);
      }

      pour(auto j = 0; j + 2 < half; j++)
        x(j) -= x(half+j+1); // x0 -= x3, x2 -= x5, x4 -= x7
                             // x

      x(half-2) -= x(n-1); // x6 -= x7
      x(half-1) -= x(n-1); // x8 -= x7
    }

    // Normalisation

    x.head(half)     *= nrm2;
    x.tail(n - half) *= nrm1;

    //pour (j = 0; j < half; j++)
    //  x[j] *= nrm2;
    //pour (j = half; j < n; j++)
    //  x[j] *= nrm1;
  }

  void ilift_step(Vecteur<T> &x, entier n)
  {
    //entier n = x.rows();
    entier half = (n + 1) >> 1;

    // Dénormalisation
    //pour (j = 0; j < half; j++)
    //  x[j] *= nrm1;
    //pour (j = half; j < n; j++)
    //  x[j] *= nrm2;
    x.head(half)     *= nrm1; // /= ?
    x.tail(n - half) *= nrm2;


    si((n & 1) == 0)
    {
      pour(auto j = 0; j + 1 < half; j++)
        x(j) += x(half + j + 1);
      x(half-1) += x(n-1);

      x(half) += (sqr3 * x(0)) / 4 - cf1 * x(0);
      pour (auto j = 1; j < half; j++)
        x(half+j) += (sqr3 * x(j)) / 4 - cf1 * x(j-1);

      pour(auto j = 0; j < half; j++)
        x(j) -= sqr3 * x(j+half);
    }
    sinon
    {
      // n=8: Split(0,1,2,3|4,5,6,7)   = (0,2,4,6|1,3,5,7)
      // half = 4
      // n=9: Split(0,1,2,3|4,5,6,7,8) = (0,2,4,6,8|1,3,5,7)
      // half = 5
      pour(auto j = 0; j + 2 < half; j++)
        x(j) += x(half + j + 1); // x0 += x3, x2 += 5, x4 += x7
      x(half-2) += x(n-1); // x6 += x7
      x(half-1) += x(n-1); // x8 += x7

      x(half) += (sqr3 * x(0)) / 4 - cf1 * x(0); // x1 += ...(x0,x0)
      pour(auto j = 1; j + 1 < half; j++)
        x(half+j) += (sqr3 * x(j)) / 4 - cf1 * x(j-1); // (x3,x2,x0), ..., (x7,x6,x4)

      pour(auto j = 0; j + 1 < half; j++)
        x(j) -= sqr3 * x(j+half);
      x(half-1) -= sqr3 * x(n-1);
    }

    merge(x, n);
  }


  };
#endif


/** @brief Transformation en ondelette discréte entière
 *  @param vec Tableau contenant le signal en fonction du temps
 *  @param N   Taille du tableau
 *  @require N = 2^k
 *  La transformation effectuée est d'ordre log2(N), c'est-à-dire
 *  qu'à la fin il n'y qu'un coefficient d'échelle (première position du
 *  tableau). Tous les résultats sont stockés "sur place", soit dans
 *  le tableau passé en paramètre. */
template<typename T>
void dwt(sptr<Ondelette<T>> wavelet, Vecteur<T> &x, entier depth)
{
  entier N = x.rows();
  entier last = N >> (depth - 1);

  pour(auto n = N; n >= last; n = n >> 1)
    wavelet->lift_step(x, n);
}


/** @brief Transformation en ondelette discréte entière inverse
 *  @param vec Tableau contenant les coefficient d'échelle et d'ondelette
 *  @param N   Taille du tableau
 *  @require N = 2^k
 *  Idem transformation directe, mais restitue les valeurs originales. */
template<typename T>
void iwt(sptr<Ondelette<T>> wavelet, Vecteur<T> &x, entier depth)
{
  entier N = x.rows();
  entier start = N >> (depth - 1);

  pour(auto n = start; n <= N; n = n << 1)
    wavelet->ilift_step(x, n);
}


#if 0
template<typename T>
sptr<Ondelette<T>> ondelette_db4()
{
  retourne std::make_shared<OndeletteDb4<T>>();
}

template<typename T>
sptr<Ondelette<T>> ondelette_biortho_3_5()
{
  retourne std::make_shared<OndeletteBiortho35<T>>();
}

template<typename T>
sptr<Ondelette<T>> ondelette_haar()
{
  retourne std::make_shared<OndeletteHaar<T>>();
}
#endif

template<typename T>
sptr<Ondelette<T>> ondelette_gen(const Lift &lift)
{
  retourne std::make_shared<OndeletteGen<T>>(lift);
}

template
void dwt<float>(sptr<Ondelette<float>> wavelet, Vecteur<float> &x, entier depth);
template
void iwt<float>(sptr<Ondelette<float>> wavelet, Vecteur<float> &x, entier depth);

//template
//sptr<Ondelette<float>> ondelette_db4();

//template
//sptr<Ondelette<float>> ondelette_biortho_3_5();

//template
//sptr<Ondelette<float>> ondelette_haar();

template
sptr<Ondelette<float>> ondelette_gen(const Lift &lift);

#if 0
/** @brief Quantification des coefficients d'ondelettes
 *  Réduit la résolution des premiers coefficents (non utilisé pour SPIHT)
 *  @note Attention: après quantification, le carré des  coefficients n'est
 *        plus proportionnel  à l'énergie qu'ils représentent dans le signal.
 *  c_uk = floor(cuk/2^(u/2) + 0.5)
 */
template<typename T>
void quantify(dsp::Matrix<T,1> &x)
{
  unsigned short i, n = x.extent(0);
  unsigned short u_2 = 1;
  unsigned short half = n >> 1;
  pour(;;)
  {
    pour(i = half; i < 2 * half; i++)
      x(i) = (T) floor(((float) x(i))/((float)u_2) + 0.5);
    u_2 = u_2 << 1;
    half = half >> 1;
    si(half == 0)
    {
      x(0) = (T) floor(((float) x(i))/((float)u_2) + 0.5);
      break;
    }
  }
}

/** @brief Opération inverse */
template<typename T>
void dequantify(dsp::Matrix<T,1> &x)
{
  unsigned short i, n = x.extent(0);
  unsigned short u_2 = 1;
  unsigned short half = n >> 1;
  pour(;;)
  {
    pour(i = half; i < 2 * half; i++)
      x(i) = u_2 * x(i);
    u_2 = u_2 << 1;
    half = half >> 1;
    si(half == 0)
    {
      x(i) = u_2 * x(i);
      break;
    }
  }
}
#endif


}

