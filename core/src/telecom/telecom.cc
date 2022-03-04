#include "tsd/telecom.hpp"
#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/vue.hpp"
#include <iostream>

using namespace tsd::vue;

namespace tsd::telecom {


  BitStream code_Barker(int n)
  {
    switch(n)
    {
    case 2:
      return BitStream("10");
    case 3:
      return BitStream("110");
    case 4:
      return BitStream("1101");
    case 5:
      return BitStream("11101");
    case 7:
      return BitStream("1110010");
    case 11:
      return BitStream("11100010010");
    case 13:
      return BitStream("1111100110101");
    default:
      echec("Code Barker (n = {}) : non supporté.", n);
      return BitStream(); // non atteint
    }
  }


/*ArrayXf pam_modulate(const BitStream &x)
{
  int n = x.lon();
  ArrayXf y(n);
  for(auto i = 0; i < n; i++)
    y(i) = x[i] ? 1.0f : -1.0f;
  return y;
}*/



  CmpBitsRes cmp_bits(const BitStream &bs1, const BitStream &bs2)
  {
    CmpBitsRes res;

    ArrayXf b1 = bs1.array(), b2 = bs2.array();

    ArrayXcf b1c = b1, b2c = b2;
    auto [b1a, b2a, dt, score] = tsd::fourier::aligne_entier(b1c, b2c);

    ArrayXf b1ar = b1a.real(), b2ar = b2a.real();

    ArrayXf b1ap = b1ar, b2ap = b2ar;


    // TODO : utiliser plutot segment dans le comptage d'erreurs
    if(b1a.rows() > 100)
    {
      auto n = b1a.rows();
      b1ap = b1ar.tail(n - 50).head(n-100).eval();
      b2ap = b2ar.tail(n - 50).head(n-100).eval();
    }

    res.nerr = 0;
    for(auto i = 0; i < b1ap.rows(); i++)
    {
      if(std::round(b1ap(i)) != std::round(b2ap(i)))
        res.nerr++;
    }

    //res.nerr = (b1ap-b2ap).abs().sum();

    if(b2a.rows() == 0)
      res.ber = NAN;
    else
      res.ber = ((float) res.nerr) / b2ap.rows();

    msg("Cmp bits : dec = {}, nerr = {}, rows = {}, ber = {:.1e}", dt, res.nerr, b2ap.rows(), res.ber);

    res.decalage  = dt;
    res.score = score;
    res.b0 = b1ar;
    res.b1 = b2ar;

    if(0)
    {
      Figures f;
      f.subplot().plot(b1, "hb", "Sequence d'entree ({} bits)", b1.rows());
      f.subplot().plot(res.b0, "hb", "Sequence d'entree recalée");
      f.subplot().plot(b2, "hb", "Sequence de sortie ({} bits)", b2.rows());
      f.subplot().plot(res.b1, "hb", "Sequence de sortie recalée");
      f.subplot().plot(res.b1 - res.b0, "hr", "Erreurs ({}, ber = {:.1e})", res.nerr, res.ber);
      f.afficher(fmt::format("Comparaison bits (dec={}, s={})", dt,score));
    }


    return res;
  }


  CmpBitsRes cmp_bits_psk(const BitStream &bs0, const BitStream &bs1, int k)
  {
    int M = (1 << k);
    CmpBitsRes r[M];
    CmpBitsRes &br = r[0];
    br.b1 = bs1.array();

    // Essai éventuellement plusieurs phases possibles
    for(auto l = 0; l < M; l++)
    {
      ArrayXi idx = symmap_binaire(bs1, k);

      for(auto i = 0; i < idx.rows(); i++)
        idx(i) = ((int) idx(i) + l) % M;

      BitStream bs;
      symdemap_binaire(bs, idx, k);
      //ArrayXf b1p = bs.vers_array();
      //ArrayXf b1p = wf->decode_symboles(symbs);

      msg("cmp bits...");
      r[l] = cmp_bits(bs0, bs);
      msg("  ...cmp bits ok.");

      msg(" test phase = {} : score = {} (nerr = {})", l, r[l].score, r[l].nerr);

      if((l == 0) || (r[l].nerr < br.nerr))
      {
        br           = r[l];
        br.dec_phase = l;
      }
      //infos("Essai dec phase = %d : corr = %.2f, nerr = %d", k, r[k].score, r[k].nerr);
    }
    //br.b0 = b0;
    return br;
  }




/*void OscillateurHarmonique::affiche_etat() const
{
  infos("OLH : freq = %f, omega = %f degrés.", freq, 2 * 180 * freq);
  infos("Rotation = %f, %f", rotation.real(), rotation.imag());
}*/



/*void OscillateurHarmonique::configure(cfloat rotation)
{
  this->rotation = rotation;
}

void OscillateurHarmonique::step(ArrayXXcf &io)
{
  step(io, io);
}*/

/*ArrayXf OscillateurHarmonique::avance_cos(unsigned int n)
{
  auto rot = avance(n);
  return rot.real();
}*/



/*ArrayXcf OscillateurHarmonique::avance(unsigned int n)
{
  return step(n);
}

ArrayXcf OLH_signal(unsigned int n, float nu)
{
  OscillateurHarmonique oh;
  oh.configure(nu);
  return oh.avance(n);
}*/

#if 0
void OscillateurHarmonique::step(const ArrayXcf &in, ArrayXcf &out)
{
  auto rot = avance(in.rows());

  out = (in * rot).eval();

# if 0
  auto nrows = in.rows();
  out.resize(nrows);

  phaseur = phaseur / std::abs(phaseur); // Periodic update

  auto *iptr = in.data();
  auto *optr = out.data();

  for(auto i = 0; i < nrows; i++)
  {
    *optr = *iptr * phaseur;
    iptr += nrows;
    optr += nrows;
    phaseur *= rotation;
    iptr -= (nrows /** ncols*/ - 1);
    optr -= (nrows /** ncols*/ - 1);
  }
# endif
}

    
void OscillateurHarmonique::step(const ArrayXXcf &in, ArrayXXcf &out)
{
  auto ncols = in.cols(), nrows = in.rows();
  auto rot = avance(nrows);

  if(out.data() != in.data())
      out.resize(nrows, ncols);

  for(auto j = 0; j < ncols; j++)
    out.col(j) = in.col(j) * rot;


# if 0
  auto ncols = in.cols(), nrows = in.rows();

  if(out.data() != in.data())
    out.resize(nrows, ncols);

  phaseur = phaseur / std::abs(phaseur); // Periodic update

  auto *iptr = in.data();
  auto *optr = out.data();

  for(auto i = 0; i < nrows; i++)
  {
    for(auto j = 0; j < ncols; j++)
    {
      *optr = *iptr * phaseur;
      iptr += nrows;
      optr += nrows;
    }
    phaseur *= rotation;
    iptr -= (nrows * ncols - 1);
    optr -= (nrows * ncols - 1);
  }
# endif
}

#if 0
class OscillateurHarmonique
{
public:
  /** @param nu: Fréquence normalisée, entre -0.5 et 0.5 */
  void configure(float nu);
  void configure(cfloat rotation);

  void step(const ArrayXcf &in, ArrayXcf &out);
  void step(const ArrayXXcf &in, ArrayXXcf &out);
  void step(ArrayXXcf &io);

  ArrayXf avance_cos(unsigned int n);
  ArrayXcf avance(unsigned int n);
  ArrayXcf step(unsigned int n);

  void affiche_etat() const;

//private:
  float freq = 1;
  cfloat phaseur = 1, rotation = 1;
};
using OLH = OscillateurHarmonique;
#endif
#endif




void diff_encode(BitStream &y, const BitStream &x)
{
  int n = x.lon(), lb = 0;

  for(auto i = 0; i < n; i++)
  {
    int bi = x[i];
    int bo = lb ^ bi;
    y.push(bo);
    lb = bo;
  }
}

void diff_decode(BitStream &y, const BitStream &x)
{
  int n = x.lon();

  if(n == 0)
    return;

  int lb = x[0];

  for(auto i = 0; i < n-1; i++)
  {
    int bi = x[i];
    int bo = lb ^ bi;
    y.push(bo);
    lb = bi;
  }
}

void decode_hard(BitStream &y, const ArrayXf &llr)
{
  int n = llr.rows();
  for(auto i = 0; i < n; i++)
    y.push(llr(i) > 0 ? 1 : 0);
}

ArrayXcf bruit_awgn(IArrayXcf &x, float σ)
{
  auto n = x.rows();
  ArrayXcf y(n);
  y.real() = x.real() + σ * randn(n);
  y.imag() = x.imag() + σ * randn(n);
  return y;
}

ArrayXf bruit_awgn(IArrayXf &x, float σ)
{
  return x + σ * randn(x.rows());
}



}



