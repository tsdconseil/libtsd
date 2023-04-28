#include "tsd/telecom.hpp"
#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"
#include "tsd/vue.hpp"
#include <iostream>

using namespace tsd::vue;

namespace tsd::telecom {


  BitStream code_Barker(entier n)
  {
    switch(n)
    {
    case 2:
      retourne BitStream("10");
    case 3:
      retourne BitStream("110");
    case 4:
      retourne BitStream("1101");
    case 5:
      retourne BitStream("11101");
    case 7:
      retourne BitStream("1110010");
    case 11:
      retourne BitStream("11100010010");
    case 13:
      retourne BitStream("1111100110101");
    default:
      echec("Code Barker (n = {}) : non supporté.", n);
      retourne BitStream(); // non atteint
    }
  }


/*ArrayXf pam_modulate(const BitStream &x)
{
  entier n = x.lon();
  ArrayXf y(n);
  pour(auto i = 0; i < n; i++)
    y(i) = x[i] ? 1.0f : -1.0f;
  retourne y;
}*/



  CmpBitsRes cmp_bits(const BitStream &bs1, const BitStream &bs2)
  {
    CmpBitsRes res;

    soit b1 = bs1.array(), b2 = bs2.array();

    soit b1c = b1.clone(), b2c = b2.clone();
    soit [b1a, b2a, dt, score] = tsd::fourier::aligne_entier(b1c, b2c);

    soit b1ar = real(b1a), b2ar = real(b2a);

    soit b1ap = b1ar.clone(), b2ap = b2ar.clone();


    // TODO : utiliser plutot segment dans le comptage d'erreurs
    // TODO : 100 ???
    si(b1a.rows() > 100)
    {
      soit n = b1a.rows();
      b1ap = b1ar.tail(n - 50).head(n-100).eval();
      b2ap = b2ar.tail(n - 50).head(n-100).eval();
    }

    res.nerr = 0;
    pour(auto i = 0; i < b1ap.rows(); i++)
    {
      si(round(b1ap(i)) != round(b2ap(i)))
        res.nerr++;
    }

    //res.nerr = (b1ap-b2ap).abs().sum();

    si(b2a.rows() == 0)
      res.ber = NAN;
    sinon
      res.ber = ((float) res.nerr) / b2ap.rows();

    msg("Cmp bits : dec = {}, nerr = {}, rows = {}, ber = {:.1e}", dt, res.nerr, b2ap.rows(), res.ber);

    res.decalage  = dt;
    res.score = score;
    res.b0 = b1ar;
    res.b1 = b2ar;

    si(0)
    {
      Figures f;
      f.subplot().plot(b1,              "hb", "Sequence d'entree ({} bits)", b1.rows());
      f.subplot().plot(res.b0,          "hb", "Sequence d'entree recalée");
      f.subplot().plot(b2,              "hb", "Sequence de sortie ({} bits)", b2.rows());
      f.subplot().plot(res.b1,          "hb", "Sequence de sortie recalée");
      f.subplot().plot(res.b1 - res.b0, "hr", "Erreurs ({}, ber = {:.1e})", res.nerr, res.ber);
      f.afficher(sformat("Comparaison bits (dec={}, s={})", dt,score));
    }


    retourne res;
  }


  CmpBitsRes cmp_bits_psk(const BitStream &bs0, const BitStream &bs1, entier k)
  {
    entier M = (1 << k);
    CmpBitsRes r[M];
    CmpBitsRes &br = r[0];
    br.b1 = bs1.array();

    // Essai éventuellement plusieurs phases possibles
    pour(auto l = 0; l < M; l++)
    {
      soit idx = symmap_binaire(bs1, k);

      pour(auto i = 0; i < idx.rows(); i++)
        idx(i) = ((entier) idx(i) + l) % M;

      BitStream bs;
      symdemap_binaire(bs, idx, k);
      //ArrayXf b1p = bs.vers_array();
      //ArrayXf b1p = wf->decode_symboles(symbs);

      msg("cmp bits...");
      r[l] = cmp_bits(bs0, bs);
      msg("  ...cmp bits ok.");

      msg(" test phase = {} : score = {} (nerr = {})", l, r[l].score, r[l].nerr);

      si((l == 0) || (r[l].nerr < br.nerr))
      {
        br           = r[l];
        br.dec_phase = l;
      }
      //infos("Essai dec phase = %d : corr = %.2f, nerr = %d", k, r[k].score, r[k].nerr);
    }
    //br.b0 = b0;
    retourne br;
  }


void diff_encode(BitStream &y, const BitStream &x)
{
  entier n = x.lon(), lb = 0;

  pour(auto i = 0; i < n; i++)
  {
    entier bi = x[i];
    entier bo = lb ^ bi;
    y.push(bo);
    lb = bo;
  }
}

void diff_decode(BitStream &y, const BitStream &x)
{
  entier n = x.lon();

  si(n == 0)
    retourne;

  entier lb = x[0];

  pour(auto i = 0; i < n-1; i++)
  {
    entier bi = x[i];
    entier bo = lb ^ bi;
    y.push(bo);
    lb = bi;
  }
}

void decode_hard(BitStream &y, const Vecf &llr)
{
  soit n = llr.rows();
  pour(auto i = 0; i < n; i++)
    y.push(llr(i) > 0 ? 1 : 0);
}

Veccf bruit_awgn(const Veccf &x, float σ)
{
  soit n = x.rows();
  Veccf y(n);
  y.set_real(real(x) + σ * randn(n));
  y.set_imag(imag(x) + σ * randn(n));
  retourne y;
}

Vecf bruit_awgn(const Vecf &x, float σ)
{
  retourne x + σ * randn(x.rows());
}



}



