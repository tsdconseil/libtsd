#include "tsd/tsd-all.hpp"

namespace tsd::filtrage {

float rifamp_impl(const Vecf &h, float ω, int type)
{
  double s = 0;
  soit n = h.rows();

  si((type == 1) ou (type == 2))
    pour(auto k = 0; k < n/2; k++)
      s += h(k) * 2 * cos(ω * (0.5 * (n-1) - k));
  sinon
    pour(auto k = 0; k < n/2; k++)
      s += h(k) * 2 * sin(ω * (0.5 * (n-1) - k));

  if(type == 1)
    s += h(n/2);

  retourne s;
}


int riftype(const Vecf &h)
{
  soit n = h.rows();
  soit pair = est_pair(n);

  // Symétrique
  si(h(0) == h(n-1))
    retourne pair ? 2 : 1;
  sinon
    retourne pair ? 4 : 3;
}


Vecf rifamp2(const Vecf &h, const Vecf &fr)
{
  int L = fr.rows();
  Vecf A(L);
  soit type = riftype(h);

  pour(auto i = 0; i < L; i++)
    A(i) = rifamp_impl(h, fr(i) * 2 * π, type);

  retourne A;
}

template<typename T>
Vecteur<T> rifamp_naif(const Vecteur<T> &h, entier L)
{
  Vecteur<T> A(L);
  soit type = riftype(h);

  pour(auto i = 0; i < L; i++)
    A(i) = rifamp_impl(h, (i * π) / L, type);

  retourne A;

# if 0
  soit n = h.rows();
  soit M = (n - 1) / 2;
  //assertion_msg((n & 1) == 1, "rifamp_naif : nombre de coefficients impair attendu (n = {})", n);


  auto eval = &amp_type_II;

  // Type I
  si(est_impair(n) && (h(0) == h(n-1)))
    retourne


  Vecteur<T> A(L);
  pour(auto i = 0; i < L; i++)
  {
    soit ω = (i * π) / L;
    // Coefficient central
    A(i) = h(M);
    pour(auto k = 0; k < M; k++)
      A(i) += 2 * h(k) * cos((M-k) * ω);
  }
  retourne A;
# endif
}


tuple<Vecf, Vecf> rifampn(const Vecf &h, entier L, bouléen symetrique)
{
  soit fr = linspace(0, 0.5 - (0.5 / L), L);
  retourne {fr, rifamp_naif(h, L)};
}

// TODO : template
tuple<Vecf, Vecf> rifamp(const Vecf &h, entier L, bouléen symetrique)
{
  // si N (nb coefs) > 2 * L (nb points de sortie désirés) :
  //  comment faire ?
  soit N = h.rows();
  soit M = (N-1) / 2;
  soit dim_fft = max(2*L, N);

  soit h2 = Veccf::zeros(dim_fft);
  h2.head(N) = h;

  soit H = fft(h2);

  // select [0,π]
  H = H.head(dim_fft/2).eval();
  // frequency grid
  soit w = linspace(0,dim_fft/2-1,dim_fft/2) * π / (dim_fft / 2);

  // Corrige le facteur de normalisation de la fonction FFT
  H *= sqrt((float)dim_fft);
  H *= polar(M * w); // type I ou II

  // type III ou IV
  si(!symetrique)
    H *= -std::complex<float>(0,1);

  assertion(H.rows() == L);

  soit fr = linspace(0, 0.5 - (0.5 / L), L);

  retourne {fr, real(H)}; // discard zero imaginary part
}






Vecf lsi(const Tabf &A, const Vecf &b, const Vecf &w, bouléen debug = non)
{
  soit n = A.rows();
  assertion_msg(n == b.rows(), "lsi: ");

  soit Ai = A;
  soit bi = b,
       x  = A.lsq(b),
       p = Vecf::ones(n);

  si(debug)
  {
    Figures f;
    f.subplot().plot(b, "g-", "Vecteur cible");
    f.subplot().plot(w, "a-", "Pondération");
    f.afficher("LSI - Vecteur cible");
  }



  pour(auto itr = 0; itr < 100; itr++)
  {
    soit e    = A.matprod(x) - b,
         Ea   = abs(e) * w;
    soit rms  = e.rms(),
         emax = (float) Ea.maxCoeff(),
         emin = (float) Ea.minCoeff();

    //msg("Vecteur d'errreur : {}", e);

    si(debug && ((itr % 10) == 0))
    {
      msg("Itération {} : rms = {}, norme inf = {}, emin = {}", itr, rms, emax, emin);
      Figures f;
      f.subplot().plot(A.matprod(x), "-g", "Ax");
      f.subplot().plot(Ea, "-r", "Erreur");
      f.subplot().plot(p, "-b", "Pondération L2");
      f.afficher(sformat("LSI - ITR {}",  itr));
    }

    si((emax - emin) / emax < 1e-5)
    {
      msg("emax - emin = {} -> exit.", emax - emin);
      break;
    }

    soit α = 0.5;
    p = (p * (1 + (α / (emax - emin)) * (Ea - emin)));
    p /= p.somme();

    // TODO : faire plus efficace !!!
    Ai = Tabf::diagonal(p).matprod(A);
    bi = b * p;

    x = Ai.lsq(bi);
  }
  msg("Fin lsi.");
  retourne x;
}





Vecf design_rif_eq(entier nc, const Vecf &D, const Vecf &W)
{
  soit n = D.rows(),
       m = (nc+1) / 2;

  msg("Design rif eq : nc = {}, n = {}, m = {}", nc, n, m);

  assertion(W.rows() == n);

  Tabf A(n, m);

  pour(auto i = 0; i < n; i++)
  {
    soit Ω = (i * π) / n;

    // TODO: utiliser rifamp

    pour(auto k = 0; k < m; k++)
    {
      si((k == m-1) && (nc & 1))
        A(i, k) = 1;
      sinon
      // Type I (nc impair)
      //si((nc % 2))
        A(i, k) = 2 * cos(Ω * (k - (nc-1)/2.0));
      //sinon
        // Type II (nc pair)
        //A(i, k) = 2 * cos(Omega * (k - (nc-1)/2.0f));
    }
  }


  soit x = lsi(A, D, W, debug_design);

  Vecf h(nc);
  h.head(nc/2) = x.head(nc/2);
  h.tail(nc/2) = x.head(nc/2).reverse();
  // Filtre de type I, coefficient central
  si(nc & 1)
    h(nc/2) = x(nc/2);

  msg("ok.");

  retourne h;
}

Vecf design_rif_eq(entier nc, const vector<SpecFreqIntervalle> &spec)
{
  soit ε = 1e-10f;
  soit l = 16 * nc;
  // Poids nul par défaut
  soit W = Vecf::zeros(l) + ε,
       D = Vecf::ones(l);

  msg("rifeq : nb intervalles de fréq : {}", nc);

  pour(auto i = 0u; i < spec.size(); i++)
  {
    soit &s = spec[i];
    soit ib = min((entier) round(2 * s.fb * l), l-1),
         ih = min((entier) round(2 * s.fh * l), l-1);

    msg("Intervalle : [{} ({}) - {} ({})], atten = {}, poids = {}", s.fb, ib, s.fh, ih, s.atten, s.poids);

    D.segment(ib, (ih - ib) + 1).setConstant(s.atten);
    W.segment(ib, (ih - ib) + 1).setConstant(s.poids + ε);

    si(i + 1 < spec.size())
    {
      soit ib2 = (entier) round(2 * spec[i+1].fb * l);
      soit cnt = (ib2 - ih) + 1;
      D.segment(ih, cnt) = linspace(s.atten, spec[i+1].atten, cnt);
    }
  }
  retourne design_rif_eq(nc, D, W);
}




}






