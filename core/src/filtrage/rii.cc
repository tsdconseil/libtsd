#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"

using namespace tsd::vue;
using namespace std;

namespace tsd::filtrage {



// wa = 1 <=> 2 tan(wd/2) = 1
// tan(wd/2) = 0.5
// wd = 0.92729..
// nuc = π * wd = 2.91

// pour la transformée bilinéaire
// wd : omega cut digital
// wa : omega cut analogique
float ωd_vers_ωa(float wd, float fe)
{
  retourne 2 * fe * tan(wd / (2 * fe));
}

float ωa_vers_ωd(float wa, float fe)
{
  retourne 2 * fe * atan(wa / (2 * fe));
}

float fd_vers_fa(float fd)
{
  retourne tan(π  * fd) / π;
}

float fa_vers_fd(float fa)
{
  retourne atan(π * fa) / π;
}

// Conversion fonction de transfert analogique vers numérique
FRat<cfloat> trf_bilineaire(const FRat<cfloat> &ha, float fe)
{
  FRat<cfloat> hz;

  soit nz = ha.numer.coefs.rows(),
       np = ha.denom.coefs.rows();
  soit [z,p] = ha.roots();

  hz.numer = Poly<cfloat>::from_roots((z + 2*fe) / (- z + 2*fe));
  hz.numer.mlt = ha.numer.mlt;
  hz.denom = Poly<cfloat>::from_roots((p + 2*fe) / (- p + 2*fe));
  hz.denom.mlt = ha.denom.mlt;


  si(np != nz)
  {
    soit p = Poly<cfloat>::from_roots(-Veccf::ones(abs(np - nz)));
    si(nz < np)
      hz.numer *= p;
    sinon
      hz.denom *= p;
  }

  cfloat gain = 1.0f;

  pour(auto i = 0; i < nz; i++)
    gain *= (2 * fe - ha.numer.coefs(i));
  pour(auto i = 0; i < np; i++)
    gain /= (2 * fe - ha.denom.coefs(i));

  hz.numer.mlt *= gain;
  retourne hz;
}

#if 0
// Conversion fonction de transfert analogique vers numérique
FRat<cfloat> trf_bilineaire_tustin(const FRat<cfloat> &ha, float K)
{
  float fe = K / 2;
  // On pourrait faire ceci :
  //  retourne ha.horner(2 * fe * (1 - z^-1) / (1 + z^-1));
  // mais cela ne preserve pas la forme racines.
  //soit n = pa.rows();
  FRat<cfloat> hz;

  //fe /= 2;

  soit nz = ha.numer.coefs.rows(), np = ha.denom.coefs.rows();



  hz.numer = Poly<cfloat>::from_roots((2*fe + ha.numer.roots()) / (2*fe - ha.numer.roots()));
  hz.numer.mlt = ha.numer.mlt;
  hz.denom = Poly<cfloat>::from_roots((2*fe + ha.denom.roots()) / (2*fe - ha.denom.roots()));
  hz.denom.mlt = ha.denom.mlt;


  si(np != nz)
  {
    auto p = Poly<cfloat>::from_roots(-ArrayXf::Ones(abs(np - nz)));
    si(nz < np)
      hz.numer *= p;
    sinon
      hz.denom *= p;
  }

  cfloat gain = 1.0f;

  pour(auto i = 0; i < ha.numer.coefs.rows(); i++)
    gain *= (2 * fe - ha.numer.coefs(i));
  pour(auto i = 0; i < ha.denom.coefs.rows(); i++)
    gain /= (2 * fe - ha.denom.coefs(i));

  //pour(auto i = 0; i < ha.numer.coefs.rows(); i++)
  //  gain *= (1.0f - hz.denom.coefs(i)) / (1.0f - hz.numer.coefs(i));

  hz.numer.mlt *= gain;
  retourne hz;
}
#endif

// On est sur le cercle unité -> donc inversion <=> interversion avec le conjugé -> sans effet!

// Passe-bas analogique normalisé vers passe-haut analogique
FRat<cfloat> pban_vers_pha(const FRat<cfloat> &ha, float ωc)
{
  tsd_assert(ha.numer.mode_racines && ha.denom.mode_racines);
  FRat<cfloat> res;

  soit [N,M] = ha.degrés();

  res.numer = Poly<cfloat>::from_roots(ha.numer.coefs.inv() * ωc);
  res.denom = Poly<cfloat>::from_roots(ha.denom.coefs.inv() * ωc);

  res.numer.mlt = ha.numer.mlt * ha.numer.coefs.produit();
  res.denom.mlt = ha.denom.mlt * ha.denom.coefs.produit();

  si(M > N)
    res.numer.coefs = vconcat(res.numer.coefs, Veccf::zeros(M - N));
  sinon si(M < N)
    res.denom.coefs = vconcat(res.denom.coefs, Veccf::zeros(N - M));

  retourne res;
}

// passe bas normalisé vers passe-bande
FRat<cfloat> pban_vers_pbda(const FRat<cfloat> &ha, float ωc, float Q)
{
  // https://en.wikipedia.org/wiki/Prototype_filter
  // w -> Q (s/wc + wc/s)

  tsd_assert(ha.numer.mode_racines && ha.denom.mode_racines);

  FRat<cfloat> res;// = ha;

  soit [N,M] = ha.degrés();

  res.numer = Poly<cfloat>::from_roots(ha.numer.coefs.inv() * ωc);
  res.denom = Poly<cfloat>::from_roots(ha.denom.coefs.inv() * ωc);

  res.numer.mlt = ha.numer.mlt * res.numer.coefs.produit();
  res.denom.mlt = ha.denom.mlt * res.denom.coefs.produit();

  si(M > N)
    res.numer.coefs = vconcat(res.numer.coefs, Veccf::zeros(M - N));
  sinon si(M < N)
    res.denom.coefs = vconcat(res.denom.coefs, Veccf::zeros(N - M));

  retourne res;
}

// Passe-bas analogique normalisé vers passe-bas analogique
FRat<cfloat> pban_vers_pba(const FRat<cfloat> &ha, float ωc)
{
  tsd_assert(ha.numer.mode_racines && ha.denom.mode_racines);
  FRat<cfloat> res;

  res.numer = Poly<cfloat>::from_roots(ha.numer.coefs * ωc);
  res.denom = Poly<cfloat>::from_roots(ha.denom.coefs * ωc);
  res.numer.mlt = ha.numer.mlt;
  res.denom.mlt = ha.denom.mlt;

  soit [M, N] = ha.degrés();

  res.numer.mlt *= pow(ωc, N - M);
  retourne res;
}







FRat<cfloat> butterworth_analogique(entier n)
{
  // Pôles analogiques, wc = wcut
  /*ArrayXcf butterworth_poles(unsigned int n)
  {
    dsp_assert(n > 0,        "Butterworh coef gen: order must > 0.");
    // Formule d'après Wikipédia
    auto k = linspace(1, n, n);
    retourne polar(π*(2*k+n-1)/(2*n));
  }*/

  // Formule d'après Wikipédia
  soit k   = linspace(1, n, n);
  soit pol = polar(π*(2*k+n-1)/(2*n));
  //soit pol = butterworth_poles(n);

  FRat<cfloat> ha;
  ha.numer.mode_racines = oui;
  ha.numer.coefs.resize(0);
  ha.denom = Poly<cfloat>::from_roots(pol);

  retourne ha;
}

// "Elliptic Functions pour Filter Design"
// H J Orchard and Alan N Willson  IEE Trans on Circuits and Systems April 97
FRat<cfloat> elliptique_analogique(entier n, float rp, float rs)
{
  vector<cfloat> poles, zeros;
  FRat<cfloat> ha;

  Veccf hp, hz;

  msg("Design elliptique : n = {}, rp = {}, rs = {}", n, rp, rs);

  si(n == 1)
  {
    hp.resize(1);
    hp(0) = -sqrt(1.0 / db2pow(rp) - 1);
  }
  sinon
  {
    soit dbn = log(10.0) / 20;
    soit apn = dbn * rp, asn = dbn * rs;

    vector<double> g;
    soit e0 = sqrt(2*exp(apn)*sinh(apn));
    g.push_back(e0 / sqrt(exp(2*asn)-1));
    soit v = g[0];
    tantque(v > 1e-150)
    {
      v = carré(v / (1 + sqrt(1-v*v)));
      g.push_back(v);
    }

    soit eg = Vecd::map(g);

    // m2 = dim de g et de e
    soit m2 = (entier) g.size();

    soit ek = Vecd::zeros(m2 + 11);
    soit m1 = m2;
    pour(auto index = 0; index <= 10; index++)
    {
      m1 = m2 + index;
      ek(m1-1) = 4 * pow(g[m2-1]/4, (1.0 * (1 << index)) / n);
      si(ek(m1-1) < 1e-14)
        break;
    }
    pour(auto i = m1 - 1; i >= 1; i--)
      ek(i-1) = 2 * sqrt(ek(i)) / (1 + ek(i));

    //msg("ek = {}", ek);

    Vecd e(m2);
    e(0) = e0;
    // Calcul des poles et des zéros
    pour(auto i = 1; i < m2; i++)
    {
      soit a = (1 + g[i]) * e(i-1) / 2;
      e(i) = a + sqrt(a*a + g[i]);
    }

    msg("e = {}", e);

    soit elast = e(m2-1);
    soit u2 = log((1 + sqrt(1 + elast * elast))/elast)/n;

    msg("u2 = {}", u2);

    // Cosinus complexe
    soit cosc = [](cdouble x) -> cdouble
    {
      retourne cdouble{cos(x.real()) * cosh(x.imag()), -sin(x.real()) * sinh(x.imag())};
    };

    pour(auto i = 0; i < n / 2; i++)
    {
      double u1 = (2 * i + 1) * π / (2 * n);
      cdouble c = cdouble{0.0, -1.0} / cosc(cdouble{-u1, u2});
      double d = 1 / cos(u1);

      // m1+1 < dim(ek)
      pour(auto j = m1; j >= 2; j--)
      {
        c = (c - ek(j-1)/c)/(1 + ek(j-1));
        d = (d + ek(j-1)/d)/(1 + ek(j-1));
      }
      poles.push_back((cfloat) (1.0 / c));
      poles.push_back((cfloat) (conj(1.0 / c)));
      zeros.push_back(cfloat(0.0, d/ek(0)));
      zeros.push_back(cfloat(0.0, -d/ek(0)));
    }


    si(n & 1)
    {
      double a = 1 / sinh(u2);
      pour(auto j = m1; j >= 2; j--)
        a = (a - ek(j)/a) / (1.0 + ek(j));
      poles.push_back((cfloat) (-1.0f / a));
    }

    hz = Veccf::map(zeros);
    hp = Veccf::map(poles);



    //si((n & 1) == 0)
    //  k /= sqrt(1 + ε)
  }

  ha.numer.mode_racines = oui;
  ha.numer = Poly<cfloat>::from_roots(hz);
  ha.denom = Poly<cfloat>::from_roots(hp);
  ha.numer.mlt = real((-hp).produit() / (-hz).produit());
  retourne ha;
}






FRat<cfloat> tchebychev_I_analogique(entier n, float δ)
{
  soit m = linspace(1, n, n);
  soit θ = (2*m - 1) * (π / (2*n));

  // Formule d'après Wikipédia
  soit ε   = sqrt(db2pow(δ) - 1);
  soit ash = asinh(1.0f/ε)/n;
  soit s   = sinh(ash),
       c   = cosh(ash);

  soit sc = s * signe(s);

  Veccf poles(n);
  pour(auto k = 0; k < n; k++)
  {
    poles(k) = - sin(θ(k)) * signe(sin(θ(k))) * sc
               + cos(θ(k)) * 1.0if * c;
  }


  FRat<cfloat> ha;
  ha.numer.mode_racines = oui;
  ha.numer.coefs.resize(0);
  ha.denom = Poly<cfloat>::from_roots(poles);

  // Etant donné une fonction de transfert analogique,
  // comment forcé le gain DC à être unitaire ?
  // Gain DC = ha(0) = produit des racines
  ha.numer.mlt *= ha.denom.horner<cfloat>(0);

  retourne ha;
}

FRat<cfloat> tchebychev_II_analogique(unsigned int n, float δ)
{
  //soit pol = tchebychev_II_poles(n, δ);

  soit m = linspace(1, n, n);
  soit θ = (2.0f*m - 1.0f) * (π/(2*n));

  // Formule d'après Wikipédia
  soit ε = 1.0f / sqrt(db2pow(δ) - 1);

  msg("techbychev II : n={}, δ={:.2f} dB, ε={:e}", n, δ, ε);

  soit ash = asinh(1.0f/ε)/n;
  soit s   = sinh(ash);
  soit c   = cosh(ash);

  Veccf pôles(n);

  pour(auto i = 0u; i < n; i++)
    pôles(i) = 1.0f / (- s * sin(θ(i)) * copysign(1.0f, s) * copysign(1.0f, sin(θ(i)))
                 + 1.0if * c * cos(θ(i)));

  Veccf zéros(n);
  pour(auto i = 0u; i < n; i++)
    zéros(i) = 1.0f / (- 1.0if * cos((float) ((π/2) * (2*(i+1.0f)-1)/n)));

  FRat<cfloat> ha;
  ha.numer = Poly<cfloat>::from_roots(zéros);
  ha.denom = Poly<cfloat>::from_roots(pôles);
  ha.numer.mlt *= ha.denom.horner<cfloat>(0) / ha.numer.horner<cfloat>(0);

  retourne ha;
}

FRat<cfloat> design_riia_laplace(entier n, TypeFiltre type, PrototypeAnalogique prototype,
    float fcut, float δ_bp, float δ_bc)
{
  soit ωa = ωd_vers_ωa(2*π*fcut, 1);

  msg("WD = {} -> WA = {}", 2*π*fcut, ωa);

  FRat<cfloat> ha;

  si(prototype == PrototypeAnalogique::BUTTERWORTH)
    ha = butterworth_analogique(n);
  sinon si(prototype == PrototypeAnalogique::TCHEBYCHEV_I)
    ha = tchebychev_I_analogique(n, δ_bp);
  sinon si(prototype == PrototypeAnalogique::TCHEBYCHEV_II)
    ha = tchebychev_II_analogique(n, δ_bc);
  sinon si(prototype == PrototypeAnalogique::ELLIPTIQUE)
    ha = elliptique_analogique(n, δ_bp, δ_bc);
  sinon
  {
    msg_erreur("TODO");
  }

  si(type == TypeFiltre::PASSE_BAS)
  {
    msg("design_riia_laplace : passe bas");
    ha = pban_vers_pba(ha, ωa);
  }
  sinon si(type == TypeFiltre::PASSE_HAUT)
  {
    msg("design_riia_laplace : passe haut");
    ha = pban_vers_pha(ha, ωa);
  }
  sinon
  {
    msg_erreur("TODO : design_riia PB / SB.");
  }

  retourne ha;
}

// K = wo / tan(w0/(2fe)) à la place de 2fe

FRat<cfloat> design_riia(entier n, TypeFiltre type, PrototypeAnalogique prototype, float fcut, float δ_bp, float δ_bc)
{
  soit ha = design_riia_laplace(n, type, prototype, fcut, δ_bp, δ_bc);
  retourne trf_bilineaire(ha, 1.0f);
}

PrototypeAnalogique parse_proto(const string &s)
{
  si(s.substr(0, 1) == "b")
    retourne PrototypeAnalogique::BUTTERWORTH;
  sinon si(s == "cheb1")
    retourne PrototypeAnalogique::TCHEBYCHEV_I;
  sinon si(s == "cheb2")
    retourne PrototypeAnalogique::TCHEBYCHEV_II;
  sinon si(s.substr(0, 5) == "ellip")
    retourne PrototypeAnalogique::ELLIPTIQUE;
  msg_erreur("Prototype analogique inconnu : {}", s);
  retourne PrototypeAnalogique::BUTTERWORTH;
}

TypeFiltre parse_tf(const string &s)
{
  si((s == "lp") || (s == "pb"))
    retourne TypeFiltre::PASSE_BAS;
  sinon si((s == "hp") || (s == "ph"))
    retourne TypeFiltre::PASSE_HAUT;
  sinon si((s == "pb") || (s == "bp"))
    retourne TypeFiltre::PASSE_BANDE;
  sinon si(s == "sb")
    retourne TypeFiltre::COUPE_BANDE;
  msg_erreur("Type filtre inconnu : {}", s);
  retourne TypeFiltre::PASSE_BAS;
}


FRat<cfloat> design_riia(entier n, const string &type,
    const string &prototype, float fcut, float δ_bp, float δ_bc)
{
  retourne design_riia(n, parse_tf(type), parse_proto(prototype), fcut, δ_bp, δ_bc);
}

FRat<float> design_biquad(const std::string type, float f, float Q, float gain_dB)
{
  BiquadSpec spec{.f = f, .Q = Q, .gain_dB = gain_dB};
  si((type == "lp") || (type == "pb"))
    spec.type = BiquadSpec::PASSE_BAS;
  sinon si((type == "hp") || (type == "ph"))
    spec.type = BiquadSpec::PASSE_HAUT;
  sinon si((type == "bp") || (type == "passe-bande"))
    spec.type = BiquadSpec::PASSE_BANDE;
  sinon si((type == "cb") || (type == "notch") || (type == "sb"))
    spec.type = BiquadSpec::COUPE_BANDE;
  sinon si((type == "plateau-bf"))
    spec.type = BiquadSpec::PLATEAU_BF;
  sinon si((type == "plateau-hf"))
    spec.type = BiquadSpec::PLATEAU_HF;
  sinon si((type == "res"))
    spec.type = BiquadSpec::RESONATEUR;
  retourne design_biquad(spec);
}

std::ostream& operator<<(std::ostream &ss, const BiquadSpec &t)
{
  switch(t.type)
  {
  case BiquadSpec::PASSE_BAS:
    ss << fmt::format("Passe-bas(fc={}, Q={})", t.f, t.Q);
    break;
  case BiquadSpec::PASSE_HAUT:
    ss << fmt::format("Passe-haut(fc={}, Q={})", t.f, t.Q);
    break;
  case BiquadSpec::PASSE_BANDE:
    ss << fmt::format("Passe-bande(fc={}, Q={})", t.f, t.Q);
    break;
  case BiquadSpec::COUPE_BANDE:
    ss << fmt::format("Coupe-bande(fc={}, Q={})", t.f, t.Q);
    break;
  case BiquadSpec::RESONATEUR:
    ss << fmt::format("Résonateur(fc={}, Q={}, g={:.1f} dB)", t.f, t.Q, t.gain_dB);
    break;
  case BiquadSpec::PLATEAU_BF:
    ss << fmt::format("plateau-bf(fc={}, g={:.1f} dB)", t.f, t.gain_dB);
    break;
  case BiquadSpec::PLATEAU_HF:
    ss << fmt::format("plateau-hf(fc={}, g={:.1f} dB)", t.f, t.gain_dB);
    break;
  }
  retourne ss;
}

// D'après https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
FRat<float> design_biquad(const BiquadSpec &spec)
{
  // Gain en amplitude
  float A = sqrt(db2pow(spec.gain_dB));

  msg("design_biquad: {}", spec);
  msg("design_biquad: gain = {} dB -> A = {}", spec.gain_dB, A);

  // Pulsation
  float ω = 2 * π * spec.f;
  float sn = sin(ω), cs = cos(ω);

  // Facteur d'ammortissement
  float α = sn / (2 * spec.Q);
  //float α = sn * sinh(0.5f * log(2.0f) * bw * omega);
  float β = sqrt(2 * A);

  float a0 = 1, a1 = 1, a2 = 1, b0 = 1, b1 = 1, b2 = 1;

  switch(spec.type)
  {
    case BiquadSpec::PASSE_BAS:
    b0 = (1 - cs) / 2;
    b1 = 1 - cs;
    b2 = (1 - cs) / 2;
    a0 = 1 + α;
    a1 = -2 * cs;
    a2 = 1 - α;
    break;
  case BiquadSpec::PASSE_HAUT:
    b0 = (1 + cs) / 2;
    b1 = -(1 + cs);
    b2 = (1 + cs) / 2;
    a0 = 1 + α;
    a1 = -2 * cs;
    a2 = 1 - α;
    break;
  case BiquadSpec::COUPE_BANDE:
    b0 = α;
    b1 = 0;
    b2 = -α;
    a0 = 1 + α;
    a1 = -2 * cs;
    a2 = 1 - α;
    break;
  case BiquadSpec::PASSE_BANDE:
    b0 = 1;
    b1 = -2 * cs;
    b2 = 1;
    a0 = 1 + α;
    a1 = -2 * cs;
    a2 = 1 - α;
    break;
  case BiquadSpec::RESONATEUR:
    b0 = 1 + (α * A);
    b1 = -2 * cs;
    b2 = 1 - (α * A);
    a0 = 1 + (α / A);
    a1 = -2 * cs;
    a2 = 1 - (α / A);
    break;
  case BiquadSpec::PLATEAU_BF:
    b0 = A * ((A + 1) - (A - 1) * cs + β * sn);
    b1 = 2 * A * ((A - 1) - (A + 1) * cs);
    b2 = A * ((A + 1) - (A - 1) * cs - β * sn);
    a0 = (A + 1) + (A - 1) * cs + β * sn;
    a1 = -2 * ((A - 1) + (A + 1) * cs);
    a2 = (A + 1) + (A - 1) * cs - β * sn;
    break;
  case BiquadSpec::PLATEAU_HF:
    b0 = A * ((A + 1) + (A - 1) * cs + β * sn);
    b1 = -2 * A * ((A - 1) + (A + 1) * cs);
    b2 = A * ((A + 1) + (A - 1) * cs - β * sn);
    a0 = (A + 1) - (A - 1) * cs + β * sn;
    a1 = 2 * ((A - 1) - (A + 1) * cs);
    a2 = (A + 1) - (A - 1) * cs - β * sn;
    break;
  default:
    echec("Type de biquad invalide.");
  }

  a1 /= a0;
  a2 /= a0;
  b0 /= a0;
  b1 /= a0;
  b2 /= a0;
  a0 = 1;

  soit z = (FRat<float>::z()).inv();

  retourne (b0 + b1 * z + b2 * z.pow(2)) / (1 + a1 * z + a2 * z.pow(2));
}




}

