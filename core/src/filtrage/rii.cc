#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"

using namespace tsd::vue;
using namespace std;

namespace tsd::filtrage {


// Pb : si définition par coefficient, dérive
/*Poly<float> bessel_polynome(unsigned int n)
{
  ArrayXf res;


  return Poly<float>(res);
}


FRat<float> bessel(unsigned int n)
{
  FRat<float> b;

  if(n == 0)
    return

}*/




// wa = 1 <=> 2 tan(wd/2) = 1
// tan(wd/2) = 0.5
// wd = 0.92729..
// nuc = π * wd = 2.91

// pour la transformée bilinéaire
// wd : omega cut digital
// wa : omega cut analogique
float ωd_vers_ωa(float wd, float fe)
{
  return 2 * fe * tan(wd / (2 * fe));
}

float ωa_vers_ωd(float wa, float fe)
{
  return 2 * fe * atan(wa / (2 * fe));
}

float fd_vers_fa(float fd)
{
  return tan(π  * fd) / π;
}

float fa_vers_fd(float fa)
{
  return atan(π * fa) / π;
}

// Conversion fonction de transfert analogique vers numérique
FRat<cfloat> trf_bilineaire(const FRat<cfloat> &ha, float fe)
{
  // On pourrait faire ceci :
  //  return ha.horner(2 * fe * (1 - z^-1) / (1 + z^-1));
  // mais cela ne preserve pas la forme racines.
  //auto n = pa.rows();
  FRat<cfloat> hz;

  //fe /= 2;

  auto nz = ha.numer.coefs.rows(), np = ha.denom.coefs.rows();



  hz.numer = Poly<cfloat>::from_roots((2*fe + ha.numer.roots()) / (2*fe - ha.numer.roots()));
  hz.numer.mlt = ha.numer.mlt;
  hz.denom = Poly<cfloat>::from_roots((2*fe + ha.denom.roots()) / (2*fe - ha.denom.roots()));
  hz.denom.mlt = ha.denom.mlt;


  if(np != nz)
  {
    auto p = Poly<cfloat>::from_roots(-ArrayXf::Ones(abs(np - nz)));
    if(nz < np)
      hz.numer = hz.numer * p;
    else
      hz.denom = hz.denom * p;
  }

  cfloat gain = 1.0f;

  for(auto i = 0; i < ha.numer.coefs.rows(); i++)
    gain *= (2 * fe - ha.numer.coefs(i));
  for(auto i = 0; i < ha.denom.coefs.rows(); i++)
    gain /= (2 * fe - ha.denom.coefs(i));

  //for(auto i = 0; i < ha.numer.coefs.rows(); i++)
  //  gain *= (1.0f - hz.denom.coefs(i)) / (1.0f - hz.numer.coefs(i));

  hz.numer.mlt *= gain;
  return hz;
}

#if 0
// Conversion fonction de transfert analogique vers numérique
FRat<cfloat> trf_bilineaire_tustin(const FRat<cfloat> &ha, float K)
{
  float fe = K / 2;
  // On pourrait faire ceci :
  //  return ha.horner(2 * fe * (1 - z^-1) / (1 + z^-1));
  // mais cela ne preserve pas la forme racines.
  //auto n = pa.rows();
  FRat<cfloat> hz;

  //fe /= 2;

  auto nz = ha.numer.coefs.rows(), np = ha.denom.coefs.rows();



  hz.numer = Poly<cfloat>::from_roots((2*fe + ha.numer.roots()) / (2*fe - ha.numer.roots()));
  hz.numer.mlt = ha.numer.mlt;
  hz.denom = Poly<cfloat>::from_roots((2*fe + ha.denom.roots()) / (2*fe - ha.denom.roots()));
  hz.denom.mlt = ha.denom.mlt;


  if(np != nz)
  {
    auto p = Poly<cfloat>::from_roots(-ArrayXf::Ones(abs(np - nz)));
    if(nz < np)
      hz.numer = hz.numer * p;
    else
      hz.denom = hz.denom * p;
  }

  cfloat gain = 1.0f;

  for(auto i = 0; i < ha.numer.coefs.rows(); i++)
    gain *= (2 * fe - ha.numer.coefs(i));
  for(auto i = 0; i < ha.denom.coefs.rows(); i++)
    gain /= (2 * fe - ha.denom.coefs(i));

  //for(auto i = 0; i < ha.numer.coefs.rows(); i++)
  //  gain *= (1.0f - hz.denom.coefs(i)) / (1.0f - hz.numer.coefs(i));

  hz.numer.mlt *= gain;
  return hz;
}
#endif

// On est sur le cercle unité -> donc inversion <=> interversion avec le conjugé -> sans effet!

// Passe-bas analogique normalisé vers passe-haut analogique
FRat<cfloat> pban_vers_pha(const FRat<cfloat> &ha, float ωc)
{
  tsd_assert(ha.numer.mode_racines && ha.denom.mode_racines);
  FRat<cfloat> res = ha;

  int N = ha.numer.coefs.rows(), M = ha.denom.coefs.rows();

  res.numer.mlt *= res.numer.coefs.prod();
  res.denom.mlt *= res.denom.coefs.prod();

  cfloat tmp = ωc;

  res.numer.coefs = tmp / ha.numer.coefs;
  res.denom.coefs = tmp / ha.denom.coefs;

  if(M > N)
    res.numer.coefs = vconcat(res.numer.coefs, ArrayXcf::Zero(M - N));
  else if(M < N)
    res.denom.coefs = vconcat(res.denom.coefs, ArrayXcf::Zero(N - M));

  return res;
}

// passe bas normalisé vers passe-bande
FRat<cfloat> pban_vers_pbda(const FRat<cfloat> &ha, float ωc, float Q)
{
  // https://en.wikipedia.org/wiki/Prototype_filter
  // w -> Q (s/wc + wc/s)

  //ha.horner(s)

  tsd_assert(ha.numer.mode_racines);
  tsd_assert(ha.denom.mode_racines);
  FRat<cfloat> res = ha;

  int N = ha.numer.coefs.rows(), M = ha.denom.coefs.rows();

  res.numer.mlt *= res.numer.coefs.prod();
  res.denom.mlt *= res.denom.coefs.prod();

  cfloat tmp = ωc;

  res.numer.coefs = tmp / ha.numer.coefs;
  res.denom.coefs = tmp / ha.denom.coefs;

  if(M > N)
    res.numer.coefs = vconcat(res.numer.coefs, ArrayXcf::Zero(M - N));
  else if(M < N)
    res.denom.coefs = vconcat(res.denom.coefs, ArrayXcf::Zero(N - M));

  return res;
}

// Passe-bas analogique normalisé vers passe-bas analogique
FRat<cfloat> pban_vers_pba(const FRat<cfloat> &ha, float ωc)
{
  tsd_assert(ha.numer.mode_racines);
  tsd_assert(ha.denom.mode_racines);
  FRat<cfloat> res = ha;
  res.numer.coefs *= ωc;
  res.denom.coefs *= ωc;

  int M = ha.numer.coefs.rows(), N = ha.denom.coefs.rows();

  res.numer.mlt *= pow(ωc, N - M);

  return res;
}







FRat<cfloat> butterworth_analogique(unsigned int n)
{
  // Pôles analogiques, wc = wcut
  /*ArrayXcf butterworth_poles(unsigned int n)
  {
    dsp_assert(n > 0,        "Butterworh coef gen: order must > 0.");
    // Formule d'après Wikipédia
    auto k = linspace(1, n, n);
    return polar(π*(2*k+n-1)/(2*n));
  }*/

  // Formule d'après Wikipédia
  auto k = linspace(1, n, n);
  auto pol = polar(π*(2*k+n-1)/(2*n));
  //auto pol = butterworth_poles(n);

  FRat<cfloat> ha;
  ha.numer.mode_racines = true;
  ha.numer.coefs.resize(0);
  ha.denom = Poly<cfloat>::from_roots(pol);

  return ha;
}

// "Elliptic Functions for Filter Design"
// H J Orchard and Alan N Willson  IEE Trans on Circuits and Systems April 97
FRat<cfloat> elliptique_analogique(int n, float rp, float rs)
{
  vector<cfloat> poles, zeros;
  FRat<cfloat> ha;

  ArrayXcf hp, hz;

  msg("Design elliptique : n = {}, rp = {}, rs = {}", n, rp, rs);

  if(n == 1)
  {
    hp.resize(1);
    hp(0) = -sqrt(1.0 / db2pow(rp) - 1);
  }
  else
  {
    double dbn = log(10.0) / 20;
    //int n3 = n/2;
    double apn = dbn * rp;
    double asn = dbn * rs;

    vector<double> g;
    double e0 = sqrt(2*exp(apn)*sinh(apn));
    g.push_back(e0 / sqrt(exp(2*asn)-1));
    double v = g[0];
    while(v > 1e-150)
    {
      v = carré(v / (1 + sqrt(1-v*v)));
      g.push_back(v);
    }

    Eigen::ArrayXd eg = Eigen::Map<Eigen::ArrayXd>(g.data(), g.size());

    // m2 = dim de g et de e
    int m2 = g.size();

    //msg("m2 = {}, g = {}", m2, eg);

    Eigen::ArrayXd ek = Eigen::ArrayXd::Zero(m2 + 11);
    int m1 = m2;
    for(auto index = 0; index <= 10; index++)
    {
      m1 = m2 + index;
      ek(m1-1) = 4 * pow(g[m2-1]/4, (1.0 * (1 << index)) / n);
      if(ek(m1-1) < 1e-14)
        break;
    }
    for(auto i = m1 - 1; i >= 1; i--)
      ek(i-1) = 2 * sqrt(ek(i)) / (1 + ek(i));

    //msg("ek = {}", ek);

    Eigen::ArrayXd e(m2);
    e(0) = e0;
    // Calcul des poles et des zéros
    for(auto i = 1; i < m2; i++)
    {
      double a = (1 + g[i]) * e(i-1) / 2;
      e(i) = a + sqrt(a*a + g[i]);
    }

    msg("e = {}", e);

    auto elast = e(m2-1);
    double u2 = log((1 + sqrt(1 + elast * elast))/elast)/n;

    msg("u2 = {}", u2);

    // Cosinus complexe
    auto cosc = [](cdouble x) -> cdouble
    {
      return cdouble{cos(x.real()) * cosh(x.imag()), -sin(x.real()) * sinh(x.imag())};
    };

    for(auto i = 0; i < n / 2; i++)
    {
      double u1 = (2 * i + 1) * π / (2 * n);
      cdouble c = cdouble{0.0, -1.0} / cosc(cdouble{-u1, u2});
      double d = 1 / cos(u1);

      // m1+1 < dim(ek)
      for(auto j = m1; j >= 2; j--)
      {
        c = (c - ek(j-1)/c)/(1 + ek(j-1));
        d = (d + ek(j-1)/d)/(1 + ek(j-1));
      }
      poles.push_back((cfloat) (1.0 / c));
      poles.push_back((cfloat) (conj(1.0 / c)));
      zeros.push_back(cfloat(0.0, d/ek(0)));
      zeros.push_back(cfloat(0.0, -d/ek(0)));
    }


    if(n & 1)
    {
      double a = 1 / sinh(u2);
      for(auto j = m1; j >= 2; j--)
        a = (a - ek(j)/a) / (1.0 + ek(j));
      poles.push_back((cfloat) (-1.0f / a));
    }

    hz = Eigen::Map<ArrayXcf>(zeros.data(), zeros.size());
    hp = Eigen::Map<ArrayXcf>(poles.data(), poles.size());



    //if((n & 1) == 0)
    //  k /= sqrt(1 + ε)
  }

  ha.numer.mode_racines = true;
  ha.numer = Poly<cfloat>::from_roots(hz);
  ha.denom = Poly<cfloat>::from_roots(hp);

  float k = real((-hp).prod() / (-hz).prod());
  ha.numer.mlt = k;
  return ha;
}






FRat<cfloat> tchebychev_I_analogique(unsigned int n, float δ)
{
  ArrayXf m = linspace(1, n, n);
  ArrayXf θ = (2.0f*m - 1.0f) * (π/(2*n));

  // Formule d'après Wikipédia
  auto ε = sqrt(db2pow(δ) - 1);
  auto ash = asinh(1.0f/ε)/n;
  auto s = sinh(ash);
  auto c = cosh(ash);

  ArrayXcf poles = - s * sin(θ) * copysign(1.0f, s) * sin(θ).sign()
        + 1.0if * c * cos(θ);

  FRat<cfloat> ha;
  ha.numer.mode_racines = true;
  ha.numer.coefs.resize(0);
  ha.denom = Poly<cfloat>::from_roots(poles);

  // Etant donné une fonction de transfert analogique,
  // comment forcé le gain DC à être unitaire ?
  // Gain DC = ha(0) = produit des racines

  ha.numer.mlt *= ha.denom.horner<cfloat>(0);

  return ha;
}

FRat<cfloat> tchebychev_II_analogique(unsigned int n, float δ)
{
  //auto pol = tchebychev_II_poles(n, δ);

  ArrayXf m = linspace(1, n, n);
  ArrayXf θ = (2.0f*m - 1.0f) * (π/(2*n));

  // Formule d'après Wikipédia
  auto ε = 1.0f / sqrt(db2pow(δ) - 1);

  msg("techbychev II : n={}, δ={:.2f} dB, ε={:e}", n, δ, ε);

  auto ash = asinh(1.0f/ε)/n;
  auto s   = sinh(ash);
  auto c   = cosh(ash);

  ArrayXcf pôles = 1.0f / (- s * sin(θ) * copysign(1.0f, s) * sin(θ).sign()
        + 1.0if * c * cos(θ));

  ArrayXcf zéros(n);
  for(auto i = 0u; i < n; i++)
    zéros(i) = 1.0f / (- 1.0if * cos((float) ((π/2) * (2*(i+1.0f)-1)/n)));

  FRat<cfloat> ha;
  ha.numer = Poly<cfloat>::from_roots(zéros);
  ha.denom = Poly<cfloat>::from_roots(pôles);
  ha.numer.mlt *= ha.denom.horner<cfloat>(0) / ha.numer.horner<cfloat>(0);

  return ha;
}

FRat<cfloat> design_riia_laplace(int n, TypeFiltre type, PrototypeAnalogique prototype,
    float fcut, float δ_bp, float δ_bc)
{
  auto ωa = ωd_vers_ωa(2*π*fcut, 1);

  msg("WD = {} -> WA = {}", 2*π*fcut, ωa);

  FRat<cfloat> ha;

  if(prototype == PrototypeAnalogique::BUTTERWORTH)
    ha = butterworth_analogique(n);
  else if(prototype == PrototypeAnalogique::TCHEBYCHEV_I)
    ha = tchebychev_I_analogique(n, δ_bp);
  else if(prototype == PrototypeAnalogique::TCHEBYCHEV_II)
    ha = tchebychev_II_analogique(n, δ_bc);
  else if(prototype == PrototypeAnalogique::ELLIPTIQUE)
    ha = elliptique_analogique(n, δ_bp, δ_bc);
  else
  {
    msg_erreur("TODO");
  }

  if(type == TypeFiltre::PASSE_BAS)
  {
    msg("design_riia_laplace : passe bas");
    ha = pban_vers_pba(ha, ωa);
  }
  else if(type == TypeFiltre::PASSE_HAUT)
  {
    msg("design_riia_laplace : passe haut");
    ha = pban_vers_pha(ha, ωa);
  }
  else
  {
    msg_erreur("TODO : design_riia PB / SB.");
  }

  return ha;
}

// K = wo / tan(w0/(2fe)) à la place de 2fe

FRat<cfloat> design_riia(int n, TypeFiltre type, PrototypeAnalogique prototype, float fcut, float δ_bp, float δ_bc)
{
  auto ha = design_riia_laplace(n, type, prototype, fcut, δ_bp, δ_bc);
  return trf_bilineaire(ha, 1.0f);
}

PrototypeAnalogique parse_proto(const string &s)
{
  if(s.substr(0, 1) == "b")
    return PrototypeAnalogique::BUTTERWORTH;
  else if(s == "cheb1")
    return PrototypeAnalogique::TCHEBYCHEV_I;
  else if(s == "cheb2")
    return PrototypeAnalogique::TCHEBYCHEV_II;
  else if(s.substr(0, 5) == "ellip")
    return PrototypeAnalogique::ELLIPTIQUE;
  msg_erreur("Prototype analogique inconnu : {}", s);
  return PrototypeAnalogique::BUTTERWORTH;
}

TypeFiltre parse_tf(const string &s)
{
  if((s == "lp") || (s == "pb"))
    return TypeFiltre::PASSE_BAS;
  else if((s == "hp") || (s == "ph"))
    return TypeFiltre::PASSE_HAUT;
  else if((s == "pb") || (s == "bp"))
    return TypeFiltre::PASSE_BANDE;
  else if(s == "sb")
    return TypeFiltre::COUPE_BANDE;
  msg_erreur("Type filtre inconnu : {}", s);
  return TypeFiltre::PASSE_BAS;
}


FRat<cfloat> design_riia(int n, const string &type,
    const string &prototype, float fcut, float δ_bp, float δ_bc)
{
  return design_riia(n, parse_tf(type), parse_proto(prototype), fcut, δ_bp, δ_bc);
}

FRat<float> design_biquad(const std::string type, float f, float Q, float gain_dB)
{
  BiquadSpec spec{.f = f, .Q = Q, .gain_dB = gain_dB};
  if((type == "lp") || (type == "pb"))
    spec.type = BiquadSpec::PASSE_BAS;
  else if((type == "hp") || (type == "ph"))
    spec.type = BiquadSpec::PASSE_HAUT;
  else if((type == "bp") || (type == "passe-bande"))
    spec.type = BiquadSpec::PASSE_BANDE;
  else if((type == "cb") || (type == "notch") || (type == "sb"))
    spec.type = BiquadSpec::COUPE_BANDE;
  else if((type == "plateau-bf"))
    spec.type = BiquadSpec::PLATEAU_BF;
  else if((type == "plateau-hf"))
    spec.type = BiquadSpec::PLATEAU_HF;
  else if((type == "res"))
    spec.type = BiquadSpec::RESONATEUR;
  return design_biquad(spec);
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
  return ss;
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

  auto z = (FRat<float>::z()).inv();

  return (b0 + b1 * z + b2 * z.pow(2)) / (1 + a1 * z + a2 * z.pow(2));
}




}

