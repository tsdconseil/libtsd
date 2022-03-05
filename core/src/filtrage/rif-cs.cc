#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include "tsd/divers.hpp"

using namespace std;


namespace tsd::filtrage {



/** Cosinus surélevé
 *
 * @param t: temps, en multiple de la période symbole
 * @param β: facteur de dépassement
 *
 * Calcule <latex>$$y(t)=\sinc(\pi t) \frac{\cos(\pi rt)}{1-4r^2t^2}$$</latex>
 *
 * D'après : Multirate signal processing for communication systems, F.J. Harris, page 89.
 */
static float cs(float t, float β)
{
  float δ = abs(t) - 1.0f / (2*β);
  if((β > 0) && (abs(δ) < 1e-7f))
    return sinc(t) * π / 4;
  else
    return sinc(t) * cos(π*β*t) / (1 - carré(2 * β * t));
}

/** Racine de cosinus sur-élevé, d'après :
 *  Multirate signal processing for communication systems, F.J. Harris, page 90. */
static float rcs(float t, float β)
{
  auto δ = abs(t) - 1.0f / (4*β);
  if(t == 0)
    return (1-β) + 4*β/π;
  else if ((β > 0) && (abs(δ) < 1e-7))
    return (β/sqrt(2)) * ((1+(2/π))*sin(π/(4*β))+(1-(2/π))*cos(π/(4*β)));
  return (sin(π*t*(1-β)) + 4 * β * t * cos(π * t * (1+β)))
          / (π * t * (1 - carré(4 * β * t)));
}



ArrayXf design_rif_gaussien_gen(int n, float σ, float dt = 0)
{
  ArrayXf t = linspace(-((float) n-1)/2, ((float) n-1)/2, n) + dt;
  ArrayXf h = (-t.square()/(2*σ*σ)).exp();
  if(h.sum() > 0)
    h /= h.sum();
  return h;
}


ArrayXf design_rif_gaussien(int n, float σ)
{
  return design_rif_gaussien_gen(n, σ, 0);
}

/*static ArrayXf gaussien_impulsion(const ArrayXf &t, float BT)
// Gaussian pulse
//
// Calling Sequence
//   y = gaussian_pulse(t,BT)
//
// Parameters
// t: time, in multiple of symbol period (1 = 1 symbol period)
// BT: Bandwidth-time product
//
// Description
//  Compute a gaussian filter with given Bandwidth Time product.
{
  float sigma = design_rif_gaussien_telecom_BT_vers_sigma(BT);
  auto n = t.rows();
  ArrayXf y = ArrayXf::Zero(n);

  //float α = BT * sqrt(2*π/log(2.0f));

  for(auto i = 0; i < n; i++)
    y(i) = exp(-t(i)*t(i)/(2 * sigma * sigma));//α * exp(-(2.0f/log(2.0f)) * pow(BT*π*t(i), 2.0));
  return y / y.sum();
}*/

float design_rif_gaussien_telecom_BT_vers_sigma(float BT)
{
  return sqrt(log(2.0f)) / (2*BT*π);
}



ArrayXf design_rif_gaussien_telecom(int n, float BT, int osf)
{
  // Prédécalage de la réponse gaussienne, de manière à avoir un filtre final centré
  float dec_nrz = (osf-1.0)/2;
  //ArrayXf t = linspace(-((float) n-1)/2+dec_nrz, ((float) n-1)/2+dec_nrz, n) / osf;

  float σ = design_rif_gaussien_telecom_BT_vers_sigma(BT);

  msg("rif_Gaussien_telecom: σ={}", σ);
  ArrayXf h1 = design_rif_gaussien_gen(n, σ * osf, dec_nrz);
      //gaussien_impulsion(t, BT);

  // BT trop important et OSF trop faible --> le filtre Gaussien est un dirac,
  // et après décalage temporel, tous les coefficients sont nuls.
  if(h1.abs().maxCoeff() == 0)
  {
    msg_avert("design_rif_gaussien_telecom: BT={}, osf={} : résolution insuffisante pour sampler correctement la forme d'onde Gaussienne", BT, osf);
    h1(n/2) = 1;
    // Le délais introduit une distortion car le signal Gaussien ne respecte pas le critère de Nyquist.
    h1 = tsd::fourier::délais(h1, -dec_nrz);
  }

  ArrayXf h0 = ArrayXf::Ones(osf) / osf;

  //h1 = h1 | ArrayXf::Zero(5 * osf);

  ArrayXf h2 = convol<float, float>(h0, h1);

  // Corrige le retard du au filtre
  //h2 = h2.tail(h2.rows() - (osf-1)/2).eval();
  //h2 = h2.head(n).eval();


  return h2 / h2.sum();
}


// D'après "Multirate signal processing for communication systems", F.J. Harris, page 90.
ArrayXf design_rif_rcs(int n, float β, float fcoupure)
{
  // fsymb = fcoupure*2
  // osf = 1/fsymb
  // => osf = 1/(2*fcut)
  return design_rif_rcs1(n, β, 1/(2*fcoupure));
}

ArrayXf design_rif_rcs1(int n, float β, float osf, char nrm)
{
  tsd_assert(n > 0);
  // n odd
  // ex n = 3 --> -1 0 1
  // n even
  // ex n = 4 --> -1.5, -0.5, 0.5, 1.5 ????
  //t = (-(n-1)/2:1:(n-1)/2) ./ osf;
  //ArrayXf t = linspace(-((float) n-1)/2, ((float) n-1)/2, n) / osf;
  //ArrayXf h = t.unaryExpr([&](float t){return rcs(t, β);});

  int m = n / 2;
  ArrayXf h(n);
  for(auto i = 0; i < n; i++)
  {
    float t = i - m;
    h(i) = rcs(t / osf, β);
  }


  if(nrm == 'e')
    h /= sqrt(h.square().sum());
  else
    h /= h.sum();

  return h;
}

ArrayXf design_rif_cs(int n, float β, float fcoupure)
{
  int m = n / 2;
  ArrayXf h(n);
  for(auto i = 0; i < n; i++)
  {
    // ncoefs = 7 --> m = 3

    // Pour t entre -m et m :
    //  - 0   -> -m
    //  - n-1 ->  m
    // n = 3 -> -1, 0, 1
    float t = i - m;
    // n = 4 -> -2, 1, 0, 1
    // La fonction cs s'annule pour t entier non nul, si fcut = 0.5 => Dirac (filtre neutre)
    h(i) = cs(t * fcoupure * 2, β);
  }

  return h / h.sum();
}

}

