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
 * D'après : Multirate signal processing pour communication systems, F.J. Harris, page 89.
 */
static float cs(float t, float β)
{
  soit δ = abs(t) - 1 / (2*β);
  si((β > 0) && (abs(δ) < 1e-7f))
    retourne sinc(t) * π / 4;
  retourne sinc(t) * cos(π*β*t) / (1 - carré(2 * β * t));
}

/** Racine de cosinus sur-élevé, d'après :
 *  Multirate signal processing pour communication systems, F.J. Harris, page 90. */
static float rcs(float t, float β)
{
  soit δ = abs(t) - 1 / (4*β);
  si(t == 0)
    retourne (1-β) + 4*β/π;
  sinon si ((β > 0) && (abs(δ) < 1e-7))
    retourne (β/sqrt(2)) * ((1+(2/π))*sin(π/(4*β))+(1-(2/π))*cos(π/(4*β)));
  retourne (sin(π*t*(1-β)) + 4 * β * t * cos(π * t * (1+β)))
          / (π * t * (1 - carré(4 * β * t)));
}



Vecf design_rif_gaussien_gen(entier n, float σ, float dt = 0)
{
  soit t = linspace(-(n-1)/2.0f, (n-1)/2.0f, n) + dt,
       h = exp(-carré(t)/(2*σ*σ));
  si(h.somme() > 0)
    h /= h.somme();
  retourne h;
}


Vecf design_rif_gaussien(entier n, float σ)
{
  retourne design_rif_gaussien_gen(n, σ, 0);
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
  soit n = t.rows();
  ArrayXf y = ArrayXf::Zero(n);

  //float α = BT * sqrt(2*π/log(2.0f));

  pour(auto i = 0; i < n; i++)
    y(i) = exp(-t(i)*t(i)/(2 * sigma * sigma));//α * exp(-(2.0f/log(2.0f)) * pow(BT*π*t(i), 2.0));
  retourne y / y.sum();
}*/

float design_rif_gaussien_telecom_BT_vers_sigma(float BT)
{
  retourne sqrt(log(2.0f)) / (2*BT*π);
}



Vecf design_rif_gaussien_telecom(entier n, float BT, entier osf)
{
  // Prédécalage de la réponse gaussienne, de manière à avoir un filtre final centré
  soit dec_nrz = (osf-1.0f)/2;
  soit σ = design_rif_gaussien_telecom_BT_vers_sigma(BT);

  msg("rif_Gaussien_telecom: σ={}", σ);
  soit h1 = design_rif_gaussien_gen(n, σ * osf, dec_nrz);

  // BT trop important et OSF trop faible --> le filtre Gaussien est un dirac,
  // et après décalage temporel, tous les coefficients sont nuls.
  si(abs(h1).valeur_max() == 0)
  {
    msg_avert("design_rif_gaussien_telecom: BT={}, osf={} : résolution insuffisante pour sampler correctement la forme d'onde Gaussienne", BT, osf);
    h1(n/2) = 1;
    // Le délais introduit une distortion car le signal Gaussien ne respecte pas le critère de Nyquist.
    h1 = tsd::fourier::délais(h1, -dec_nrz);
  }

  soit h0 = Vecf::ones(osf) / osf,
       h2 = convol<float, float>(h0, h1);

  // Corrige le retard du au filtre
  //h2 = h2.tail(h2.rows() - (osf-1)/2).eval();
  //h2 = h2.head(n).eval();


  retourne h2 / h2.somme();
}


// D'après "Multirate signal processing pour communication systems", F.J. Harris, page 90.
Vecf design_rif_rcs(entier n, float β, float fcoupure)
{
  // fsymb = fcoupure*2
  // osf = 1/fsymb
  // => osf = 1/(2*fcut)
  retourne design_rif_rcs1(n, β, 1/(2*fcoupure));
}

Vecf design_rif_rcs1(entier n, float β, float osf, char nrm)
{
  assertion(n > 0);
  // n odd
  // ex n = 3 --> -1 0 1
  // n even
  // ex n = 4 --> -1.5, -0.5, 0.5, 1.5 ????
  //t = (-(n-1)/2:1:(n-1)/2) ./ osf;
  //ArrayXf t = linspace(-((float) n-1)/2, ((float) n-1)/2, n) / osf;
  //ArrayXf h = t.unaryExpr([&](float t){retourne rcs(t, β);});

  soit m = n / 2;

  soit h = Vecf::int_expr(n, IMAP(
      rcs(((float) i - m) / osf, β)));


  si(nrm == 'e')
    h /= sqrt(square(h).somme());
  sinon
    h /= h.somme();

  retourne h;
}

Vecf design_rif_cs(entier n, float β, float fcoupure)
{
  soit m = n / 2;
  Vecf h(n);
  pour(auto i = 0; i < n; i++)
  {
    // ncoefs = 7 --> m = 3

    // pour t entre -m et m :
    //  - 0   -> -m
    //  - n-1 ->  m
    // n = 3 -> -1, 0, 1
    soit t = i - m;
    // n = 4 -> -2, 1, 0, 1
    // La fonction cs s'annule pour t entier non nul, si fcut = 0.5 => Dirac (filtre neutre)
    h(i) = cs(t * fcoupure * 2, β);
  }

  retourne h / h.somme();
}

}

