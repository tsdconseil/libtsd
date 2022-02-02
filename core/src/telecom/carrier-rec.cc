#include "tsd/telecom/carrier-rec.hpp"
#include "tsd/telecom.hpp"
#include "tsd/figure.hpp"
#include <cmath>
#include <cassert>

using namespace tsd::vue;

const auto CREC_MODE_SAFE = false;

namespace tsd::telecom {

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


/** Calcul d'exponentielle à partir d'une LUT */
struct LUTOsc
{
  int n;
  ArrayXcf table;

  LUTOsc(int n)
  {
    this->n = n;
    table.resize(n);
    for(auto i = 0; i < n; i++)
      table(i) = std::polar(1.0, i * 2 * π / n);
  }
  cfloat step(float θ)
  {
    θ = modulo_2π(θ);
    int x0 = floor((n * θ) / (2 * π));
    return table(x0);
  }
};





/** Filtre de boucle du second ordre */
struct LF2: FiltreBoucle
{
  float A, BL, η, γ, ρ, θ, μ, last_ped;

  LF2(float BL, float η)
  {
    this->BL  = BL;
    this->η   = η;
    A     = 1; // PED gain at origin = 1 (ped supposed to be normalized)
    γ  = (16 * η * η * BL)  / (A * (1 + 4 * η * η));
    ρ   = (4 * BL) / (1 + 4 * η * η);
    reset();
  }
  void reset()
  {
    θ = μ = last_ped = 0;
  }
  float step(float x)
  {
    θ += μ;
    μ += γ * ((1 + ρ) * x - last_ped);
    last_ped = x;
    return θ;
  }
};

struct LF1: FiltreBoucle
{
  float α, θ;

  LF1(float τ)
  {
    α = tsd::filtrage::rii1_tc_vers_coef(τ);
    reset();
  }
  float step(float x)
  {
    θ += α * x;
    return θ;
  }
  void reset()
  {
    θ = 0;
  }
};



sptr<FiltreBoucle> filtre_boucle_ordre_1(float τ)
{
  return std::make_shared<LF1>(τ);
}

sptr<FiltreBoucle> filtre_boucle_ordre_2(float BL, float η)
{
  return std::make_shared<LF2>(BL, η);
}

Ped ped_costa(int M)
{
  if((M != 2) && (M != 4))
  {
    msg_erreur("ped_costa: fonctionne seulement en BPSK (M=2) ou QPSK (M=4), ici M = {}.", M);
    return Ped();
  }
  // QPSK : m^4 = 1
  if(M == 2)
    // BPSK : m² = 1
    // sin(phi) * cos(phi) = 0.5 sin(2phi) ~ phi
    // Note: c'est la même chose qu'une squaring loop
    // en prenant slmt la partie imaginaire.
    return [](cfloat x) {return x.real() * x.imag();};
  else if(M == 4)
  {
    // QPSK costa loop locks to square constellation
    // And we expect a "losange" constellation.
    return [](cfloat x){
    cfloat z = x * std::polar(1.0f, π_f/4);
    return z.imag() * sign(z.real()) - z.real() * sign(z.imag());
    };
  }
  return Ped();
}
Ped ped_ploop(int M)
{
  return [M](cfloat x){
  // Apparement, qd l'amplitude est légérement supérieure à 1, il y amplification exponentielle du gain ?
  //
  return (std::pow(x, M)).imag() / M;
  };
  //return (std::arg(std::pow(x, M))) / M;
}
Ped ped_tloop(int M)
{
  return [M](cfloat x) -> float {
    if(x == 0.0f)
      return 0;
    return std::arg(std::pow(x, M)) / M;
  };
}
Ped ped_decision(sptr<FormeOnde> wf)
{
  return [wf](cfloat x){
    auto c = x * conj(wf->lis_symbole(wf->symbole_plus_proche(x)));
    if(c == 0.0f)
      return 0.0f;
    return std::arg(c);
  };
}

Ped ped_init(PedType type, sptr<FormeOnde> wf)
{
  if(type == PedType::AUTO)
  {
    if(wf->infos.est_psk)
    {
      //if(wf.M >= 4) // QPSK, 8PSK, ...
        //scr.ted = ted_qpsk();
      type = PedType::POWER_LOOP;
    }
    else if(wf->infos.est_ask)
    {
      type = PedType::TAN_LOOP;
    }
    else
    {
      type = PedType::DEC_LOOP;
    }
    /*else
    {
      msg_erreur("ped_init() : ne sait pas quel PED créer pour cette modulation.");
      return Ped();
    }*/
  }

  int M = wf->infos.M;
  if(wf->infos.est_ask)
    M = 2;

  if(type == PedType::COSTA)
    return ped_costa(M);
  else if(type == PedType::POWER_LOOP)
    return ped_ploop(M);
  else if(type == PedType::TAN_LOOP)
    return ped_tloop(M);
  else if(type == PedType::DEC_LOOP)
    return ped_decision(wf);
  msg_erreur("PED : type inconnu ({} - {}).", (int) type, type);
  return Ped();
}

#if 0
struct PedCosta: Ped
{
  PedCosta(unsigned int M)
  {
    nom = "costa";
    this->M = M;
    require_agc = true;
    if((M != 2) && (M != 4))
      erreur("psk_costa: expect BPSK or QPSK.");
  }
  float calcule(const std::complex<float> &x)
  {
    // QPSK : m^4 = 1
    if(M == 2)
      // BPSK : m² = 1
      // sin(phi) * cos(phi) = 0.5 sin(2phi) ~ phi
      // Note: c'est la même chose qu'une squaring loop
      // en prenant slmt la partie imaginaire.
      return x.real() * x.imag();
    else if(M == 4)
    {
      // QPSK costa loop locks to square constellation
      // And we expect a "losange" constellation.
      cfloat z = x * std::polar(1.0f, pi/4);
      return z.imag() * sign(z.real()) - z.real() * sign(z.imag());
    }
    return 0;
  }
};

struct PedPLoop: Ped
{
  PedPLoop(unsigned int M)
  {
    nom = "power loop";
    this->M = M;
    require_agc = true;
  }
  float calcule(const std::complex<float> &x)
  {
    // Apparement, qd l'amplitude est légérement supérieure à 1, il y amplification exponentielle du gain ?
    //
    return (std::pow(x, M)).imag() / M;
    //return (std::arg(std::pow(x, M))) / M;
  }
};

struct PedDec: Ped
{
  sptr<FormeOnde> wf;
  PedDec(sptr<FormeOnde> wf)
  {
    nom = "decision based loop";
    this->wf = wf;
    require_agc = true;
  }

  float calcule(const std::complex<float> &x)
  {
    auto c = x * conj(wf->lis_symbole(wf->symbole_plus_proche(x)));
    if(c == 0.0f)
      return 0.0f;
    return std::arg(c);
  }
};


struct PedTLoop: Ped
{
  PedTLoop(unsigned int M)
  {
    nom = "tan loop";
    this->M = M;
    require_agc = false;
  }
  float calcule(const std::complex<float> &x)
  {
    //atan(imag(z .^ ped.M),real(z .^ ped.M)) / ped.M;
    if(x == 0.0f)
      return 0;
    return std::arg(std::pow(x, M)) / M;
  }
};


sptr<Ped> ped_basee_decision(sptr<FormeOnde> wf)
{
  return std::make_shared<PedDec>(wf);
}

sptr<Ped> ped_init(PedType type, unsigned int M)
{
  if(type == PedType::COSTA)
    return std::make_shared<PedCosta>(M);
  else if(type == PedType::POWER_LOOP)
    return std::make_shared<PedPLoop>(M);
  else if(type == PedType::TAN_LOOP)
    return std::make_shared<PedTLoop>(M);
  msg_erreur("PED : type inconnu ({} - {}).", (int) type, type);
  return sptr<Ped>();
}
#endif


#if 0
// A SUPPRIMER
struct CarrierRec: Filtre<cfloat, cfloat, CarrierRecConfig>
{
  float θ = 0;
  float rssi = 1.0f, g = 1.0f;


  CarrierRec(const CarrierRecConfig &config)
  {
    configure(config);
  }

  int configure(const CarrierRecConfig &config)
  {
    this->config = config;
    tsd_assert(config.ped);
    tsd_assert(config.lf);
    if(config.ped->require_agc)
    {
      g     = tsd::filtrage::iir1_tc_vers_coef(config.ped->agc_tc);
      rssi  = 1;
    }
    return 0;
  }

  void step(const Eigen::Ref<const ArrayXcf> x, ArrayXcf &y)
  {
    auto n = x.rows();
    y.resize(n);

    ArrayXf pe(n);
    ArrayXf vtheta(n);
    ArrayXcf ploop(n);
    ArrayXcf vagc(n);
    ArrayXf vrssi = ArrayXf::Ones(n);

    for(auto i = 0; i < n; i++)
    {
      vtheta(i) = θ;
      y(i) = x(i) * std::polar(1.0f, -θ);
      auto zagc = y(i);
      if(config.ped->require_agc)
      {
          rssi = (1-g) * rssi + g * std::abs(zagc);
          vrssi(i) = rssi;
          zagc = zagc / (rssi + 1e-10f);
      }
      vagc(i) = zagc;


      ArrayXf dphi(1);
      dphi(0) = config.ped->calcule(zagc);
      pe(i) = dphi(0);
      ploop(i) = config.ped->calcule(x(i));

      θ = config.lf->step(dphi)(0);
    }

    if(config.debug_actif)
    {
      Figure f("Recouvrement de porteuse");
      auto nl = 5, nc = 1, ids = 1;
      f.subplot(nl,nc,ids++);
      f.plot(vrssi, "-b", "RSSI");
      f.subplot(nl,nc,ids++);
      f.plot(pe * 180 / pi, "-b", "Erreur de phase (deg.)");
      f.subplot(nl,nc,ids++);
      f.plot(vtheta * 180 / pi, "-b", "Phase (deg.)");
      f.subplot(nl,nc,ids++);
      f.plot_iq(vagc, ".b", "AGC");
      f.subplot(nl,nc,ids++);
      f.plot_iq(ploop, ".b", "PLOOP");
      f.afficher();
    }

    if(y.hasNaN())
    {
      erreur("Carrier rec : NaN.");
    }

  }

};


sptr<Filtre<cfloat, cfloat, CarrierRecConfig>> carrier_rec_init(const CarrierRecConfig &config)
{
  return std::make_shared<CarrierRec>(config);
}
#endif

/*int DEPQuadratique::configure(float power)
{
  this->power = power;
  return 0;
}

void DEPQuadratique::step(cfloat x, float &dep)
{
  //dep = std::arg(x);
  dep = std::arg(std::pow(x, power)) / power;
  //dep = std::arg(x) * power;
}

void DEPQuadratique::step(ArrayXcf &x, float &df)
{
  auto n = x.rows();

  float alpha = 0.95;

  for(auto i = 0; i < n; i++)
  {
    float dep;
    step(x(i), dep);

    float d = dep - dernier_dep;

    if(d > pi)
      d = 2 * pi - d;
    if(d < -pi)
      d = d + 2 * pi;

    omega = alpha * omega + (1.0 - alpha) * d;

    dernier_dep = dep;
  }

  df = omega / (2 * pi);
}*/

std::tuple<float,float> localise_pic_frequence(const ArrayXcf &x)
{
  SuiviPicFrequence suivi;
  suivi.configure(x.rows());
  float f, snr;
  suivi.step(x, f, snr);
  return std::make_tuple(f, snr);
}

int SuiviPicFrequence::configure(unsigned int N)
{
  this->N = N;
  fft_plan = tsd::fourier::fftplan_création(N);
  return 0;
}

void SuiviPicFrequence::step(const ArrayXcf &x, float &freq_detectee, float &snr)
{
  ArrayXcf X;
  fft_plan->step(x, X, true);

  int i2 = 0;
  auto a2 = X.abs2();
  float y2 = a2.maxCoeff(&i2);

  snr = y2 / a2.mean();

  // Interpolation barycentrique (https://dspguru.com/dsp/howtos/how-to-interpolate-fft-peak/)
  int i1 = (i2 - 1 + N) % N;
  int i3 = (i2 + 1) % N;

  auto y1 = a2(i1);
  auto y3 = a2(i3);

  //    i = N/2 <=> f = fs/2
  // => i       <=> i * fs / N
  if(i2 >= (int) N/2)
    i2 = i2 - N; // 0.6 => -0.4

  freq_detectee = ((float) i2) / N;
  tsd_assert(std::abs(freq_detectee) <= 0.5);

  auto d = (y3 - y1) / (y1 + y2 + y3 + 1e-30f);

  tsd_assert(std::abs(d) <= 1);

  freq_detectee += d / N;
}

template<typename T>
struct RPLL: Filtre<T, T, RPLLConfig>
{
  sptr<Filtre<std::complex<T>, std::complex<T>, PLLConfig>> cpll;

  sptr<SourceGen<cfloat>> osc, ol;
  sptr<FiltreGen<float>> lf;
  sptr<FiltreGen<cfloat>> rif;

  float theta = 0;

  RPLL(const RPLLConfig &cfg)
  {
    Configurable<RPLLConfig>::configure(cfg);
  }

  // On peut faire un suivi à l'ordre :
  // - Ordre 1 : erreur de phase uniquement
  // - Ordre 2 : erreur de fréquence
  int configure_impl(const RPLLConfig &c)
  {
    cpll = cpll_création(c.pll_interne);

    //cpll->configure(c.pll_interne);

    msg("Oscillateur : freq = {}", -c.freq);

    osc = source_ohc(-c.freq);
    //ol  = source_ohc(0);
    //else
    //ol = source_ohc(c.freq);


    ArrayXf coefs = tsd::filtrage::design_rif_cs(c.ncoefs_bb, 0.1, c.bp / 2);
    rif = tsd::filtrage::filtre_rif<float, cfloat>(coefs);

    return 0;
  }
  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    auto &config = Configurable<RPLLConfig>::config;
    auto n = x.rows();

    ArrayXcf sosc, x1, x2;

    // Transposition bande de base
    sosc = osc->step(n);
    x1 = sosc * x;

    // Filtrage passe-bas
    if(config.filtre_bb_actif)
      x2 = rif->step(x1);
    else
      x2 = x1;

    ArrayXcf x3 = cpll->step(x2);

    // Regénération de la porteuse
    //if(config.sortie_porteuse)
    {
      // Regénére le signal à la fréquence porteuse
      ArrayXcf y1 = sosc.conjugate().head(n) * x3;
      // Si signal complexe attendu en sortie
      //if constexpr (std::is_same<T, cfloat>::value)
        //y = y1;
      // Si signal réel attendu en sortie
      //else
        y = y1.real();
    }
    /*else
    {
      if constexpr (std::is_same<T, cfloat>::value)
        y = mem_cor;
      else
      {
        y = mem_cor.real();
        //erreur("PLL : freq attendue nulle : le signal résultant est forcément complexe.");
      }
    }*/




    if(config.debug)
    {
#     if 0
      Figure f("PLL (1)");

      f.subplot(321);
      f.plot(real(x), "b-", "Signal entree");
      f.subplot(322);
      f.plot_psd(x);
      f.subplot(323);
      f.plot(x1.real(), "b-", "transposition (I)");
      f.plot(x1.imag(), "g-", "transposition (Q)");
      f.subplot(324);
      f.plot_psd(x1);
      f.subplot(325);
      f.plot(x2.real(), "b-", "passe-bas (I)");
      f.plot(x2.imag(), "g-", "passe-bas (Q)");
      f.subplot(326);
      f.plot_psd(x2);

      f.afficher();

      Figure f2("PLL (2)");
      f2.subplot(311);
      f2.plot(verr * 180 / pi, "r-", "Erreur de phase (degres)");
      f2.subplot(312);
      f2.plot(vtheta * 180 / pi, "b-", "Phase (degres)");



      //Figure fig("Démo PLL");

      f2.subplot(313);
      f2.plot(real(x), "b-", "Signal a suivre");
      auto c = f2.plot(real(y), "g-", "PLL");
      c.def_epaisseur(2);
      f2.afficher();
#     endif
    }
  }
};

// Création d'un signal périodique accroché sur un autre,
// tel que passé en entrée
template<typename T>
struct CPLL: Filtre<T, T, PLLConfig>
{
  LUTOsc lut_osc;

  //sptr<SourceGen<cfloat>> ol;
  sptr<FiltreBoucle> lf;

  float theta = 0;

  CPLL(const PLLConfig &cfg): lut_osc(256)
  {
    Configurable<PLLConfig>::configure(cfg);
  }

  // On peut faire un suivi à l'ordre :
  // - Ordre 1 : erreur de phase uniquement
  // - Ordre 2 : erreur de fréquence
  int configure_impl(const PLLConfig &c)
  {
    auto &config = Configurable<PLLConfig>::config;

    if(!config.ped)
      config.ped = [](cfloat x){return std::arg(x);};

    msg("Configuration cpll : fréq oscillateur = {}.", -c.freq);
    //osc = source_ohc(-c.freq);

    /*if(config.sortie_porteuse)
      ol  = source_ohc(0);
    else*/
    //ol = source_ohc(c.freq);

    if(config.loop_filter_order == 1)
      lf = filtre_boucle_ordre_1(config.tc);
    else
      lf = filtre_boucle_ordre_2(c.bp/* BL */, 1.0f/* η */);

    return 0;
  }
  void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
  {
    auto &config = Configurable<PLLConfig>::config;
    auto n = x.rows();
    y.resize(n);

    ArrayXf pe, vtheta;
    if(config.debug)
    {
      pe.resize(n);
      vtheta.resize(n);
    }

    for(auto i = 0; i < n; i++)
    {
      // TODO: utiliser un oscillateur LUT ici
      //y(i) = x(i) * lut_osc.step(-theta);//;
      y(i) = x(i) * std::polar(1.0f, -theta);

      auto err = config.ped(y(i));

      if(config.debug)
      {
        vtheta(i) = theta;
        pe(i) = err;
      }
      theta = lf->step(err);
    }

    if(config.debug)
    {
      Figures f;

      f.subplot().plot_iq(x, "ab", "Entrée");
      f.subplot().plot_iq(y, "ag", "Sortie");

      //f.subplot(nl,nc,ids++);
      //f.plot(vrssi, "-b", "RSSI");
      f.subplot().plot(pe * 180 / pi, "-b", "Erreur de phase (deg.)");
      f.subplot().plot(vtheta * 180 / pi, "-b", "Phase (deg.)");

      //f.subplot(nl,nc,ids++);
      //f.plot_iq(ploop, ".b", "PED");
      f.afficher("Recouvrement de porteuse");
    }

    if constexpr(CREC_MODE_SAFE)
    {
      if(y.hasNaN())
        msg_erreur("Carrier rec : NaN.");
    }

#   if 0
    auto n = x.rows();


    ArrayXf vtheta(n), verr(n);
    ArrayXcf mem_nco(n), mem_cor(n);

    // On a maintenant un signal, qui serait constant et égal à 1 si on était accroché.
    for(auto i = 0; i < n; i++)
    {
      // TODO : utiliser plutôt un oscillateur lut.
      //
      cfloat rot = ol->step(1)(0);
      //

      //cfloat rot = std::polar(1.0f, theta);
      cfloat cor = x2(i) * conj(rot);
      float  err = config.detecteur_erreur_phase(cor);

      //auto err = std::arg(e1);
      θ = lf->step(err);

      // Pour débug
      mem_nco(i)    = rot;
      mem_cor(i)    = cor;
      verr(i)       = err;
      vtheta(i)     = θ;
    }

    if(config.sortie_porteuse)
      y = ?;
    else
      y = mem_cor;




    if(config.debug)
    {
      Figure f("PLL (1)");

      f.subplot(321);
      f.plot(real(x), "b-", "Signal entree");
      f.subplot(322);
      f.plot_psd(x);
      f.subplot(323);
      f.plot(x1.real(), "b-", "transposition (I)");
      f.plot(x1.imag(), "g-", "transposition (Q)");
      f.subplot(324);
      f.plot_psd(x1);
      f.subplot(325);
      f.plot(x2.real(), "b-", "passe-bas (I)");
      f.plot(x2.imag(), "g-", "passe-bas (Q)");
      f.subplot(326);
      f.plot_psd(x2);

      f.afficher();

      Figure f2("PLL (2)");
      f2.subplot(311);
      f2.plot(verr * 180 / pi, "r-", "Erreur de phase (degres)");
      f2.subplot(312);
      f2.plot(vtheta * 180 / pi, "b-", "Phase (degres)");



      //Figure fig("Démo PLL");

      f2.subplot(313);
      f2.plot(real(x), "b-", "Signal a suivre");
      auto c = f2.plot(real(y), "g-", "PLL");
      c.def_epaisseur(2);
      f2.afficher();
    }
#   endif
  }
};



sptr<Filtre<float, float, RPLLConfig>> rpll_création(const RPLLConfig &cfg)
{
  return std::make_shared<RPLL<float>>(cfg);
}

sptr<Filtre<cfloat, cfloat, PLLConfig>> cpll_création(const PLLConfig &cfg)
{
  return std::make_shared<CPLL<cfloat>>(cfg);
}


}

