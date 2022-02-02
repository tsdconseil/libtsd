#pragma once

// TODO : FICHIER OBSOLETE

#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"

namespace tsd::telecom {



class SuiviPicFrequence
{
public:
  int configure(unsigned int N);
  void step(const ArrayXcf &x, float &freq_detectee, float &snr);
private:
  sptr<tsd::fourier::FFTPlan> fft_plan;
  unsigned int N;
};


extern std::tuple<float,float> localise_pic_frequence(const ArrayXcf &x);


/** Détecteur d'erreur de phase quadratique */
/*class DEPQuadratique
{
public:
  // Par défaut, BPSK standard
  // Mettre 4.0 pour de la QPSK, etc.
  int configure(float power = 2.0);
  void step(cfloat x, float &dep);
  void step(ArrayXcf &x, float &df);

private:
  float power = 2.0;
  float dernier_dep = 0.0;
  float omega = 0.0;
};*/


struct PLL1Config
{
  /** Time constant */
  float tc = 10;
};

struct PLL2Config
{
  /** Loop normalized bandwidth */
  float BL;
  /** Damping factor */
  float η;
};

//extern sptr<FiltreConfigurable<float, float, PLL1Config>> pll_ordre_1();
//extern sptr<FiltreConfigurable<float, float, PLL2Config>> pll_ordre_2();



// Il faut supprimer cette classe !!!
class PLLOrdre2
{
public:

  PLLOrdre2()
  {
    configure({1,1});
  }

  void configure(const PLL2Config &config)
  {
    this->config = config;
    // PED gain at origin = 1 (ped supposed to be normalized)
    float A = 1;
    float η2 = config.η * config.η;

    γ   = (16 * η2 * config.BL) / (A * (1 + 4 * η2)) ;
    rho = (4 * config.BL) / (1 + 4 * η2);
    
    // Si A = 1, gamma = 4 * eta² * rho

    theta    = 0;
    mu       = 0;
    last_ped = 0;

    msg("Carrier rec. loop: BW = {} Hz, eta = {}", config.BL, config.η);
    msg("                   gamma = {}, rho = {}.", γ, rho);
  }

  /** @param ped Phase error detector, in radians */
  float process(float ped)
  {
    theta += mu;
    mu    += γ * ((1 + rho) * ped - last_ped);
    last_ped = ped;
    return theta;
  }

  /** @param freq Initial frequency, in Hz.
   *  @param phase Initial phase, in radians. */
  void set_initial_state(float freq, float phase)
  {
    mu    = freq * 2 * π;
    theta = phase;
  }

  float get_current_freq() const
  {
    return mu / (2 * π);
  }

  float get_current_phase() const
  {
    return theta;
  }

private:

  /** Coefficients for carrier recovery loop filter */
  float γ, rho;
  float theta, mu, last_ped;
  PLL2Config config;
};

}

