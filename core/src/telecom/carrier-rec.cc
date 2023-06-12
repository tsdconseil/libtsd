#include "tsd/telecom/carrier-rec.hpp"
#include "tsd/telecom.hpp"
#include "tsd/vue.hpp"
#include <cmath>
#include <cassert>

using namespace tsd::vue;

const auto CREC_MODE_SAFE = non;

namespace tsd::telecom {



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
    retourne θ;
  }
};

struct LF1: FiltreBoucle
{
  float α, θ;

  LF1(float τ)
  {
    α = tsd::filtrage::lexp_tc_vers_coef(τ);
    reset();
  }
  float step(float x)
  {
    θ += α * x;
    retourne θ;
  }
  void reset()
  {
    θ = 0;
  }
};



sptr<FiltreBoucle> filtre_boucle_ordre_1(float τ)
{
  retourne std::make_shared<LF1>(τ);
}

sptr<FiltreBoucle> filtre_boucle_ordre_2(float BL, float η)
{
  retourne std::make_shared<LF2>(BL, η);
}

Ped ped_costa(entier M)
{
  si((M != 2) && (M != 4))
  {
    msg_erreur("ped_costa: fonctionne seulement en BPSK (M=2) ou QPSK (M=4), ici M = {}.", M);
    retourne {};
  }
  // QPSK : m^4 = 1
  si(M == 2)
    // BPSK : m² = 1
    // sin(phi) * cos(phi) = 0.5 sin(2phi) ~ phi
    // Note: c'est la même chose qu'une squaring loop
    // en prenant slmt la partie imaginaire.
    retourne [](cfloat x)
    {
      retourne x.real() * x.imag();
    };

  assertion(M == 4);

  // QPSK costa loop locks to square constellation
  // And we expect a "losange" constellation.
  retourne [](cfloat x){
    soit z = x * std::polar(1.0f, π_f/4);
    retourne z.imag() * signe(z.real()) - z.real() * signe(z.imag());
  };
}
Ped ped_ploop(entier M)
{
  retourne [M](cfloat x)
  {
    // Attention, nécessite une CAG en amont
    retourne (pow(x, M)).imag() / M;
  };
}
Ped ped_tloop(entier M)
{
  retourne [M](cfloat x) -> float
  {
    si(x == 0.0f)
      retourne 0;
    retourne std::arg(pow(x, M)) / M;
  };
}
Ped ped_decision(sptr<FormeOnde> wf)
{
  retourne [wf](cfloat x)
  {
    soit c = x * conj(wf->lis_symbole(wf->symbole_plus_proche(x)));
    si(c == 0.0f)
      retourne 0.0f;
    retourne std::arg(c);
  };
}

Ped ped_init(PedType type, sptr<FormeOnde> wf)
{
  si(type == PedType::AUTO)
  {
    si(wf->infos.est_psk)
      type = PedType::POWER_LOOP;
    sinon si(wf->infos.est_ask)
      type = PedType::TAN_LOOP;
    sinon
      type = PedType::DEC_LOOP;
  }

  soit M = wf->infos.M;
  si(wf->infos.est_ask)
    M = 2;

  si(type == PedType::COSTA)
    retourne ped_costa(M);
  sinon si(type == PedType::POWER_LOOP)
    retourne ped_ploop(M);
  sinon si(type == PedType::TAN_LOOP)
    retourne ped_tloop(M);
  sinon si(type == PedType::DEC_LOOP)
    retourne ped_decision(wf);
  échec("PED : type inconnu ({}).", (entier) type);
  retourne {};
}


tuple<float,float> localise_pic_frequence(const Veccf &x)
{
  SuiviPicFrequence suivi;
  suivi.configure(x.rows());
  float f, snr;
  suivi.step(x, f, snr);
  retourne {f, snr};
}

void SuiviPicFrequence::configure(entier N)
{
  this->N = N;
  fft_plan = tsd::fourier::tfrplan_création(N);
}

void SuiviPicFrequence::step(const Veccf &x, float &freq_detectee, float &snr)
{
  soit X = fft_plan->step(x);
  soit a2 = abs2(X);
  soit [y2, i2] = a2.max();

  snr = y2 / a2.moyenne();

  // Interpolation barycentrique (https://dspguru.com/dsp/howtos/how-to-interpolate-fft-peak/)
  soit i1 = (i2 - 1 + N) % N;
  soit i3 = (i2 + 1) % N;

  soit y1 = a2(i1);
  soit y3 = a2(i3);

  //    i = N/2 <=> f = fs/2
  // => i       <=> i * fs / N
  si(i2 >= (entier) N/2)
    i2 = i2 - N; // 0.6 => -0.4

  freq_detectee = ((float) i2) / N;
  assertion(abs(freq_detectee) <= 0.5);

  soit d = (y3 - y1) / (y1 + y2 + y3 + 1e-30f);

  assertion(abs(d) <= 1);

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
  void configure_impl(const RPLLConfig &c)
  {
    cpll = cpll_création(c.pll_interne);

    //cpll->configure(c.pll_interne);

    msg("Oscillateur : freq = {}", -c.freq);

    osc = source_ohc(-c.freq);
    //ol  = source_ohc(0);
    //sinon
    //ol = source_ohc(c.freq);


    soit coefs = tsd::filtrage::design_rif_cs(c.ncoefs_bb, 0.1, c.bp / 2);
    rif = tsd::filtrage::filtre_rif<float, cfloat>(coefs);
  }
  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit &config = Configurable<RPLLConfig>::config;
    soit n = x.rows();

    Veccf sosc, x1, x2;

    // Transposition bande de base
    sosc = osc->step(n);
    x1 = sosc * x;

    // Filtrage passe-bas
    si(config.filtre_bb_actif)
      x2 = rif->step(x1);
    sinon
      x2 = x1;

    soit x3 = cpll->step(x2);

    // Regénération de la porteuse
    //si(config.sortie_porteuse)
    {
      // Regénére le signal à la fréquence porteuse
      soit y1 = sosc.conjugate().head(n) * x3;
      // si signal complexe attendu en sortie
      //si constexpr (std::is_same<T, cfloat>::value)
        //y = y1;
      // si signal réel attendu en sortie
      //sinon
        y = real(y1);
    }
    /*sinon
    {
      si constexpr (std::is_same<T, cfloat>::value)
        y = mem_cor;
      sinon
      {
        y = mem_cor.real();
        //erreur("PLL : freq attendue nulle : le signal résultant est forcément complexe.");
      }
    }*/




    si(config.debug)
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
      soit c = f2.plot(real(y), "g-", "PLL");
      c.def_epaisseur(2);
      f2.afficher();
#     endif
    }
  }
};




// Création d'un signal périodique accroché sur un autre,
// tel que passé en entrée
template<typename T>
struct CPLL: Filtre<T, T, PLLConfig>, AFigures
{
  sptr<FiltreBoucle> lf;
  OLUT lut;

  float θ = 0;
  Ped ped;

  CPLL(const PLLConfig &cfg)
  {
    Configurable<PLLConfig>::configure(cfg);
  }

  // Suivi à l'ordre :
  // - 1 : erreur de phase uniquement
  // - 2 : erreur de fréquence
  void configure_impl(const PLLConfig &c)
  {
    soit &config = Configurable<PLLConfig>::config;
    config = c;

    msg("Configuration cpll : fréq oscillateur = {}.", -c.freq);

    si(!config.ped)
    {
      msg("cpll : ped par défaut.");
      ped = [](cfloat x){retourne std::arg(x);};
    }
    sinon
      ped = config.ped;

    si(config.loop_filter_order == 1)
      lf = filtre_boucle_ordre_1(config.tc);
    sinon
      lf = filtre_boucle_ordre_2(c.bp/* BL */, 1.0f/* η */);
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    soit &config = Configurable<PLLConfig>::config;
    soit n = x.rows();
    y.resize(n);

    Vecf pe, vθ;
    si(config.debug)
    {
      pe.resize(n);
      vθ.resize(n);
    }

    assertion(ped);

    pour(auto i = 0; i < n; i++)
    {
      y(i) = x(i) * lut.step(-θ);

      soit err = ped(y(i));

      si(config.debug)
      {
        vθ(i) = θ;
        pe(i) = err;
      }
      θ = lf->step(err);
    }

    si(config.debug)
    {
      soit &figs = AFigures::figures;
      figs.clear();
      figs.subplot().plot_iq(x, "Ab", "Entrée");
      figs.gcf().def_rdi({-1.1f,-1.1f,2.2f,2.2f});
      figs.subplot().plot_iq(y, "Ag", "Sortie");
      figs.gcf().def_rdi({-1.1f,-1.1f,2.2f,2.2f});
      figs.subplot().plot(pe * 180 / π_f, "-b", "Erreur de phase (deg.)");
      figs.subplot().plot(vθ * 180 / π_f, "-b", "Phase (deg.)");
      figs.afficher("Recouvrement de porteuse");
    }

    si constexpr(CREC_MODE_SAFE)
    {
      si(y.hasNaN())
        msg_erreur("Carrier rec : NaN.");
    }
  }
};



sptr<Filtre<float, float, RPLLConfig>> rpll_création(const RPLLConfig &cfg)
{
  retourne std::make_shared<RPLL<float>>(cfg);
}

sptr<Filtre<cfloat, cfloat, PLLConfig>> cpll_création(const PLLConfig &cfg)
{
  retourne std::make_shared<CPLL<cfloat>>(cfg);
}


}

