// implémentation filtre CIC.

#include "tsd/filtrage.hpp"

using namespace std;
using namespace tsd::vue;


namespace tsd::filtrage {


template<typename T, typename Ti>
struct FiltreCIC: FiltreGen<T>
{
  CICConfig config;
  char mode; // 'd' ou 'i'
  Vecteur<Ti> mem_diff;
  float gain = 1;

  FiltreCIC(const CICConfig &config, char mode = 'd')
  {
    this->config  = config;
    this->mode    = mode;

    si((mode != 'd') && (mode != 'u') && (mode != 'i'))
      échec("cic_init: le mode doit être 'd' ou 'i'.");
    si(config.M != 1)
      échec("cic_init: seulement M = 1 est supporté pour l'instant.");
    mem_diff.setZero(config.N);
    soit RM = config.R*config.M, N = config.N;
    si(mode == 'd')
      gain = 1.0f / pow(RM,N);
    sinon
      gain = ((float) config.R) / pow(RM, N);
  }

  void step(const Vecteur<T> &x, Vecteur<T> &y)
  {
    // Conversion évenutuelle en virgule fixe (entiers)
    soit xi = x.template as<Ti>();
    soit N = config.N, R = config.R;

    ////////////////////////////////////
    /// DECIMATION  ////////////////////
    ////////////////////////////////////
    si(mode == 'd')
    {
      // (1) Integration
      pour(auto i = 0; i < N; i++)
        xi = cumsum(xi);

      // (2) Decimation
      soit xd = sousech(xi, R);

      // (3) Peignes
      pour(auto i = 0; i < N; i++)
      {
        soit n = xd.rows();
        Vecteur<Ti> xp(n + 1);
        xp(0) = mem_diff(i);
        xp.segment(1, n) = xd;
        mem_diff(i) = xd(n-1);
        xd = diff(xp);
      }
      y = xd.template as<float>() * gain;
      retourne;
    }

    ////////////////////////////////////
    /// INTERPOLATION  /////////////////
    ////////////////////////////////////
    // (1) Peignes
    pour(auto i = 0; i < N; i++)
    {
      soit n = xi.rows();
      Vecteur<Ti> xp(n + 1);
      xp(0) = mem_diff(i);
      xp.segment(1,n) = xi;
      mem_diff(i) = x(n-1);
      xi = diff(xp);
    }


    // (2) Upsampling (insertion de zéros)
    soit n = xi.rows();
    soit tmp = Vecteur<Ti>::zeros(n * R);

    pour(auto i = 0; i < n; i++)
      tmp(i * R) = xi(i);

    xi = tmp;

    // (3) Integration
    pour(auto i = 0; i < N; i++)
      xi = cumsum(xi);

    y = xi.template as<float>() * gain;
  }
};


template<typename T, typename Ti>
sptr<FiltreGen<T>> filtre_cic(const CICConfig &config, char mode)
{
  retourne make_shared<FiltreCIC<T,Ti>>(config, mode);
}

FRat<float> design_cic(const CICConfig &config)
{
  assertion_msg((config.R > 0) && (config.N >= 0) && (config.M > 0),
      "Design CIC: paramètres invalides (R={}, N={}, M={})", config.R, config.N, config.M);

  soit RM = config.R * config.M;
  FRat<float> h{Vecf::ones(RM) / RM, Vecf::ones(1)};
  retourne h.eval_inv_z().pow(config.N);
}

CICAnalyse cic_analyse(const CICConfig &config, float fe, float f1)
{
  CICAnalyse res;
  res.config = config;

  soit R = config.R;

  // Fréquence de sortie
  soit fs = fe / R;

  // Par défaut, mesure l'atténuation au niveau du Nyquist de sortie
  si(f1 == 0)
    f1 = fs / 2;

  soit h = design_cic(config);

  msg("h = {}", h);

  soit [fr,mag] = frmag(h, 4096);
  fr *= fe;
  soit lmag = mag2db(mag + 1e-30);

  // Trace entre 0 et f1
  soit idf1 = trouve_premier(fr > f1),
       idfs = trouve_premier(fr > fs/2);
  assertion((idf1 >= 0) && (idfs >= 0));

  msg("CIC analyse : R={}, fe={} Hz, fs={} Hz, f1={} Hz, idfs={}, idf1={}",
      R, fe, fs, f1, idfs, idf1);


  soit atten  = lmag(idfs),
       att    = lmag(idf1),
       attmax = lmag.head(idf1).valeur_min();
  msg("Attenuation pour f > fs/2 : {:2f} dB.", atten);
  msg("Attenuation at {:.2f} Hz: {:.2f} dB.", f1, att);
  msg("Attenuation max. entre 0 et {:.2f} Hz: {:.2f} dB.", f1, attmax);
  msg("E.g. en échelle linéaire : * {:.3f}", db2mag(attmax));

  // Calcul repliement
  soit f      = fr.head(idfs);
  soit nrep = clamp((entier) floor((mag.rows() - idfs) / idfs), 0, 4);

  soit m_alias = Tabf::zeros(idfs,nrep);
  pour(auto i = 0; i < nrep; i++)
  {
    m_alias.col(i) = lmag.segment(idfs*(i+1), idfs);
    si((i % 2) == 0)
      m_alias.col(i).reverseInPlace();
  }

  {
    soit &g = res.fig_repliement;
    g.clear();
    soit c = g.plot(f, lmag.head(idfs), "g-", "Baseband spectrum");
    c.def_couleur({0,120,0});
    pour(auto i = 0; i < nrep; i++)
    {
      soit c = g.plot(f, m_alias.col(i), "r-", "Aliasing #{}", i);
      soit v = (200.0f*i) / nrep;
      c.def_couleur({200,v,v});
    }

    g.def_pos_legende("se");
    g.titres("Output of CIC filter, including decimation",
             "Fréquence",
             "Atténuation (dB)");
  }

  {
    soit &g = res.fig_spectre_global;
    g.clear();
    soit c = g.plot(fr,lmag);
    c.def_remplissage(oui, non);//oui, -150);
    g.plot(fs/2, atten, "sr");
    //g.def_rdi({-0.05, -150, 0.6, 150});
    g.titres("CIC filter / global view", "Fréquence", "Atténuation (dB)");
  }

  {
    soit &g = res.fig_spectre_bf;
    g.clear();
    soit c = g.plot(fr.head(idf1),lmag.head(idf1));
    c.def_remplissage(oui, non);//oui, -150);
    g.plot(fs/2, atten, "sr");
    //g.def_rdi({-0.05/R, -150, 3*(0.05+0.5)/R, 150});
    g.titres("CIC filter / view centered on the passband", "Fréquence", "Atténuation (dB)");
  }

  res.h     = h;
  res.fr    = fr;
  res.mag   = mag;
  res.nbits = (entier) ceil(config.N * log2(R) - 1);
  msg("Nombre de bits additionnels nécessaire pour l'implémentation : {}.", res.nbits);
  retourne res;
}


// TODO: supprimer et remplacer directement par repfreq / frmag
Vecf cic_freq(const CICConfig &config, const Vecf &f)
{
  soit n   = f.rows(), RM  = config.R * config.M;
  soit mag = Vecf::ones(n);

  pour(auto i = 0; i < n; i++)
  {
    soit d = RM * sin(π*f(i));
    si(d != 0)
      mag(i) = pow((float) abs(sin(RM*π*f(i)) / d), (float) config.N);
  }
  retourne mag;
}


CICComp design_cic_comp(const CICConfig &config, float fe, entier R2, float fc, entier ncoefs)
{
  // Pour l'instant, R2 non utilisé !
  CICComp res;
  soit analyse = cic_analyse(config, fe, fc);
  soit R = config.R;

  // Fréquence d'échantillonnage en sortie de filtre CIC
  soit fs = fe / R;
  msg("CIC comp : R = {}, fe = {:.2f} Hz, fs (CIC) = {:.2f} Hz, fc={:.2f} Hz.", R, fe, fs, fc);

  // Locate frequency response after decimation
  // (after decimation, representable frequencies are 0 to Fout/2)
  soit idx = trouve_premier(analyse.fr > fs / 2);
  soit fr2 = analyse.fr.segment(0,idx-1),
       fm2 = analyse.mag.segment(0,idx-1);


  //msg("CIC COM: idx = {}, fm = {}", idx, analyse.mag);

  // Ideal response of the FIR compensation filter = 1 / CIC response
  soit idealc = 1.0f / fm2;

  // Et, aussi, coupe @ fc Hz
  soit id = 1 + trouve_premier(analyse.fr > fc);

  assertion((id > 1) && (id + 1 < analyse.fr.rows()));

  // Reduce the steepness of required frequency profile
  idealc(id-1)  = idealc(id-2) / 2;
  idealc(id)    = idealc(id-1) / 2;
  idealc(id+1)  = idealc(id)   / 2;
  idealc.tail(idealc.rows() - (id+2)).setZero();


  // Make a fir filter from frequency sampling technique
  //msg("CIC COM: idealc = {}", idealc);
  res.h = design_rif_freq(ncoefs, idealc);


  //res.fig_rif_freq = design_rif_freq_analyse(ncoefs, idealc);

  // Compute frequency response of the FIR compensation filter
  soit hc = FRat<float>::rif(res.h);
  soit [frc,fmc] = frmag(res.h);

  // Fonction de transfert du filtre de compensation
  // decim(R) -> hc = hc(z^R)
  soit hc_lb = hc.horner(FRat<float>::z_power(config.R));

  // Fonction de transfert globale :
  // h(z) -> decim(R) -> hc(z) = h(z) * hc(z^R)
  soit hg = analyse.h * hc_lb;

  // Réponse filtre global
  soit [frg,fmg]      = frmag(hg);
  // Réponse filtre de compensation
  soit [frg2,fmhc_lb] = frmag(hc_lb);

  /////////////////////////////////////////////////////
  // Tracé des différentes réponses
  /////////////////////////////////////////////////////
  soit f = res.fig_reponse_ideale;
  f.clear();
  f.plot(fr2, idealc);

  f = res.fig_spectre_global;
  f.clear();
  f.plot(analyse.fr, mag2db(analyse.mag), "g-", "CIC filter");
  f.plot(frg2*fe, mag2db(fmhc_lb), "m-", "Compensation filter");
  f.plot(frg*fe, mag2db(fmg), "b-", "Global filter");
  f.titres("Decimation / wideband view", "Fréquence", "Attenuation (dB)");

  f = res.fig_spectre_bf;
  f.clear();
  f.plot(fr2,mag2db(fm2), "g-", "CIC filter");
  f.plot(frc*fs,mag2db(fmc), "m-","Compensation filter");
  soit fmaxz = 2 * fc;
  idx = trouve_premier(frg*fe > fmaxz);
  f.plot(frg.head(idx)*fe,mag2db(fmg.head(idx)),"b-","Global filter");
  f.titres("Decimation / passband view", "Fréquency", "Atténuation (dB)");

  res.fig_comp_rimp = plot_filtre(res.h);

  retourne res;
}

namespace hidden {
soit filtre_cic1 = filtre_cic<float, entier>;
}


}



