// implémentation filtre CIC.

#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/filtrage/frat.hpp"

using namespace std;
using namespace tsd::vue;


namespace tsd::filtrage {




  /** @brief CIC filtering class (for interpolation or decimation)
   *  @bref Filtre CIC (interpolation ou décimation) */
template<typename T, typename Ti>
struct FiltreCIC: FiltreGen<T>
{
  CICConfig config;
  char mode;
  Vecteur<Ti> mem_diff;
  float gain;

  FiltreCIC(const CICConfig &config, char mode = 'd')
  {
    this->config  = config;
    this->mode    = mode;

    if((mode != 'd') && (mode != 'u') && (mode != 'i'))
      echec("cic_init: le mode doit être 'd' ou 'i'.");
    if(config.M != 1)
      echec("cic_init: seulement M = 1 est supporté pour l'instant.");
    mem_diff.setZero(config.N);
    if(mode == 'd')
      gain = 1 / pow(config.R*config.M,config.N);
    else
      gain = config.R / pow(config.R*config.M,config.N);
  }

    void step(const Eigen::Ref<const Vecteur<T>> x, Vecteur<T> &y)
    {
      Vecteur<Ti> xi = x.template cast<Ti>();

      ////////////////////////////////////
      /// DECIMATION  ////////////////////
      ////////////////////////////////////
      if(mode == 'd')
      {
        // (1) Integration
        for(auto i = 0; i < config.N; i++)
          xi = cumsum(xi);

        // (2) Decimation
        Vecteur<Ti> xd = sousech(xi, config.R);

        // (3) Peignes
        for(auto i = 0; i < config.N; i++)
        {
          Vecteur<Ti> xp(xd.rows() + 1);
          xp(0) = mem_diff(i);
          xp.segment(1, xd.rows()) = xd;
          mem_diff(i) = xd(xd.rows()-1);
          xd = diff(xp);
        }
        y = xd.template cast<float>() * gain;
        return;
      }
      ////////////////////////////////////
      /// INTERPOLATION  /////////////////
      ////////////////////////////////////
      // (1) Peignes
      for(auto i = 0; i < config.N; i++)
      {
        Vecteur<Ti> xp(xi.rows() + 1);
        xp(0) = mem_diff(i);
        xp.segment(1,xi.rows()) = xi;
        mem_diff(i) = x(xi.rows()-1);
        xi = diff(xp);
      }

      // (2) Upsampling
      Vecteur<Ti> tmp = Vecteur<Ti>::Zero(xi.rows()*config.R);
      auto map = Eigen::Map<Vecteur<Ti>,0,Eigen::InnerStride<>>(tmp.data(), tmp.rows() / config.R, Eigen::InnerStride<>(config.R));
      map = xi;
      xi = tmp;

      // (3) Integration
      for(auto i = 0; i < config.N; i++)
        xi = cumsum(xi);

      y = xi.template cast<float>() * gain;
    }
  };


template<typename T, typename Ti>
sptr<FiltreGen<T>> filtre_cic(const CICConfig &config, char mode)
{
  return make_shared<FiltreCIC<T,Ti>>(config, mode);
}

FRat<float> design_cic(const CICConfig &config)
{
  tsd_assert((config.R > 0) && (config.M > 0) && (config.N >= 0));

  FRat<float> h;
  h.numer.coefs.setOnes(config.R * config.M);
  h.numer.coefs /= (config.R * config.M);
  return h.horner(FRat<float>::z().inv()).pow(config.N);
}

CICAnalyse cic_analyse(const CICConfig &config, float Fin, float Fint)
{
  CICAnalyse res;
  res.config = config;

  if(Fint == 0)
    Fint = Fin / config.R;

  auto h = design_cic(config);
  auto [fr,mag] = frmag(h, 4096);
  fr *= Fin;

  ArrayXf lmag = 20.0f*(mag+1e-10f).log10();

  // Trace entre 0 et Fint*2
  float fmaxz = Fint * 2;
  int idx = trouve(fr > fmaxz)[0];


  msg("CIC analyse : R={}, Fin={} Hz, fmaxz={} Hz, idx={}", config.R, Fin, fmaxz, idx);

  auto frz = fr.head(idx);
  auto fmz = mag.head(idx);

  auto Fout = Fin / config.R;
  msg("R = {}, Fin = {:.2f} Hz, Fout = {:.2f} Hz.", config.R, Fin, Fout);

  auto fc2 = Fout/2;
  idx = trouve(fr > fc2)[0];
  auto atten = 20*log10(mag(idx));
  msg("Attenuation for f > fout/2 : {:2f} dB.", atten);
  auto att = 20*log10(mag((int)ceil(Fint*1024/(Fin/2))));
  auto attmax = 20*log10(mag.head(ceil(Fint*1024/(Fin/2)))).minCoeff();
  msg("Attenuation at {:.2f} Hz: {:.2f} dB.", Fint, att);
  msg("Attenuation max. between 0 et {:.2f} Hz: {:.2f} dB.", Fint, attmax);
  msg("E.g. in linear scale : * {:.3f}", pow(10.0f,attmax/20));

  // Calcul repliement
  Fout = Fin / config.R;
  auto Fnyqout = Fout / 2;
  int id = trouve(fr > Fnyqout)[0];
  auto f = fr.head(id);
  auto m_main = lmag.head(id);

  auto nrep = (int) max(0.0f, (float) floor((mag.rows() - id) / id));
  nrep = min(nrep, 4);

  ArrayXXf m_alias = ArrayXXf::Zero(id,nrep);
  for(auto i = 0; i < nrep; i++)
  {
    //auto id2 = trouve(fr > (i+1)*Fnyqout);
    m_alias.col(i) = lmag.segment(id*(i+1), id);
    if((i % 2) == 0)
      m_alias.col(i).reverseInPlace();
  }


  {
    auto &g = res.fig_repliement;
    g.clear();
    auto c = g.plot(f * config.R, m_main, "g-", "Baseband spectrum");
    c.def_couleur({0,120,0});
    for(auto i = 0; i < nrep; i++)
    {
      auto c = g.plot(f * config.R, m_alias.col(i), "r-", "Aliasing #{}", i);
      auto v = (200.0f*i) / nrep;
      c.def_couleur({200,v,v});
    }

    g.def_pos_legende("se");
    g.titres("Output of CIC filter, including decimation",
             "Fréquence",
             "Atténuation (dB)");
    //g.afficher("Analyse CIC (repliement)");
  }

  {
    auto &g = res.fig_spectre_global;
    g.clear();
    auto c = g.plot(fr,lmag);
    c.def_remplissage(true, true, -150);
    g.plot(Fout/2, atten, "sr");
    g.def_rdi({-0.05, -150, 0.6, 150});
    g.titres("CIC filter / global view", "Fréquence", "Atténuation (dB)");
    //g.afficher("Analyse CIC (globale)");
  }
  {
    auto &g = res.fig_spectre_bf;
    auto c = g.plot(frz,20*(fmz+1e-10).log10());
    c.def_remplissage(true, true, -150);
    g.plot(Fout/2, atten, "sr");
    g.def_rdi({-0.05f/config.R, -150.0f, 3.0f*(0.05f+0.5f)/config.R, 150.0f});
    g.titres("CIC filter / view centered on the passband", "Fréquence", "Atténuation (dB)");
    //g.afficher("Analyse CIC (bande passante)");
  }

  res.h     = FRat<float>(h);
  res.fr    = fr;
  res.mag   = mag;
  res.nbits = (int) ceil(config.N * log2(config.R) - 1);
  msg("Number of additionnal bits needed for implementation: {}.", res.nbits);
  return res;
}



ArrayXf cic_freq(const CICConfig &config, const ArrayXf &f)
{
  ArrayXf mag = ArrayXf::Ones(f.rows());
  for(auto i = 0; i < f.rows(); i++)
  {
    float d = config.R*config.M*sin(π*f(i));
    if(d != 0)
      mag(i) = pow((float) abs(sin(config.R*config.M*π*f(i)) / d), (float) config.N);
  }
  return mag;
}


CICComp design_cic_comp(const CICConfig &config, float Fin, int R2, float Fcut, int ncoefs)
{
  CICComp res;
  //auto [h,bits,fr,fm] = cic_analyse(config,Fin,Fcut);
  auto analyse = cic_analyse(config, Fin, Fcut);

  // Frequency at output of CIC
  auto Fout = Fin / config.R;
  msg("R = {}, Fout = {:.2f} Hz, Fin = {:.2f} Hz.", config.R, Fout, Fin);

  // Locate frequency response after decimation
  // (after decimation, representable frequencies are 0 to Fout/2)
  auto idx = trouve(analyse.fr > Fout / 2)[0];
  msg("index = {} / {}.", idx, length(analyse.fr));
  auto fr2 = analyse.fr.segment(0,idx-1);
  auto fm2 = analyse.mag.segment(0,idx-1);

  msg("idx = {}", idx);

  //scf(2); clf();

  // Ideal response of the FIR compensation filter = 1 / CIC response
  ArrayXf idealc = 1.0f / fm2;

  // And, also, cut @ Fint Hz

  auto id = 1 + trouve(analyse.fr > Fcut)[0];

  msg("Fcut = {} Hz, id = {}", Fcut, id);

  tsd_assert((id > 1) && (id + 1 < analyse.fr.rows()));

  idealc(id-1)  = idealc(id-2) / 2;
  idealc(id)    = idealc(id-1) / 2; // Reduce the steepness of required frequency profile
  idealc(id+1)  = idealc(id) / 2;
  idealc.tail(idealc.rows() - (id+2)).setZero();

  {
    auto &f = res.fig_reponse_ideale;
    f.clear();
    f.plot(fr2, idealc);
    //f.titre("Réponse idéale");
  }



  // Make a fir filter from frequency sampling technique
  res.h = design_rif_freq(ncoefs, idealc);

  //    idx = trouve(fr * Fin > Fout / (2*R2), 1);
  //    fr2 = fr(1:idx);
  //    fm2 = fm(1:idx);

  // Round to 16 bits coefficients
  //hd = round(hd * 2^15);

  // Number of coefficients
 // auto nhd = hd.rows();
  //printf("FIR compensation filter (%d coefs, 16 bits quant):\n", nhd);
  //disp(hd);

  // Compute frequency response of the FIR compentation filter
  //hd = hd / 2^15;
  auto hc = FRat<float>::rif(res.h);
  auto [frc,fmc] = frmag<float>(res.h, 2048);

  // Fonction de transfert du filtre de compensation
  // decim(R) -> hc = hc(z^R)
  auto hc_lb = hc.horner(FRat<float>::z_power(config.R));

  // Fonction de transfert globale :
  // h(z) -> decim(R) -> hc(z) = h(z) * hc(z^R)
  auto hg = analyse.h * hc_lb;

  // Réponse filtre global
  auto [frg,fmg] = frmag(hg, 2048);
  // Réponse filtre de compensation
  auto [frg2,fmhc_lb] = frmag(hc_lb, 2048);

  /////////////////////////////////////////////////////
  // ( ) Tracé des différentes réponses
  /////////////////////////////////////////////////////
  auto f = res.fig_spectre_global;
  f.clear();
  f.plot(analyse.fr,20*log10(analyse.mag+1e-10f),"g-","CIC filter");
  f.plot(frg2*Fin,20*log10(fmhc_lb+1e-10), "m-", "Compensation filter");
  f.plot(frg*Fin,20*log10(fmg+1e-10),"b-","Global filter");
  f.titres("Decimation / wideband view", "Fréquence", "Attenuation (dB)");

  f = res.fig_spectre_bf;
  f.clear();
  f.plot(fr2,20*log10(fm2+1e-10),"g-","CIC filter");
  f.plot(frc*Fout,20*log10(fmc+1e-10),"m-","Compensation filter");
  auto fmaxz = 2 * Fcut;
  idx = trouve(frg*Fin > fmaxz)[0];
  auto frgz  = frg.head(idx);
  auto fmgz = fmg.head(idx);
  f.plot(frgz*Fin,20*log10(fmgz+1e-10),"b-","Global filter");
  f.titres("Decimation / passband view", "Fréquency", "Atténuation (dB)");

  f = res.fig_comp_rimp;
  f.clear();
  f.axes().supprime_decorations();
  f.plot(res.h, "|bo");//, "Filtre de compensation");

  return res;
}

namespace hidden {
auto filtre_cic1 = filtre_cic<float, int>;
//auto filtre_cic2 = filtre_cic<cfloat, int>;
}


}



