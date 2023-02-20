#include "tsd/telecom.hpp"
#include "tsd/filtrage/spline.hpp"
#include "tsd/vue.hpp"

using namespace tsd::vue;

namespace tsd::telecom {

const auto CLKREC_MODE_SAFE = non;

struct TedMM: Ted
{
  TedMM()
  {
    osf  = 1;
    npts = 2;
  }
  float calcule(cfloat x0, cfloat x1, cfloat x2)//const ArrayXcf &x)
  {
    //ArrayXf xr = x.real(), xi = x.imag();
    //auto dr = xr.sign(), di = xi.sign(); // Décision (A FAIRE : appeler un slicer en fonction de la modulation)
    //retourne -(dr(0) * xr(1) - dr(1) * xr(0)) / 2.0f - (di(0) * xi(1) - di(1) * xi(0)) / 2.0f;

    //auto dr = x2.real().sign(), di = x2.imag().sign(); // Décision (A FAIRE : appeler un slicer en fonction de la modulation)
    //retourne -(dr(0) * xr(1) - dr(1) * xr(0)) / 2.0f - (di(0) * xi(1) - di(1) * xi(0)) / 2.0f;
    retourne 0;
  }
};

struct TedEL: Ted
{
  TedEL()
  {
    osf  = 2;
    npts = 4;
  }
  float calcule(cfloat x0, cfloat x1, cfloat x2)//const ArrayXcf &x)
  {
    //ArrayXf xr = x.real(), xi = x.imag();
    //retourne -xr(1) * (xr(2) - xr(0))/2 - xi(1) * (xi(2) - xi(0)) / 2;
    //retourne -xr(1) * (xr(2) - xr(0))/2 - xi(1) * (xi(2) - xi(0)) / 2;

    // TODO
    retourne 0;
  }
};


#if 0
static float gardner_reel(float x0, float x1, float x2)
{
  // si petit écart de temps :
  // 2 * tau
  retourne ((std::copysign(1.0f, x2) - std::copysign(1.0f, x0)) * x1);
  //retourne (x2 - x0) * x1;
  //retourne ((std::copysign(1.0f, x2) - std::copysign(1.0f, x0)) * (x1 - (x0 + x2) / 2.0f));
}
#endif


// Pb avec TedGardner en QAM (surement pareil en M-ASK)
struct TedGardner: Ted
{
  TedGardner()
  {
    osf  = 2;
    npts = 3;
  }
  float calcule(cfloat x0, cfloat x1, cfloat x2)//const ArrayXcf &x)
  {
    retourne real((x2 - x0) * conj(x1));

    //float err_I = gardner_reel(x(0).real(), x(1).real(), x(2).real());
    //float err_Q = gardner_reel(x(0).imag(), x(1).imag(), x(2).imag());
    //retourne err_I + err_Q;
  }
};


sptr<Ted> ted_init(TedType type)
{
  si(type == TedType::GARDNER)
    retourne std::make_shared<TedGardner>();
  sinon si(type == TedType::MM)
    retourne std::make_shared<TedMM>();
  sinon si(type == TedType::EARLY_LATE)
    retourne std::make_shared<TedEL>();

  msg_erreur("Ted init");
  retourne sptr<Ted>();
}





struct ClockRec: FiltreGen<cfloat>
{

  float phase;
  entier K1 = 0, K2 = 0, cnt = 0;
  Veccf fenetre_x;
  float gain = 1.0f;
  entier index_fenetre_x = 0;

  // pour la TED
  cfloat x0 = 0.0f, x1 = 0.0f, x2 = 0.0f;

  ClockRecConfig config;

  ClockRec(const ClockRecConfig &config)
  {
    this->config = config;
    phase = ((float) config.osf) / 2;
    K1   = config.osf;

    si(!config.itrp)
      echec("clock rec : intepolateur non spécifié.");

    si(!config.ted)
      echec("clock rec : ted non spécifié.");

    K2   = config.ted->osf;
    tsd_assert(K2 > 0);

    cnt  = 0;
    // Sliding window pour the interpolator
    fenetre_x   = Veccf::zeros(config.itrp->K);
    // Sliding window pour the TED
    //fenetre_ted = ArrayXcf::Zero(config.ted->npts);
    // conversion tc en période d'échantillonnage
    this->gain = K1 * (1 - exp(-1/(config.tc * K1)));

    x0 = x1 = x2 = 0.0f;

    msg("clock rec init: K1 (osf) = {}, K2 = {}, npts ted = {}, npts itrp = {}. phase initiale = {}, tc = {}", config.osf, config.ted->osf, config.ted->npts, config.itrp->K, phase, config.tc);
    msg("Fréquence d'entrée : [K1 = {}] * fsymb", K1);
    msg("Fréquence de travail TED : [K2 = {}] * fsymb", config.ted->osf);
    msg("Fréquence de sortie : fsymb");
  }

  /*inline void maj_fenetre(ArrayXcf &wnd, cfloat x)
  {
    auto n = wnd.rows();
    wnd.head(n-1) = wnd.tail(n-1).eval();
    wnd(n-1) = x;
  }*/

  void step(const Veccf &x, Veccf &y)
  {
    // 3 fréquences différentes :
    // F d'entrée
    // F de fonctionnement du TED
    // F de sortie

    // Fréquence de base = Fsymbole (1 Hz)
    // Fréquence d'entrée : osf * Fsymbole
    // Fréquence de travail de la TED : K2 * Fsymbole

    //ArrayXf E, MU, DEC;

    //ArrayXf coarse_rssi = ArrayXf::Ones(x.rows());

    //si(argn(2) > 2)
    //    coarse_rssi = varargin(1);
    //end;

    //std::vector<cfloat> res(n);



    /*si(coarse_rssi.rows() != x.rows())
    {
      erreur("clock_rec_process: x and coarse_rssi must be of the same length.");
      retourne;
    }*/

    // Phase est un compteur exprimé en pas de Tsymb / K1 = Tsymb / OSF
    //

    entier n = x.rows();

    entier nout_max = 2 + n / config.osf;
    entier oindex = 0;

    Veccf res(nout_max);
    Vecf vphase, verr, vphase0, vitrp, vdec;
    std::vector<entier> idx, idx_inter;
    std::vector<float> t_inter, v_inter;
    std::vector<float> t_inter2, v_inter2;

    si(config.debug_actif)
    {
      vitrp = vphase0 = vdec = verr = vphase = Vecf::zeros(n);
    }


    float ph0 = 0;

    pour(auto i = 0; i < n; i++)
    {
      si(config.debug_actif)
      {
        vphase(i)  = phase;
        vphase0(i) = ph0;
      }

      //msg("i = {} / {}, phase = {}", i, n, phase);

      //maj_fenetre(fenetre_x, x(i));

      fenetre_x(index_fenetre_x) = x(i);
      index_fenetre_x = (index_fenetre_x + 1) % config.itrp->K;

      // Requiert: phase >= 1
      phase--;
      si(phase > 1)
        continue;

      // Requiert: phase >= 0
      si constexpr (CLKREC_MODE_SAFE)
      {
        si(phase < 0)
        {
          msg_erreur("clock rec : phase négative ({}). Incrément phase = {}", phase, ((float) K1) / K2);
          phase = 0;
        }
      }


      /// x entrant échantilloné ofs * fsymb
      // -> sortie d'interpolateur :
      //    2 * fsymb, après accrochage : aligné

      // Ici on est à la fréquence de la TED
      soit interpol = config.itrp->step(fenetre_x, index_fenetre_x, phase);

      si constexpr (CLKREC_MODE_SAFE)
      {
        si(std::isnan(interpol.real()) || std::isnan(interpol.imag()))
          msg_erreur("Itrp : nan, fen itrp = {}.", fenetre_x);
      }

      phase += ((float) K1) / K2; // Lecture au rythme de la TED

      // Idée : pour la TED : récupérer la fenêtre x
      // (prendre un échantillon / K2/K1)
      //maj_fenetre(fenetre_ted, interpol);
      x0 = x1;
      x1 = x2;
      x2 = interpol;

      //idx_inter.push_back(i);


      si(config.debug_actif)
      {
        vitrp(i) = interpol.real();
        t_inter.push_back(i + phase - config.itrp->delais * ((float) K1) / K2);
        v_inter.push_back(interpol.real());
      }

      si(cnt == K2 - 1)
      {
        tsd_assert(oindex < nout_max);
        res(oindex++) = interpol;

        // N'appelle pas la TED à chaque fois
        // (au même rythme que les éch. de sortie)
        //soit e = config.ted->calcule(x0, x1, x2);//fenetre_ted);//  / coarse_rssi(i));

        soit e = real((x2 - x0) * conj(x1));

        si constexpr (CLKREC_MODE_SAFE)
        {
          si(std::isnan(e))
          {
            msg_erreur("Ted : nan, fenetre = {}, {}, {}", x0, x1, x2);
            e = 0;
          }
        }

        // Filtre RII du premier ordre
        // mu est exprimé en : nombre de samples d'entrée
        // e : en multiple de la période symbole
        soit dec = gain * e;

        //infos("e = %f, dec = %f", e, dec);

        // Décalage maximum = 0.25 symboles
        dec = std::clamp(dec, -K1/4.0f, K1/4.0f);

        //infos("clamp = %f", dec);
        si(config.debug_actif)
        {
          t_inter2.push_back(i + phase - config.itrp->delais * ((float) K1) / K2);
          v_inter2.push_back(interpol.real());
          verr(i) = 100 * e;
          vdec(i) = 100 * dec;
        }

        phase -= dec;
        ph0   -= dec;
        cnt   = 0;
      }
      sinon
        cnt++;
      //cnt = (cnt + 1) % K2;
    }


    //y = Eigen::Map<ArrayXcf>(res.data(), res.size());
    y = res.head(oindex);

    //infos("clock rec : %d in --> %d out", x.rows(), y.rows());

    si(config.debug_actif)
    {
      Figure f2;
      f2.plot(vitrp);
      f2.afficher("itrp");

      Figures figs;
      //soit nl = 3, nc = 1, ids = 1;
      soit f = figs.subplot();
      f.plot(x, "", "x");

      soit a = f.canva();

      a.set_couleur({100,0,0});
      //
      //pour(soit i: idx_inter)


      pour(auto i = 0u; i < t_inter.size(); i++)
      {
        a.marqueur({t_inter[i], v_inter[i]}, tsd::vue::Marqueur::CERCLE, 7);
      }


      a.set_remplissage(oui, {150,0,0});
      //pour(auto i: t_inter2)//idx)
      pour(auto i = 0u; i < t_inter2.size(); i++)
      {
        a.marqueur({t_inter2[i], v_inter2[i]}, tsd::vue::Marqueur::CERCLE, 7);
      }


      a.set_couleur({100,100,0});
      //pour(auto i: idx)
      pour(auto i : t_inter2)
      {
        a.ligne(i, -1, i, 1.0);
      }

      //infos("Infos GCA clk rec :");
      //a.dump_infos();

      f = figs.subplot();
      vphase0 *= (100.0f / K1);
      soit [vpmin, vpmax] = vphase0.valeurs_minmax();

      f.plot(vphase0, "-m", "Phase (% de période symbole)");

      //f.subplot(nl,nc,ids++);
      //f.plot(vphase  * (100.0f / K1), "-b", "Phase (avec éch)");


      a = f.canva();
      pour(auto i: t_inter2)
      {
        a.set_couleur({100,100,0});
        a.ligne(i, vpmin, i, vpmax);
      }

      f = figs.subplot();
      f.plot(verr / K1, "-r", "Erreur (% de période symbole)");

      //f.subplot(nl,nc,ids++);
      //f.plot(vdec / K1, "-r", "Décalage (% de période symbole)");


      figs.afficher("Clock recovery");

    }
  }
};





sptr<FiltreGen<cfloat>> clock_rec_init(const ClockRecConfig &config)
{
  retourne std::make_shared<ClockRec>(config);
}


// MENGALI, page 374
struct ClockRec2: FiltreGen<cfloat>
{

  float phase;
  entier K1 = 0;
  float gain = 1.0f;
  Veccf fenetre_x, fenetre_dx;


  sptr<FiltreGen<cfloat>> fa, fda;


  ClockRecConfig config;

  ClockRec2(const ClockRecConfig &config)
  {
    this->config = config;
    phase = ((float) config.osf) / 2;
    this->K1   = config.osf;

    si(!config.itrp)
      echec("clock rec : intepolateur non spécifié.");

    // Sliding window pour the interpolator
    fenetre_x  = Veccf::zeros(config.itrp->K);
    fenetre_dx = Veccf::zeros(config.itrp->K);
    // conversion tc en période d'échantillonnage
    this->gain = K1 * (1 - exp(-1/(config.tc * K1)));

    fa = tsd::filtrage::filtre_rif<float,cfloat>(config.h_fa);

    entier n = config.h_fa.rows();
    Vecf dh(n+1);
    dh(0) = config.h_fa(0);
    pour(auto i = 1; i < n; i++)
      dh(i) = config.h_fa(i) - config.h_fa(i-1);
    dh(n) = -config.h_fa(n-1);

    si(config.debug_actif)
    {
      Figures f;
      f.subplot(211).plot(config.h_fa, "b-o", "Filtre adapté");
      f.subplot(212).plot(dh, "g-o", "Dérivée du filtre adapté");
      f.afficher("Filtre adapté clk rec");
    }

    fda = tsd::filtrage::filtre_rif<float,cfloat>(dh);



    msg("clock rec 2 init: K1 (osf) = {}, npts itrp = {}. phase initiale = {}, tc = {}", config.osf, config.itrp->K, phase, config.tc);
    msg("Fréquence d'entrée : [K1 = {}] * fsymb", K1);
    msg("Fréquence de sortie : fsymb");
  }

  void maj_fenetre(Veccf &wnd, cfloat x)
  {
    soit n = wnd.rows();
    wnd.head(n-1) = wnd.tail(n-1).eval();
    wnd(n-1) = x;
  }

  void step(const Veccf &x, Veccf &y)
  {
    std::vector<cfloat> res;

    soit xf  = fa->step(x);
    soit xdf = fda->step(x);

    // Phase est un compteur exprimé en pas de Tsymb / K1 = Tsymb / OSF
    //

    soit n = xf.rows();

    tsd_assert(n == xdf.rows());

    Vecf vphase(n), verr = Vecf::zeros(n), vphase0(n), vitrp = Vecf::zeros(n), vdec = Vecf::zeros(n);
    std::vector<entier> idx, idx_inter;

    std::vector<float> t_inter, v_inter;


    float ph0 = 0;

    pour(auto i = 0; i < n; i++)
    {
      vphase(i) = phase;
      vphase0(i) = ph0;

      //msg("i = {} / {}, phase = {}", i, n, phase);

      maj_fenetre(fenetre_x,  xf(i));
      maj_fenetre(fenetre_dx, xdf(i));

      // Requiert: phase >= 1
      phase--;
      si(phase > 1)
        continue;

      // Requiert: phase >= 0
      si(phase < 0)
      {
        msg_erreur("clock rec : phase négative ({}). Incrément phase = {}", phase, K1);
        phase = 0;
      }

      // Ici on est à la fréquence de la TED
      soit interpol    = config.itrp->step(fenetre_x, 0, phase);
      soit interpol_dx = config.itrp->step(fenetre_dx, 0, phase);

      si(std::isnan(interpol.real()) || std::isnan(interpol.imag()))
      {
        //msg_erreur("Itrp : nan, fen itrp = {}.", fenetre_x);
        msg_erreur("Itrp : nan");
      }

      phase += K1; // Lecture au rythme symbole

      si(config.debug_actif)
      {
        vitrp(i) = interpol.real();
        t_inter.push_back(i + phase - config.itrp->delais * ((float) K1));
        v_inter.push_back(interpol.real());
      }


      res.push_back(interpol);

      // Appele pas la TED à chaque fois
      // (au même rythme que les éch. de sortie)
      //auto e = config.ted->calcule(fenetre_ted);//  / coarse_rssi(i));

      // Non -> interpol --> décision

      soit e = (interpol * interpol_dx).real();

      /////////////////////////////////////////////////////////////////
      // PB : suppose que la carrier rec a été faite avant !!!
      //float symb = std::copysign(1.0f, interpol.real());
      //soit e = (symb * interpol_dx).real();

      si(std::isnan(e))
      {
        msg_erreur("Ted : nan");//, x = {}, dx = {}", interpol, interpol_dx);
        e = 0;
      }

      // Filtre IIR du premier ordre
      // mu est exprimé en : nombre de samples d'entrée
      // e : en multiple de la période symbole
      soit dec = gain * e;

      verr(i) = 100 * e;

      //infos("e = %f, dec = %f", e, dec);

      // Décalage maximum = 0.25 symboles
      dec = std::clamp(dec, -K1/4.0f, K1/4.0f);

      //infos("clamp = %f", dec);

      vdec(i) = 100 * dec;

      phase -= dec;
      ph0 -= dec;
    }

    y.resize(res.size());
    memcpy(y.data(), res.data(), res.size() * sizeof(cfloat));

    //y = Eigen::Map<ArrayXcf>(res.data(), res.size());

    //infos("clock rec : %d in --> %d out", x.rows(), y.rows());

    si(config.debug_actif)
    {
      Figure f2;
      f2.plot(vitrp);
      f2.afficher("itrp");

      Figures figs;
      soit f = figs.subplot();
      f.plot(x, "", "x");

      soit a = f.canva();

      a.set_couleur({100,0,0});
      //
      //pour(auto i: idx_inter)

      a.set_remplissage(oui, {150,0,0});
      pour(auto i = 0u; i < t_inter.size(); i++)
        a.marqueur({t_inter[i], v_inter[i]}, tsd::vue::Marqueur::CERCLE, 7);
      a.set_couleur({100,100,0});
      pour(auto i : t_inter)
        a.ligne(i, -1, i, 1.0);

      //infos("Infos GCA clk rec :");
      //a.dump_infos();

      f = figs.subplot();
      vphase0 *= (100.0f / K1);
      soient [vpmin, vpmax] = vphase0.valeurs_minmax();

      f.plot(vphase0, "-m", "Phase (% de période symbole)");

      //f.subplot(nl,nc,ids++);
      //f.plot(vphase  * (100.0f / K1), "-b", "Phase (avec éch)");


      a = f.canva();
      pour(auto i: t_inter)
      {
        a.set_couleur({100,100,0});
        a.ligne(i, vpmin, i, vpmax);
      }

      figs.subplot().plot(verr / K1, "-r", "Erreur (% de période symbole)");

      //f.subplot(nl,nc,ids++);
      //f.plot(vdec / K1, "-r", "Décalage (% de période symbole)");


      figs.afficher("Clock recovery");
    }
  }
};

sptr<FiltreGen<cfloat>> clock_rec2_init(const ClockRecConfig &config)
{
  retourne std::make_shared<ClockRec2>(config);
}




}
