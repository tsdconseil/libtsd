#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include <iostream>

using namespace tsd::vue;


namespace tsd::filtrage {




template<typename T>
void plot_plz(Figure &fig, const FRat<T> &h)
{
  //msg("Calcul racines, nb coefs numer = {}, denom = {}...", (int) h.numer.coefs.rows(), (int) h.denom.coefs.rows());
  //std::cout << "Numer : " << h.numer.coefs.transpose() << std::endl;
  auto zeros = h.numer.roots();


  //std::cout << "Dénom : " << h.denom.coefs.transpose() << std::endl;
  auto poles = h.denom.roots();

  //msg("racines ok.");

  //std::cout << "Zéros : " << zeros.transpose() << std::endl;
  //std::cout << "Poles : " << poles.transpose() << std::endl;

  fig.axes().set_isoview(true);
  //fig.canva().set_couleur(tsd::vue::Couleur{0,0,0});
  //fig.canva().set_epaisseur(1);

  fig.def_pos_legende("se");
  fig.def_rdi_min({-1.1, -1.1, 2.1, 2.1});
  if(zeros.rows() > 0)
    fig.plot_iq(zeros, "*b", "Zéros");
  if(poles.rows() > 0)
    fig.plot_iq(poles, "sr", "Pôles");

  fig.canva().set_couleur(tsd::vue::Couleur{0,0,0});
  fig.canva().set_epaisseur(1);
  fig.canva().cercle(0, 0, 1.0);
  fig.titre("Pôles et zéros");
}



template<typename T>
ArrayXf repimp(const FRat<T> &h, int npts)
{
  if(npts == -1)
  {
    npts = h.numer.coefs.rows();

    if(!h.est_fir())
      npts *= 20; // Réponse infinie
  }

  ArrayXcf x = ArrayXcf::Zero(npts);
  x(0) = 1;

  auto y = filtrer(h, x);
  if(!h.est_fir())
  {
    while((y.rows() > 1) && (abs(y(y.rows()-1)) < y.abs().mean() * 0.01))
      y = y.segment(0,y.rows()-2).eval();
  }
  return y.real();
}

template<typename T> ArrayXf repimp(const Vecteur<T> &h, int npts)
{
  return repimp(FRat<T>::rif(h), npts);
}



template<typename T>
  tsd::vue::Figures affiche_filtre(const FRat<T> &h, float fe)
{
  Figures f;
  auto [fr, mag] = frmag(h);

  ArrayXf y = repimp(h);

  f.subplot().plot(y, "|bo", "Réponse impulsionnelle");
  f.gcf().axes().supprime_decorations();

  ArrayXf lmag = 20*log10(mag);

  // Nettoyage des valeurs négligeables
  for(auto i = 0; i < lmag.rows(); i++)
    if(abs(lmag(i)) < 1e-5)
      lmag(i) = 0;

  auto c = f.subplot().plot(fr, lmag, "", "Réponse fréquentielle (log)");
  c.def_remplissage(true, false);
  return f;
}

Figures affiche_filtre(const ArrayXf &h, float fe)
{
  return affiche_filtre(FRat<float>::rif(h), fe);
}

Figures analyse_filtre(const ArrayXf &h, float fe)
{
  return analyse_filtre(FRat<float>::rif(h), fe);
}

// A faire : analyse stabilité...
template<typename T>
Figures analyse_filtre(const FRat<T> &h, float fe)
{
  // Nb points d'analyse
  auto N = 512u;
  // Calcul de la réponse fréquentielle
  auto [fr, hzm] = frmag(h, N);

  ArrayXf pzm, gzm;
  std::tie(std::ignore, pzm) = frphase(h, N);
  std::tie(std::ignore, gzm) = frgroup(h, N);

  Figures fs;

  {
    auto f = fs.subplot(241);
    f.plot(fr*fe,hzm);
    f.titres("Réponse fréquentielle", "Fréquence (Hz)", "Atténuation");

    f = fs.subplot(242);
    f.plot(fr*fe,20.0f*hzm.log10());
    f.titres("Réponse fréquentielle (log)", "Fréquence (Hz)", "Atténuation (dB)");

    f = fs.subplot(243);
    f.plot(fr*fe,pzm * 180.0f / π);
    f.titres("Phase", "Fréquence (Hz)", "Phase (degrés)");

    f = fs.subplot(244);
    f.plot(fr*fe,gzm);
    f.titres("Temps de groupe", "Fréquence (Hz)", "Délais");
  }


  auto y = repimp(h);

  auto npts = h.numer.coefs.rows();

  if(!h.est_fir())
    npts *= 20; // Réponse infinie

  //msg("rimp...");
  {
    auto f = fs.subplot(245);

    //ArrayXcf x = ArrayXcf::Zero(npts);
    //x(0) = 1;

    // Pb si Frat complex !!

    //msg("filtrer...");
    //auto y = filtrer(h, x);
    //msg("ok.");
    //if(!h.est_fir())
    //{
    //  while((y.rows() > 1) && (abs(y(y.rows()-1)) < y.abs().mean() * 0.01))
    //    y = y.segment(0,y.rows()-2).eval();
    //}

    // 3)
    f.plot(y.real(), y.rows() < 30 ? "s-b" : "-b");
    f.titres("Réponse impulsionnelle", "Echantillons");
  }

  {
    auto f = fs.subplot(246);

    ArrayXcf x = ArrayXf::Ones(2 * npts);
    x(0) = 1;

    auto y = filtrer(h, x);
    f.plot(y.real(), y.rows() < 30 ? "s-b" : "-b");
    f.titres("Réponse échelon", "Echantillons");
  }

  //msg("plz..");
  {
    auto f = fs.subplot(247);
    plot_plz(f, h);
  }


  return fs;


  //f.afficher();

  /*
  subplot(236);
  [tg,fr] = group(128,sl);
  plot(fr,tg);
  xtitle("Temps de groupe")
  */

}


template Figures analyse_filtre(const FRat<float> &h, float fe);
template Figures analyse_filtre(const FRat<cfloat> &h, float fe);
template Figures affiche_filtre(const FRat<float> &h, float fe);
template Figures affiche_filtre(const FRat<cfloat> &h, float fe);
template ArrayXf repimp(const Vecteur<float> &h, int npts);
template ArrayXf repimp(const Vecteur<cfloat> &h, int npts);

}
