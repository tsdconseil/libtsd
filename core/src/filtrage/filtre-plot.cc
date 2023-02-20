#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include <iostream>

using namespace tsd::vue;
using namespace std;

namespace tsd::filtrage {


template<typename T>
  void plot_plz(tsd::vue::Figure &fig, const Vecteur<T> &h)
{
  plot_plz(fig, FRat<T>::rif(h));
}


template<typename T>
void plot_plz(Figure &fig, const FRat<T> &h)
{
  soit [zéros,pôles] = h.roots();

  fig.axes().set_isoview(oui);
  fig.axes().get_config().grille_mineure.afficher = non;
  fig.def_pos_legende("se");
  fig.def_rdi_min({-1.1, -1.1, 2.1, 2.1});

  si(zéros.rows() > 0)
    fig.plot_iq(zéros, "sb", "zéros");
  si(pôles.rows() > 0)
    fig.plot_iq(pôles, "sr", "pôles");

  soit c = fig.canva_pre();
  soit gris = Couleur{200,200,200,128};
  c.set_couleur(gris);
  c.set_epaisseur(3);
  c.set_remplissage(oui, Couleur::Vert.eclaircir(0.05));
  c.cercle(0, 0, 1.0);
  c.set_remplissage(non);
  c.set_epaisseur(1);
  pour(auto k = 0.2; k <= 0.8; k += 0.2)
  {
    c.set_couleur(gris);
    c.cercle(0, 0, k);
  }
  c.set_align(Align::CENTRE, Align::CENTRE);
  //c.set_couleur()
  pour(auto θ = 0.0; θ <= 315; θ += 45)
  {
    Pointf p{cos(deg2rad(θ)), sin(deg2rad(θ))};
    c.ligne({0, 0}, p);
    Pointf p2 = p;
    p2.x *= 1.15;
    p2.y *= 1.15;
    c.texte(p2, fmt::format("{}°", (int) θ));
  }

  //fig.titre("Pôles et zéros");
}



template<typename T>
  tsd::vue::Figures affiche_filtre(const FRat<T> &h, float fe)
{
  Figures f;
  soit [fr, mag] = frmag(h);
  soit y = repimp(h);

  soit c = f.subplot().plot(real(y), y.rows() < 50 ? "o|g" : "-g");
  c.def_epaisseur(2);
  f.gcf().titres("Réponse impulsionnelle", "Echantillons");

  soit lmag = mag2db(mag);

  // Nettoyage des valeurs négligeables
  //pour(auto i = 0; i < lmag.rows(); i++)
    //si(abs(lmag(i)) < 1e-5)
      //lmag(i) = 0;

  string axe_fréq = "Fréquence";
  si(fe != 1)
    axe_fréq += " (Hz)";

  c = f.subplot().plot(fr, lmag);
  c.def_remplissage(oui, non);
  c.def_epaisseur(2);
  f.gcf().titres("Réponse fréquentielle (log)", axe_fréq, "Atténuation (dB)");

  retourne f;
}



Figures plot_filtre(const Vecf &h, bouléen complet, float fe)
{
  retourne plot_filtre(FRat<float>::rif(h), complet, fe);
}

// A faire : analyse stabilité...
template<typename T>
Figures plot_filtre(const FRat<T> &h, bouléen complet, float fe)
{
  si(!complet)
    retourne affiche_filtre(h, fe);

  // Nb points d'analyse
  soit N = 512;
  // Calcul de la réponse fréquentielle
  soit [fr, hzm] = frmag(h, N);

  Vecf pzm, gzm;
  tie(ignore, pzm) = frphase(h, N);
  tie(ignore, gzm) = frgroup(h, N);

  Figures fs;

  string axe_fréq = "Fréquence";
  si(fe != 1)
    axe_fréq += " (Hz)";

  {
    soit f = fs.subplot(241);
    f.plot(fr*fe, hzm);
    f.titres("Réponse fréquentielle", axe_fréq, "Atténuation");

    f = fs.subplot(242);
    f.plot(fr*fe, log10(hzm) * 20);
    f.titres("Réponse fréquentielle (log)", axe_fréq, "Atténuation (dB)");

    f = fs.subplot(243);
    f.plot(fr*fe, pzm * 180.0f / π);
    f.titres("Phase", axe_fréq, "Phase (degrés)");

    f = fs.subplot(244);
    f.plot(fr*fe,gzm);
    f.titres("Temps de groupe", axe_fréq, "Délais");
  }


  soit y = repimp(h);
  soit npts = h.numer.coefs.rows();

  si(!h.est_rif())
    npts *= 20; // Réponse infinie

  {
    soit f = fs.subplot(245);
    f.plot(y, y.rows() < 30 ? "o|b" : "-b");
    f.titres("Réponse impulsionnelle", "Echantillons");
  }

  {
    soit f = fs.subplot(246);
    soit x = Veccf::ones(2 * npts);
    soit y = filtrer(h, x);
    f.plot(real(y), y.rows() < 30 ? "s-b" : "-b");
    f.titres("Réponse échelon", "Echantillons");
  }

  {
    soit f = fs.subplot(247);
    plot_plz(f, h);
  }

  retourne fs;
}


template Figures plot_filtre(const FRat<float> &h, bouléen, float fe);
template Figures plot_filtre(const FRat<cfloat> &h, bouléen, float fe);
template void plot_plz(tsd::vue::Figure &fig, const Vecteur<float> &h);

}
