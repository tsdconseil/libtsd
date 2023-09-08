#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include <iostream>

using namespace tsd::vue;
using namespace std;

namespace tsd::filtrage {


void plot_plz(Figure fig, const Design &d, bouléen cmap)
{
  Veccf zéros, pôles;

  si(!d.est_rif && d.est_complexe)
    std::tie(zéros, pôles) = d.frat_c.roots(); // Expect root representation
  sinon si(!d.est_rif)
    std::tie(zéros, pôles) = d.frat.roots();
  sinon
    std::tie(zéros, pôles) = FRat<float>::rif(d.coefs).roots();

  soit h = d.complex_frat();

  fig.axes().set_isoview(oui);
  soit &cfg = fig.axes().get_config();
  cfg.grille_mineure.afficher = non;
  si(cmap)
  {
    cfg.couleur_avant_plan = Couleur::Gris.eclaircir(0.5);
    cfg.axe_horizontal.afficher = non;
    cfg.axe_vertical.afficher   = non;
  }
  fig.def_pos_legende("se");
  fig.def_rdi_min({-1.1, -1.1, 2.1, 2.1});


  si(cmap)
  {
    soit nx = 101, ny = 101;
    Tabf Z(nx,ny);

    pour(auto ix = 0; ix < nx; ix++)
    {
      pour(auto iy = 0; iy < ny; iy++)
      {
        cfloat z;
        z.real(-1.5+(3.0*ix)/nx);
        z.imag(-1.5+(3.0*iy)/ny);
        Z(ix,iy) = abs(h.horner(z));
      }
    }

    // Gain = +/- 10 dB par rapport à 0 dB
    // Soit [0,1 - 10]
    // [max-10 dB, max]

    Z = Z.cwiseMax(0.001f);
    Z = Z.cwiseMin(1000.0f);
    Z = log10(Z); // -> [-3, 3]
    Z = (Z + 3) / 6; // [0, 1]

    fig.canva_pre().plot_cmap(Z, {-1.5, -1.5, 3.0, 3.0}, cmap_parse("pm"));
  }

  string cz = "b", cp = "r";

  si(cmap)
  {
    cz = "w";
    cp = "y";
  }

  si(!cmap)
  {
    si(zéros.rows() > 0)
      fig.plot_iq(zéros, "s" + cz, est_fr() ? "zéros" : "zeros");
    si(pôles.rows() > 0)
      fig.plot_iq(pôles, "s" + cp, est_fr() ? "pôles" : "poles");
  }

  soit max_mag = max(abs(zéros).valeur_max(), abs(pôles).valeur_max());

  si(fig.axes().a_décorations())
  {
    soit c = fig.canva_pre();
    soit gris = Couleur{200,200,200,128};
    c.set_couleur(gris);
    c.set_epaisseur(3);
    si(!cmap)
      c.set_remplissage(oui, Couleur::Jaune.eclaircir(0.05));
    c.cercle(0, 0, 1.0);
    c.set_remplissage(non);
    c.set_epaisseur(1);
    pour(auto k = 0.2; k <= 0.8; k += 0.2)
    {
      c.set_couleur(gris);
      c.cercle(0, 0, k);
    }
    c.set_align(Align::CENTRE, Align::CENTRE);
    pour(auto θ = 0.0; θ <= 315; θ += 45)
    {
      Pointf p{cos(deg2rad(θ)), sin(deg2rad(θ))};
      c.ligne({0, 0}, p);
      p.x *= max_mag * 1.15;
      p.y *= max_mag * 1.15;
      c.texte(p, sformat("{}°", (entier) θ));
    }
  }
}

void plot_rech(Figure f, const Design &d)
{
  soit y = repech(d);
  soit c = f.plot(real(y), y.rows() < 128 ? "o|g" : "-g");
  c.def_epaisseur(2);
  f.titres(est_fr() ? "Réponse à un échelon" : "Step response",
      est_fr() ? "Echantillons" : "Samples");
}

void plot_rimp(Figure f, const Design &d)
{
  soit y = repimp(d);
  soit c = f.plot(real(y), y.rows() < 128 ? "o|g" : "-g");
  c.def_epaisseur(2);
  f.titres(est_fr() ? "Réponse impulsionnelle" : "Impulse response",
      est_fr() ? "Echantillons" : "Samples");
}

void plot_frmag(Figure f, const Design &d, bouléen mode_log, float fe)
{
  soit [fr, mag] = frmag(d);
  si(mode_log)
    mag = mag2db(mag);

  string axe_fréq = est_fr() ? "Fréquence" : "Frequency";
  si(fe != 1)
    axe_fréq += " (Hz)";

  string suffixe = mode_log ? " (log)" : " (lin)",
         unit    = mode_log ? " (dB)" : "";

  soit c = f.plot(fr, mag);
  c.def_remplissage(oui, non);
  c.def_epaisseur(2);
  f.titres((est_fr() ?
      "Réponse fréquentielle" : "Frequency response") + suffixe,
      axe_fréq, (est_fr() ? "Atténuation" : "Attenuation") + unit);
}

void plot_frphase(Figure f, const Design &d, float fe)
{
  soit [fr, phase] = frphase(d);

  string axe_fréq = est_fr() ? "Fréquence" : "Frequency";
  si(fe != 1)
    axe_fréq += " (Hz)";

  soit c = f.plot(fr, 180 * phase / π_f);
  c.def_remplissage(oui, non);
  c.def_epaisseur(2);
  f.titres((est_fr() ?
      "Réponse en phase" : "Phase response"),
      axe_fréq, "Phase (°)");
}


tsd::vue::Figures affiche_filtre(const Design &d, float fe)
{
  Figures f;

  plot_rimp(f.subplot(), d);
  plot_frmag(f.subplot(), d, oui, fe);

  retourne f;
}

// A faire : analyse stabilité...
Figures plot_filtre(const Design &d, bouléen complet, float fe)
{
  si(!complet)
    retourne affiche_filtre(d, fe);

  // Nb points d'analyse
  soit N = 512;
  // Calcul de la réponse fréquentielle
  soit [fr, hzm] = frmag(d, N);

  Vecf pzm, gzm;
  tie(ignore, pzm) = frphase(d, N);
  tie(ignore, gzm) = frgroup(d, N);

  Figures fs;

  string axe_fréq = "Fréquence";
  si(fe != 1)
    axe_fréq += " (Hz)";

  {
    //soit f = fs.subplot(241);

    plot_frmag(fs.subplot(241), d, non, fe);

    //f.plot(fr*fe, hzm);
    //f.titres("Réponse fréquentielle", axe_fréq, "Atténuation");

    plot_frmag(fs.subplot(242), d, oui, fe);

    //f = fs.subplot(242);
    //f.plot(fr*fe, log10(hzm) * 20);
    //f.titres("Réponse fréquentielle (log)", axe_fréq, "Atténuation (dB)");

    soit f = fs.subplot(243);
    f.plot(fr*fe, pzm * 180.0f / π);
    f.titres("Phase", axe_fréq, "Phase (degrés)");

    f = fs.subplot(244);
    f.plot(fr*fe,gzm);
    f.titres("Temps de groupe", axe_fréq, "Délais");
  }




  plot_rimp(fs.subplot(245), d);

  {
    soit y1 = repimp(d);
    soit f = fs.subplot(246);
    soit x = Veccf::ones(2 * y1.rows());
    soit y = filtrer(d, x);
    f.plot(real(y), y.rows() < 30 ? "s-b" : "-b");
    f.titres("Réponse échelon", "Echantillons");
  }

  plot_plz(fs.subplot(247), d);

  retourne fs;
}

}
