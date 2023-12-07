#include "tsd/vue/misc.hpp"


namespace tsd::vue {


  void FigureEvenements::maj(const FigureEvenementModele &content)
  {
    this->content = content;
  }

  void FigureEvenementModele::Ligne::concaténation()
  {
    vector<Evt> evts2;
    entier état = 0;
    Evt evt_en_cours;

    pour(auto i = 0; i < (entier) evts.size(); i++)
    {
      si(état == 0)
      {
        evt_en_cours = evts[i];
        état = 1;
      }
      sinon
      {
        // Evénements contigus
        si((abs((evts[i].t0 - evts[i-1].t1).nb_secondes()) < 1e-3)
           && (evts[i].label   == evts[i-1].label)
           && (evts[i].couleur == evts[i-1].couleur))
        {
          evt_en_cours.t1 = evts[i].t1;
        }
        sinon
        {
          // Nouvel événement
          evts2.push_back(evt_en_cours);
          evt_en_cours = evts[i];
        }
      }
    }
    si(état == 1)
      evts2.push_back(evt_en_cours);
    evts = evts2;
  }

  void FigureEvenements::rendre(Canva canva) const
  {
    entier nlignes = content.lignes.size();
    soit vertical = content.orientation_verticale;

    infos.lignes.resize(nlignes);
    pour(auto idl = 0; idl < nlignes; idl++)
    {
      infos.lignes[idl].evts.resize(content.lignes[idl].evts.size());
    }

    Axes axes;

    soit cfa = axes.get_config();
    cfa.axe_horizontal.afficher   = !vertical;
    cfa.axe_vertical.afficher     = vertical;
    cfa.grille_majeure.afficher   = non;
    cfa.grille_mineure.afficher   = non;

    cfa.axe_temporel.actif = oui;
    cfa.axe_temporel.t0    = content.t0;

    axes.configure(cfa);

    soit dim = canva.get_allocation();

    // Canva en pixels
    soit cp = canva.vue(Rect{0, 0, dim.l, dim.h});

    entier marge;// Marge pour les titres
    soit dt = (content.t1 - content.t0).nb_secondes();

    Canva c2, c3, cp2;

    si(vertical)
    {
      marge = 40;
      c2 = cp.clip(
          Rect{0, marge, dim.l, dim.h - marge},
          Rectf{0.0f, /*content.t0*/0.0f, 1.0f*nlignes, dt});
      c3 = axes.rendre(c2,
          Rectf{0, /*content.t0 + dt*/dt, 1.0f*nlignes, -dt});
      cp.set_align(Align::CENTRE, Align::DEBUT);
      pour(auto i = 0; i < nlignes; i++)
      {
        Pointf p = c3.v2c({i+0.5f, 0.0f});
        cp.texte(p.x, 0, content.lignes[i].nom, (dim.l / nlignes) - 2, marge);
        p = c3.v2c(Point{nlignes-i-1, 0});
        cp.ligne(p.x, 0, p.x, dim.h);
      }
    }
    sinon
    {
      marge = 200;

      Rect  zone_plot{marge, 0, dim.l - marge, dim.h};
      Rectf rdi_plot{/*content.t0*/0.0f, 0.0f, dt, 1.0f*nlignes};

      // Canva pixel intérieur
      cp2 = cp.clip(zone_plot, Rect{0, 0, dim.l - marge, dim.h});

      // Canva en unité utilisateur
      c2 = cp.clip(zone_plot, rdi_plot);

      c3 = axes.rendre(c2, //rdi_plot
           Rectf{/*content.t0*/0.0f, nlignes-1.0f, dt, -1.0f*nlignes});

      cp.set_align(Align::DEBUT, Align::CENTRE);
      pour(auto i = 0; i < nlignes; i++)
      {
        soit p = c3.v2c({0.0f, i-0.5f});
        cp.texte(0, p.y, content.lignes[i].nom, marge, ((dim.h / nlignes) * 3) / 4);
        p = c3.v2c(Pointf{0, i});//nlignes-i-1});
        cp.ligne(0, p.y, dim.l, p.y);
      }
    }


    pour(auto idl = 0; idl < nlignes; idl++)
    {
      soit &evts = content.lignes[idl].evts;

      c3.set_dim_fonte(0.5);
      c3.set_align(Align::CENTRE, Align::CENTRE);

      soit ide = 0;
      pour(auto &e: evts)
      {
        soit p0 = (e.t0-content.t0).nb_secondes(),
             p1 = (e.t1-content.t0).nb_secondes(),
             tmid = (p0 + p1) / 2.0f,
             ddt = p1 - p0;

        si(vertical)
        {
          c3.set_remplissage(oui, e.couleur);
          Rectf r{idl, p0, 1, ddt};
          c3.rectangle(r);
          c3.set_remplissage(non);

          si(!e.hide_label)
            c3.texte({idl+0.5f, tmid}, e.label, {1.0f, ddt});

          infos.lignes[idl].evts[ide].rdi = c2.v2c(c3.v2c(r));

        }
        sinon
        {
          Rectf r{p0, idl-1.0f, ddt, 1.0f};

          infos.lignes[idl].evts[ide].rdi = c3.v2c(r);
          infos.lignes[idl].evts[ide].rdi.x += marge;

          c3.set_remplissage(oui, e.sélectionné ? e.couleur.assombrir(0.7) : e.couleur);
          c3.set_epaisseur(e.sélectionné ? 3 : 1);
          c3.rectangle(r);
          c3.set_epaisseur(1);

          c3.set_remplissage(non);
          si(!e.hide_label)
          {
            c3.set_orientation(Orientation::VERTICALE);
            c3.texte({tmid, idl-0.5f}, e.label, {ddt, 1.0f});
            c3.set_orientation(Orientation::HORIZONTALE);
          }
        }
        ide++;
      }
    }
  }


}
