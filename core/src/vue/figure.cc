#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"
#include "tsd/fourier.hpp"
#include <deque>
#include <ranges>


using namespace std;

#define DBG(AA)


namespace tsd::vue {


/** @brief Système de coordonnées pour l'affichage d'une surface 2d. */
ParamGrille::ParamGrille(entier nc, entier nr, const Pointf &p0, const Pointf &p1)
{
  this->nc = nc;
  this->nr = nr;
  this->x0 = p0.x;
  this->y0 = p0.y;

  ẟx = (p1.x-x0)/(nc-1.0f);
  ẟy = (p1.y-y0)/(nr-1.0f);
}

Dimf ParamGrille::dim() const
{
  retourne {ẟx,ẟy};
}


Pointf ParamGrille::ctr(entier ix, entier iy) const
{
  retourne {x0+ix*ẟx, y0+iy*ẟy};
}
Pointf ParamGrille::tl(entier ix, entier iy) const
{
  retourne {x0+(ix-0.5f)*ẟx, y0+(iy-0.5f)*ẟy};
}
Pointf ParamGrille::br(entier ix, entier iy) const
{
  retourne {x0+(ix+0.5f)*ẟx, y0+(iy+0.5f)*ẟy};
}
Pointf ParamGrille::bl(entier ix, entier iy) const
{
  retourne {x0+(ix-0.5f)*ẟx, y0+(iy+0.5f)*ẟy};
}
Pointf ParamGrille::tr(entier ix, entier iy) const
{
  retourne {x0+(ix+0.5f)*ẟx, y0+(iy-0.5f)*ẟy};
}



static vector<string> parse_liste_chaines(cstring str, char separateur)
{
  vector<string> res;
  const char *s = str.c_str();
  string current;
  soit n = str.size();

  pour(auto i = 0u; i < n; i++)
  {
    si(s[i] != separateur)
    {
      char tmp[2];
      tmp[0] = s[i];
      tmp[1] = 0;
      current += string(tmp);
    }
    sinon
    {
      res.push_back(current);
      current = "";
    }
  }
  res.push_back(current);
  retourne res;
}


extern void stdo_attente_affichage_fait();

static bouléen mode_impression = non;

void set_mode_impression(bool mode_impression_)
{
  mode_impression = mode_impression_;
}

bouléen get_mode_impression()
{
  retourne mode_impression;
}

void Rendable::afficher(cstring titre, const Dim &dim) const
{
  stdo.affiche(shared_from_this(), titre, dim);
}

Image Rendable::genere_image(const Dim &dim_, const Couleur &arp) const
{
  Dim dim = dim_;
  si((dim.l <= 0) || (dim.h <= 0))
    dim = dimensions_idéales();
  si((dim.l <= 0) || (dim.h <= 0))
    dim = {1024, 1024};

  Canva c;
  c.set_allocation(dim);

  DBG(msg("Rendable/genere_image : dim={}, rendu -> canva...", dim);)

  rendre(c);

  DBG(msg("Rendable/genere_image : rendu -> image...");)

  retourne c.rendre(dim, Rect{0, 0, dim.l, dim.h}, arp);
}

/** @brief Enregistrement sous la forme d'un fichier image */
void Rendable::enregistrer(cstring chemin_fichier, const Dim &dim) const
{
  genere_image(dim).enregister(chemin_fichier);
}




struct Figure::Courbe::Impl
{
  Vecf x, y;

  Vecf ymin, ymax; // pour fill ymin ymax

  Tabf Z; // pour dessin surface

  Tabf couleurs_points_rvb;
  Vecf σ;

  Rectf rdi_z;

  string nom, format;

  Couleur couleur               = {0,0,0,180};

  entier epaisseur              = 3,
         dim_marqueur           = 5;

  bouléen remplissage           = non,
          remplissage_vers_ymin = oui,
          accu                  = non;

  float remplissage_vmin        = 0;

  Marqueur marqueur             = Marqueur::AUCUN;

  Trait trait                   = Trait::AUCUN;

  Impl()
  {
    si(mode_impression)
    {
      epaisseur = 5;
    }
  }
};



struct Figure::Impl: Rendable
{
  bouléen log_x = non, log_y = non,
          a_rdi_min = non, a_rdi = non;

  string nom, titre;

  // TODO : enlever ce mutable
  mutable Axes axes;

  Rectf rm_rdi;

  vector<Figure::Courbe> courbes;

  // Canva utilisateur (affiché au dessus de la figure)
  Canva canva_utilisateur;

  // Canva utilisateur (affiché avant la figure)
  Canva canva_utilisateur_pre;

  mutable float xmin = 0, xmax = 0, ymin = 0, ymax = 0;

  Dim dimensions_idéales() const
  {
    si(mode_impression)
      //retourne {1500, 400};
    retourne {2250, 600};
    sinon
      retourne {1500, 600};
  }

  Impl()
  {
    soit &cfg = axes.get_config();
    si(mode_impression)
    {
      cfg.legende.dim = 1.2;
      cfg.titre_echelle = 1.4;

      cfg.axe_horizontal.valeurs_echelle      = 1.5;
      cfg.axe_horizontal.valeurs_hauteur_max  = 40;
      cfg.axe_horizontal.valeurs_décalage     = 8;

      cfg.axe_vertical.valeurs_echelle        = 1.5;
      cfg.axe_vertical.valeurs_hauteur_max    = 40;
      cfg.axe_vertical.valeurs_décalage       = 8;
    }
    sinon
    {
      cfg.titre_echelle = 0.8;
      //cfg.axe_horizontal.valeurs_echelle      = 0.5;
      //cfg.axe_vertical.valeurs_echelle        = 1.5;
    }
    canva_utilisateur.set_allocation({1,1});
    canva_utilisateur_pre.set_allocation({1,1});
  }


  entier calc_minmax(const ConfigAxes &config, float &xmin, float &ymin, float &xmax, float &ymax) const
  {
    DBG(msg("calcul min/max ({} courbes)...", courbes.size()));

    pour(auto &c: courbes)
    {
      soit &ci = *(c.impl);

      si(ci.Z.rows() > 0)
      {
        soit ẟx = ci.rdi_z.l, ẟy = ci.rdi_z.h;

        xmin = ci.rdi_z.x - 0.5 * ẟx / (ci.Z.rows() - 1);
        ymin = ci.rdi_z.y - 0.5 * ẟy / (ci.Z.cols() - 1);
        xmax = ci.rdi_z.x + ẟx + 0.5 * ẟx / (ci.Z.rows() - 1);
        ymax = ci.rdi_z.y + ẟy + 0.5 * ẟy / (ci.Z.cols() - 1);
      }

      si(ci.x.rows() > 0)
      {
        soit ci_ymin = 0.0f, ci_ymax = 0.0f;
        soit first_min = oui, first_max = oui;

        pour(auto i = 0; i < ci.y.rows(); i++)
        {
          soit yi = ci.y(i);
          si(!isnan(yi) && !(isinf(yi)))
          {
            si(!config.axe_vertical.echelle_logarithmique || (yi > 0))
            {
              soit s = (ci.σ.rows() > 0) ? ci.σ(i) : 0.0f;
              soit y1 = yi + s, y2 = yi - s;

              si(config.axe_vertical.echelle_logarithmique)
              {
                // attention, y2 peut être négatif !!!
                y2 = pow(10.0f, 2 * log10(ci.y(i)) - log10(y1));
                si(isnan(y2) || isinf(y2))
                  y2 = 0;
              }

              si(!config.axe_vertical.echelle_logarithmique || (y2 > 0))
                si(first_min || (y2 < ci_ymin))
                {
                  ci_ymin = y2;
                  first_min = non;
                }
              si(!config.axe_vertical.echelle_logarithmique || (y1 > 0))
                si(first_max || (y1 > ci_ymax))
                {
                  ci_ymax = y1;
                  first_max = non;
                }
            }
          }
        }

        si((ci.trait == Trait::HISTO) || (ci.trait == Trait::BATON))
          ci_ymin = min(ci_ymin, 0.0f);

        DBG(msg("ci_ymin = {}, ci_ymax = {}", ci_ymin, ci_ymax));

        si(&c == &(courbes[0]))
        {
          xmin = ci.x.valeur_min();
          xmax = ci.x.valeur_max();
          ymin = ci_ymin;
          ymax = ci_ymax;
        }
        sinon
        {
          xmin = min(xmin, (float) ci.x.valeur_min());
          xmax = max(xmax, (float) ci.x.valeur_max());
          ymin = min(ymin, ci_ymin);
          ymax = max(ymax, ci_ymax);
        }


        DBG(msg("xmin={},ymin={}, xmax={},ymax={}", xmin, ymin, xmax, ymax));

        /*si(ci.sigma.rows() > 0)
        {
          ymin = min(ymin, (yn2 - ci.sigma).minCoeff());
          ymax = max(ymax, (yn2 + ci.sigma).maxCoeff());
        }*/
        DBG(msg("ci ymin = {}, ci ymax = {}, ymin = {}, ymax = {}", ci_ymin, ci_ymax, ymin, ymax);)
      }
      si(ci.trait == Trait::BATON)
      {
        ymin = min(ymin, ymax / 10.0f);
        ymax = max(ymax, ymin / 10.0f);
      }

    }

    si(courbes.size() == 0)
    {
      xmin = rm_rdi.x;
      xmax = rm_rdi.x + rm_rdi.l;
      ymin = rm_rdi.y;
      ymax = rm_rdi.y + rm_rdi.h;
    }

    si(a_rdi_min)
    {
      //infos("a rdi min : xmin = %f", s.rm_xmin);
      xmin = min(xmin, rm_rdi.x);
      xmax = max(xmax, rm_rdi.x + rm_rdi.l);
      ymin = min(ymin, rm_rdi.y);
      ymax = max(ymax, rm_rdi.y + rm_rdi.h);
    }
    sinon si(a_rdi)
    {
      xmin = rm_rdi.x;
      xmax = rm_rdi.x + rm_rdi.l;
      ymin = rm_rdi.y;
      ymax = rm_rdi.y + rm_rdi.h;
    }

    si(   isnan(xmin) || isnan(xmax) || isnan(ymin) || isnan(ymax)
       || isinf(xmin) || isinf(xmax) || isinf(ymin) || isinf(ymax))
    {
      msg_avert("xmin={},xmax={},ymin={},ymax={}",xmin,xmax,ymin,ymax);
      retourne -1;
    }

    si(xmax == xmin)
    {
      si(xmax == 0)
      {
        xmin -= 1;
        xmax += 1;
      }
      sinon
      {
        xmin = xmin - 0.1 * xmin;
        xmax = xmax + 0.1 * xmax;
      }
    }
    si(ymax == ymin)
    {
      si(ymax == 0)
      {
        ymin -= 1;
        ymax += 1;
      }
      sinon
      {
        ymin = ymin - 0.1 * abs(ymin);
        ymax = ymax + 0.1 * abs(ymax);
      }
    }
    //infos("xmin=%f,xmax=%f,ymin=%f,ymax=%f",xmin,xmax,ymin,ymax);
    retourne 0;
  }

  void extension_rdi(float &xmin, float &ymin, float &xmax, float &ymax) const
  {
    const auto &config = axes.get_config();
    soit dx = xmax - xmin, dy = ymax - ymin;
    si(dx <= 0)
      dx = 1;
    si(dy <= 0)
      dy = 1;

    soit delta_x = config.axe_horizontal.afficher ? 0.1f : 0.02f,
         delta_y = config.axe_vertical.afficher   ? 0.1f : 0.02f;

    si(config.axe_horizontal.echelle_logarithmique)
    {
      soit r = xmax / xmin;
      xmin = xmin / pow(r, 0.1);
      xmax = xmax * pow(r, 0.1);
    }
    sinon
    {
      xmin -= dx * delta_x;
      xmax += dx * delta_x;
    }


    si(config.axe_vertical.echelle_logarithmique)
    {
      // PB ICI !!!
      double r = ymax / ymin;
      double r2 = pow(r, 0.1);

      si(!isinf(r2))
      {
        DBG(msg("Extension RDI log (1) : ymin={},ymax={},r2={}", ymin, ymax, r2));

        ymin /= r2;

        //si(r2 > 10)
        //  r2 = 10;

        ymax *= r2;

        DBG(msg("Extension RDI log (2) : ymin={},ymax={}", ymin, ymax));
      }
    }
    sinon
    {
      ymin -= dy * delta_y;
      si(calc_titre().empty())
        ymax += dy * delta_y;
      sinon
        ymax += dy * (delta_y + 0.15);
    }
  }



  entier rendre00() const
  {
    soit config = axes.get_config();
    si(calc_minmax(config, xmin, ymin, xmax, ymax))
      retourne -1;

    si(!a_rdi)
      extension_rdi(xmin, ymin, xmax, ymax);
    retourne 0;
  }

  string calc_titre() const
  {
    si((courbes.size() == 1) && (titre.empty()))
      retourne courbes[0].impl->nom;
    retourne titre;
  }

  // TODO : gestion cas d'échec calc_minmax
  Canva rendre0(Canva canva_) const
  {
    soit config = axes.get_config();

    //config.legende.position = pos_cartouche;

    config.series.clear();

    si((!courbes.empty()) && (!courbes[0].impl->nom.empty()))
    {
      pour(auto &c: courbes)
      {
        ConfigSerie cs;
        cs.nom          = c.impl->nom;
        cs.couleur      = c.impl->couleur;
        cs.remplissage  = c.impl->remplissage;
        cs.trait        = c.impl->trait;
        cs.marqueur     = c.impl->marqueur;
        si(!cs.nom.empty())
          config.series.push_back(cs);
      }
    }

    config.legende.afficher = oui;
    si((courbes.size() == 1) &&
        ((config.titre.empty())
          || (config.titre == courbes[0].impl->nom)))
    {
      config.legende.afficher = non;
    }

    config.titre = calc_titre();

    //si(config.legende_dim)
    //config.legende_dim  = 6;
    axes.configure(config);

    DBG(msg("rendu axes, ymin={},ymax={},diff={}...", ymin, ymax, ymin-ymax));

    // Pb en mode logarithmique, si ymin = 1e-30, ymax = 2, on ne peut pas coder correctement ymin-ymax
    // Dans axes, le calcul suivant est fait :
    // ymax + (ymin-ymax) = 0....


    // TODO : création d'un PTI ici
    Canva canva = axes.rendre(canva_, xmin, xmax, ymin, ymax);
    DBG(msg("ok.");)

    retourne canva;
  }

  void rendre1(Canva canva_, Canva canva) const
  {
    DBG(msg("dessin courbes..."));

    canva_utilisateur_pre.forward(canva);

    pour(auto &c: courbes)
    {
      soit &ci = *(c.impl);
      soit &x = ci.x;
      soit &y = ci.y;

      soit npts = x.rows();
      assertion(npts == y.rows());


      DBG(msg("courbe : {} points.", npts));

      canva.set_couleur(ci.couleur);
      canva.set_epaisseur(ci.epaisseur);

      si(ci.Z.rows() > 0)
      {
        //soit nr = ci.Z.rows(), nc = ci.Z.cols();
        soit lst_opt = parse_liste_chaines(ci.format, ',');

        Tabf Zn;

        string str_cmap = "jet";

        si(lst_opt.size() >= 1)
          str_cmap = lst_opt[0];

        soit normaliser = oui;

        si((lst_opt.size() >= 2) && (lst_opt[1] == "o"))
          normaliser = non;

        si(normaliser)
        {
          soit [vmin, vmax] = ci.Z.valeurs_minmax();
          DBG(msg("plot image: avant normalisation : {} -- {}", vmin, vmax));
          Zn = (ci.Z - vmin) / (vmax - vmin);
          DBG(msg("plot image: après normalisation : {} -- {}", Zn.minCoeff(), Zn.valeur_max()));
        }
        sinon
          Zn = ci.Z;

        soit cmap = tsd::vue::cmap_parse(str_cmap);

#       if 0
        ParamGrille pg(nr, nc, ci.rdi_z.tl(), ci.rdi_z.br());

        DBG(msg("plot img : tl = {}, br = {}", ci.rdi_z.tl(), ci.rdi_z.br());)

        //canva.ligne(ci.rdi_z.tl(), ci.rdi_z.br());

        canva.active_contours(non);
        canva.active_remplissage(oui);

        pour(auto i = 0; i < nr; i++)
        {
          pour(auto j = 0; j < nc; j++)
          {
            soit coul = cmap->couleur(Zn(i,j));

            //si(((i % 5) == 0) && ((j % 5) == 0))
            ///{
            //coul = Couleur::rand();
            //coul.alpha = 128;
            canva.set_couleur_remplissage(coul);
            canva.rectangle(pg.tl(i, j), pg.br(i, j));
            //}
          }
        }
        canva.active_contours(oui);
        canva.active_remplissage(non);
#       endif

        canva.plot_cmap(Zn, ci.rdi_z, cmap);
      }

      si(c.impl->σ.rows() > 0)
      {
        si(c.impl->σ.rows() != npts)
        {
          msg_erreur("Affichage courbe avec σ : le nombre de points doit être identique pour σ et y.");
          continue;
        }
        soit cr = ci.couleur.eclaircir(0.3);
        canva.set_remplissage(oui, cr);
        canva.set_couleur(cr);
        pour(auto i = 0; i + 1 < npts; i++)
        {
          // A clarifier
          soit cvt = [&](float y, float s) -> tuple<float,float>
          {
            float y1 = y + s;
            // y2 = y^2 / (y + s)
            float y2 = pow(10.0f, 2 * log10(y) - log10(y1));
            si(isnan(y2) || isinf(y2))
              y2 = 0;
            retourne {y2, y1};
          };

          // attention, y2 peut être négatif !!!

          soit [y0m,y0p] = cvt(y(i), c.impl->σ(i));
          soit [y1m,y1p] = cvt(y(i+1), c.impl->σ(i+1));

          si((y0p != 0) && (y1p != 0))
            canva.remplissage_vertical(x(i), y0m, x(i+1),  y1m, y0p, y1p);
        }
        canva.set_remplissage(non, cr);
        canva.set_couleur(ci.couleur);
      }

      si(ci.trait == Trait::MINMAX)
      {
        soit cr = ci.couleur.eclaircir(0.3);
        canva.set_couleur(ci.couleur);
        canva.set_remplissage(oui, cr);
        pour(auto i = 0; i + 1 < npts; i++)
          canva.remplissage_vertical(x(i), ci.ymin(i), x(i+1), ci.ymin(i+1), ci.ymax(i), ci.ymax(i+1));
        canva.set_remplissage(non, cr);
      }
      sinon si((ci.trait == Trait::LIGNE)
              || (ci.trait == Trait::DOTTED))
      {

        //////////////////////////////////////////////////////////////////////////////////
        // 1er pré-traitement :
        // on regarde les points qui sont vraiment visibles, on peut ignorer les autres
        //
        // 2eme pré-traitement :
        // si il reste beaucoup de points, on fait une décimation spéciale
        // (calcul min / max sur des intervalle donnés)

        Vecf xr1, yr1;

        soit a_couleurs = ci.couleurs_points_rvb.rows() > 0;

        si((npts > 2000) && !a_couleurs)
        {
          bouléen monotone = oui;
          entier idmin = -1, idmax = -1;

          pour(auto i = 0; i + 1 < npts; i++)
          {
            // Vérification suite monotone, sinon on fait pas ça !
            si(x(i) > x(i+1))
            {
              monotone = non;
              break;
            }
            si((x(i) >= xmin) && (idmin == -1))
              idmin = i;
            si((x(i) > xmax) && (idmax == -1))
              idmax = i;
          }
          si((idmin >= 0) && (idmax >= 0))
          {
            entier npts3 = idmax - idmin + 1;
            // si monotone, et si la réduction en vaut la peine
            si(monotone && (npts3 < 0.5 * npts)) // Un peu arbitraire
            {
              //msg("monotone, npts={}, npts3 = {}, idmin = {}, idmax={}, xmin={}, xmax={}", npts, npts3, idmin, idmax, xmin, xmax);
              xr1 = x.segment(idmin, npts3);
              yr1 = y.segment(idmin, npts3);
              x = xr1;
              y = yr1;
              npts = npts3;
            }
          }
        }


        //msg("plot : npts={}", npts);
        Vecf xr, yrmin, yrmax, yrmoy;
        bouléen dessine_sigma = non;
        si((npts >= 2000) && !a_couleurs)
        {
          dessine_sigma = oui;
          entier R = floor(npts / 1000);
          entier npts2 = npts / R;

          //msg("plot : npts={} -> R={}, npts2={}", npts, R, npts2);

          xr.resize(npts2);
          yrmin.resize(npts2);
          yrmax.resize(npts2);
          yrmoy.resize(npts2);

          pour(auto i = 0; i < npts2; i++)
          {
            xr(i) = x(i * R);
            yrmin(i) = yrmax(i) = y(i * R);
            yrmoy(i) = y.segment(i * R, R).moyenne();
            pour(auto j = 1; j < R; j++)
            {
              yrmax(i) = max(yrmax(i), y(i*R+j));
              yrmin(i) = min(yrmin(i), y(i*R+j));
            }
          }
          npts = npts2;
          x = xr;
          y = yrmoy;

          // Fait après remplissage éventuel
          /*Couleur cr1 = axes.get_config().est_mat_sur_clair() ? ci.couleur.eclaircir(0.7) : ci.couleur.assombrir(0.7);
          canva.set_couleur(cr1);
          pour(auto i = 0; i + 1 < npts; i++)
          {
            si((x(i+1) <= xmax) && (x(i) >= xmin))
              canva.remplissage_vertical(x(i), yrmax(i), x(i+1),  yrmax(i+1), yrmin(i), yrmin(i+1));
          }
          canva.set_remplissage(non);*/
        }


        si(ci.remplissage)
        {
          Couleur cr1 = axes.get_config().est_mat_sur_clair() ?
              ci.couleur.eclaircir(0.05) : ci.couleur.assombrir(0.05);
          canva.set_couleur(cr1);
          //auto cr = ci.couleur.eclaircir(0.3);

          soit vmin = ci.remplissage_vers_ymin ?
                      ymin + (ymax-ymin)/200 : ci.remplissage_vmin;

          pour(auto i = 0; i + 1 < npts; i++)
          {
            si((x(i+1) <= xmax) && (x(i) >= xmin) && (y(i) >= ymin) && (y(i+1) >= ymin))
            {
              canva.remplissage_vertical(x(i), y(i), x(i+1),  y(i+1), vmin, vmin);
            }
          }
          canva.set_remplissage(non);
          canva.set_couleur(ci.couleur);
        }

        si(dessine_sigma)
        {
          Couleur cr1 = axes.get_config().est_mat_sur_clair() ?
              ci.couleur.eclaircir(0.2) : ci.couleur.assombrir(0.2);
          canva.set_couleur(cr1);
          pour(auto i = 0; i + 1 < npts; i++)
          {
            si((x(i+1) <= xmax) && (x(i) >= xmin))
              canva.remplissage_vertical(x(i), yrmax(i), x(i+1),  yrmax(i+1), yrmin(i), yrmin(i+1));
          }
          canva.set_remplissage(non);
        }

        bouléen dotted = ci.trait == Trait::DOTTED;
        canva.set_dotted(dotted);
        canva.set_couleur(ci.couleur);
        pour(auto i = 0; i + 1 < npts; i++)
        {
          si((x(i+1) <= xmax) && (x(i) >= xmin) && (y(i) >= ymin) && (y(i+1) >= ymin)
              && (!isnan(y(i)))
              && (!isnan(y(i+1))))
          {
            si((ci.couleurs_points_rvb.rows() > 0) && (i < ci.couleurs_points_rvb.cols()))
            {
              soit couleur = Couleur{
              ci.couleurs_points_rvb(2, i),
              ci.couleurs_points_rvb(1, i),
              ci.couleurs_points_rvb(0, i)};

              canva.set_couleur(couleur);
            }


            canva.ligne(x(i), y(i), x(i+1), y(i+1));
          }
        }
        canva.set_dotted(non);
      }
      sinon si(c.impl->trait == Trait::BATON)
      {
        pour(auto i = 0; i < npts; i++)
        {
          si(!isinf(y(i)))
          {
            canva.ligne(x(i), y(i), x(i), 0.0f);
            //canva.marqueur({x(i), y(i)}, Marqueur::CERCLE, 5);
          }
        }
      }
      sinon si(c.impl->trait == Trait::HISTO)
      {
        si(x.dim() > 1)
        {
          canva.set_remplissage(oui, ci.couleur);
          soit dx = x(1) - ci.x(0);
          pour(auto i = 0; i < npts; i++)
          {
            si(!isinf(y(i)))
              canva.rectangle(x(i), y(i), x(i) + dx, 0);
          }
          canva.set_remplissage(non, ci.couleur);
        }
      }




      // Comment faire pour accumuler ?
      //  Idée : adapter par rapport à la densité maximale
      //         -> il faudrait pouvoir dessiner en flottant et normaliser
      // Ou alors, on ajoute 1 unité à chaque itération
      //  -> puis on normalise

      si(ci.accu)
      {
        //msg("figure : rendu accu ({} points)", npts);
        Veccf pts(npts);
        canva.set_couleur(ci.couleur);
        pour(auto i = 0; i < npts; i++)
        {
          pts(i).real(ci.x(i));
          pts(i).imag(ci.y(i));
        }
        canva.dessine_accu(pts);
      }

      si(ci.marqueur == Marqueur::AUCUN)
        continue;

      pour(auto i = 0; i < npts; i++)
      {
        si(ci.accu)
          continue;

        Couleur couleur = ci.couleur;

        si((ci.couleurs_points_rvb.rows() > 0) && (i < ci.couleurs_points_rvb.cols()))
          couleur = Couleur{
          ci.couleurs_points_rvb(2, i),
          ci.couleurs_points_rvb(1, i),
          ci.couleurs_points_rvb(0, i)};

        canva.set_couleur(couleur);

        soit x = ci.x(i), y = ci.y(i);
        canva.set_remplissage(oui, couleur);
        canva.marqueur({x, y}, ci.marqueur, ci.dim_marqueur);
        canva.set_remplissage(non);
      }
    }
    axes.post_rendu(canva_);
    canva_utilisateur.forward(canva);
  }

  void rendre(Canva canva_) const
  {
    // Détection ici ????
    si(rendre00())
      retourne;
    soit canva = rendre0(canva_);
    rendre1(canva_, canva);
  }

  Figure::Courbe plot(const Vecf &x, const Vecf &y, cstring format,
      cstring titre);
};




void Figure::Courbe::def_remplissage(bouléen actif, bouléen vers_vmin, float vmin)
{
  impl->remplissage           = actif;
  impl->remplissage_vers_ymin = !vers_vmin;
  impl->remplissage_vmin      = vmin;
}

void Figure::Courbe::def_epaisseur(entier ep)
{
  impl->epaisseur = ep;
}

void Figure::Courbe::def_couleur(const tsd::vue::Couleur &coul)
{
  impl->couleur = coul;
}

void Figure::Courbe::def_σ(const Vecf &σ)
{
  impl->σ = σ;
}

void Figure::Courbe::def_légende(cstring titre)
{
  impl->nom = titre;
}

void Figure::Courbe::def_dim_marqueur(entier dim)
{
  impl->dim_marqueur = dim;
}

// Définit la couleur de chacun des points de la courbe
void Figure::Courbe::def_couleurs(const Vecf &c, cstring cmap_nom)
{
  Tabf rvb(3, c.rows());

  soit cmap = tsd::vue::cmap_parse(cmap_nom);

  pour(auto i = 0; i < c.rows(); i++)
  {
    float R = 0, V = 0, B = 0;
    cmap->calc(c(i), R, V, B);
    rvb(0, i) = R;
    rvb(1, i) = V;
    rvb(2, i) = B;
  }
  impl->couleurs_points_rvb = 255 * rvb;
}

Figure::Figure(cstring nom)
{
  impl = make_shared<Impl>();
  impl->nom = nom;
}

Axes Figure::axes()
{
  retourne impl->axes;
}

void Figure::clear()
{
  impl->courbes.clear();
  impl->canva_utilisateur.clear();
  impl->canva_utilisateur_pre.clear();
  impl->titre = "";
  impl->nom = "";
  soit &c = impl->axes.get_config();
  c.axe_horizontal.label  = "";
  c.axe_vertical.label    = "";
}

struct Figures::Impl: Rendable
{
  entier n = -1, m = -1, pos = 1;
  deque<Figure> subplots;

  tuple<entier, entier> get_nm() const;

  Dim dimensions_idéales() const
  {
    soit [nrows, ncols] = get_nm();

    si(mode_impression)
      retourne {ncols * 1500, nrows * 400};
    sinon
      retourne {ncols * 1500, nrows * 600};
  }

  Impl(entier n = 1, entier m = 1)
  {
    this->n = n;
    this->m = m;
    pos = 1;
  }

  Figure *gcs()
  {
    si(subplots.empty())
      subplots.push_back(Figure());

    si((pos-1) >= (entier) subplots.size())
    {
      msg_erreur("gcs() : pos = {}, subplots.size() = {}", pos, subplots.size());
      retourne &(subplots[0]);
    }
    retourne &(subplots[pos - 1]);
  }

  void rendre(Canva canva) const
  {
    entier nsubs = subplots.size();

    soit [nrows, ncols] = get_nm();

    si(nrows * ncols == 0)
      retourne;

    assertion((ncols >= 1) && (nrows >= 1));
    assertion_msg(nsubs <= (ncols * nrows), "nsubs = {}, m = {}, n = {}", nsubs, ncols, nrows);

    soit rdi = canva.get_rdi();
    entier mx = rdi.l / ncols, my = rdi.h / nrows;

    DBG(msg("parcours subplots..."));
    entier i = 0;
    pour(auto &s: subplots)
    {
      soit col = i % ncols, row = i / ncols;
      soit px = col * mx, py = row * my;
      // Ajoute une marge verticale
      soit hauteur = my - 15;
      si(hauteur > 0)
        s.rendre(canva.clip(Rect{px, py, mx, hauteur}, Rect{0, 0, mx, hauteur}));
      i++;
    }
    DBG(msg("ok.");)
  }


};

sptr<const Rendable> Figures::rendable() const
{
  retourne impl;
}


void Figures::afficher(cstring titre, const Dim &dim) const
{
  ARendable::afficher(titre, dim);
}

Figures::Figures(entier n, entier m)
{
  impl = make_shared<Impl>(n, m);
}

void Figures::clear()
{
  *impl = Impl();
}

Figure Figures::gf(entier sel)
{
  si(sel >= (entier) impl->subplots.size())
    échec("Figures::gcf({}) : index invalide (n figures = {})", sel, impl->subplots.size());
  retourne impl->subplots[sel];
}

Figure Figures::subplot(entier n, entier m, entier pos_)
{
  soit pos = pos_;
  si(((impl->n != n) || (impl->m != m)) && (impl->subplots.size() > 0))
  {
    impl->subplots.clear();
  }
  impl->n   = n;
  impl->m   = m;

  si(pos <= 0)
  {
    impl->subplots.push_back({});
    pos = impl->subplots.size();
  }

  si(pos > n * m)
  {
    msg_erreur("subplot : n = {}, m = {}, pos = {} (pos prm = {}), nbsubplots={}.", n, m, pos, pos_, impl->subplots.size());
  }

  tantque(pos > (entier) impl->subplots.size())
    impl->subplots.push_back({});

  impl->pos = pos;
  retourne impl->subplots[pos-1];
}

Figure Figures::gcf()
{
  retourne impl->subplots.back();
}

vector<Figure> Figures::subplots(int n)
{
  vector<Figure> l;
  pour(auto i = 0; i < n; i++)
    l.push_back(subplot());
  retourne l;
}

Figure Figures::subplot(entier i)
{
  si(i < 0)
  {
    impl->n = impl->m  = -1;
    impl->subplots.push_back({});
    impl->pos = impl->subplots.size();
    retourne impl->subplots.back();
  }

  soit n   = i / 100,
       m   = (i - 100 * n) / 10,
       pos = i % 10;

  retourne subplot(n, m, pos);
}


vector<Figure::Courbe> &Figure::courbes()
{
  retourne impl->courbes;
}

void Figure::def_rdi_min(const Rectf &rdi)
{
  assertion(impl);
  impl->a_rdi_min = oui;
  impl->rm_rdi    = rdi;
}

// TODO : ordre pas cohérent
void Figure::def_rdi(const Rectf &rdi)
{
  impl->a_rdi   = oui;
  impl->rm_rdi  = rdi;
}

Rectf Figure::get_rdi() const
{
  retourne impl->rm_rdi;
}

Figure::Courbe Figure::plot_minmax(const Vecf &x, const Vecf &y1, const Vecf &y2)
{
  Courbe res;
  res.impl = make_shared<Courbe::Impl>();
  res.impl->x     = x;
  res.impl->ymin  = y1;
  res.impl->ymax  = y2;
  res.impl->y     = (y1 + y2) / 2;
  res.impl->trait = Trait::MINMAX;

  impl->courbes.push_back(res);
  retourne res;
}

Figure::Courbe Figure::plot_iq_int(const Veccf &z, cstring format, cstring titre)
{
  axes().set_isoview(oui);
  soit x = real(z), y = imag(z);
  retourne impl->plot(x, y, format, titre);
}


Figure::Courbe Figure::plot_img(const Rectf &rdi, const Tabf &Z, cstring format)
{
  Courbe c;
  c.impl = make_shared<Courbe::Impl>();
  c.impl->Z       = Z;
  c.impl->format  = format;
  c.impl->rdi_z   = rdi;
  impl->courbes.push_back(c);

  retourne c;
}

Figure::Courbe Figure::plot_img(const Tabf &Z, cstring format)
{
  retourne plot_img({0, 0, Z.rows() - 1, Z.cols() - 1}, Z, format);
}



Figure::Courbe Figure::plot(const float &x, const float &y, cstring format)
{
  retourne plot(Vecf::valeurs({x}), Vecf::valeurs({y}), format);
}

Figure::Courbe Figure::plot_int(const Vecf &y, cstring format, cstring titre)
{
  retourne impl->plot(Vecf{}, y, format, titre);
}




Figure::Courbe Figure::plot_psd_int(const Veccf &y, float fe, cstring format, cstring titre)
{
  soit n = y.rows();

  si(n == 0)
    retourne {};

  soit max_coef_imag = abs(imag(y)).valeur_max();
  soit est_reel = (max_coef_imag == 0);
  soit fen = tsd::filtrage::fenetre("hn", n, non);
  soit ycf = y * fen.as_complex();

  soit Y = pow2db(abs2(tsd::fourier::fft(ycf)));

  Vecf freq;

  si(est_reel)
  {
    Y = Y.head(n / 2);
    freq = linspace(0, fe/2, Y.rows());
  }
  sinon
  {
    Y = tsd::fourier::fftshift(Y);
    freq = linspace(-fe/2, fe/2, Y.rows());
  }

  impl->axes.get_config().axe_horizontal.label = "Fréquence";

  retourne plot(freq, Y, format, titre);
}

Figure::Courbe Figure::plot_int(const Vecf &x, const Vecf &y, cstring format, cstring titre)
{
  retourne impl->plot(x, y, format, titre);
}

Figure::Courbe Figure::Impl::plot(const Vecf &x, const Vecf &y_, cstring format_, cstring titre)
{
  Figure::Courbe res;
  res.impl = make_shared<Figure::Courbe::Impl>();
  res.impl->nom = titre;

  string format = format_;

  si(format.empty())
  {
    entier nc = courbes.size();
    static const vector<string> fdef = {"b-", "g-", "r-", "y-", "c-", "m-", "a-", "f-",
    "B-", "V-", "R-", "J-", "C-", "W-", "O-", "M-", "k-"};
    format = fdef[nc % 17];
  }

  soit n = y_.rows();
  si(n > 1e6)
    msg_avert("plot : {} échantillons.", n);

  soit y = y_.clone();

  si(n == 0)
  {
    DBG(msg_avert("plot : y.rows() == 0. Titre = [{} / {}]", nom, res.impl->nom);)
    retourne res;
  }

  soit s = this;

  Vecf X;
  si(x.rows() == 0)
  {
    X = linspace(0, n - 1, n);
    DBG(msg("plot: x non spécifié -> X = {}", X));
  }
  sinon
    X = x.clone();

  si(X.rows() != n)
  {
    msg_erreur("Figure::plot(x,y,...) : x.rows() = {}, y.rows() = {}, y_.rows() = {}", X.rows(), y.rows(), y_.rows());
    retourne res;
  }
  assertion(X.rows() == y.rows());

  res.impl->x       = X;
  res.impl->y       = y;
  res.impl->format  = format;

  soit present = [](cstring s, char c)
  {
    retourne s.find(c) != string::npos;
  };


  struct CodeCouleur {char code; Couleur couleur;};
  static const CodeCouleur codes[] =
  {
      {'b', Couleur::BleuSombre},
      {'g', Couleur::VertSombre},
      {'r', Couleur::RougeSombre},
      {'y', Couleur::JauneSombre},
      {'c', Couleur::CyanSombre},
      {'m', Couleur::VioletSombre},
      {'a', Couleur::OrangeSombre},
      {'f', Couleur::MarronSombre},


      {'B', Couleur::Bleu},
      {'V', Couleur::Vert},
      {'R', Couleur::Rouge},
      {'J', Couleur::Jaune},
      {'C', Couleur::Cyan},
      {'W', Couleur::Violet},
      {'O', Couleur::Orange},
      {'M', Couleur::Marron},



      {'k', Couleur::Noir},
      {'w', Couleur::Blanc}

  };

  entier nc = courbes.size();
  res.impl->couleur = codes[nc % 16].couleur;

  pour(auto c: codes)
    si(present(format, c.code))
      res.impl->couleur = c.couleur;

  res.impl->couleur.alpha = 128;

  res.impl->accu = present(format, 'A');

  struct CodeTrait {char code; Trait trait;};

  // https://help.scilab.org/docs/6.0.2/en_US/LineSpec.html
  CodeTrait codest[] =
  {
      {'-', Trait::LIGNE},
      {'|', Trait::BATON},
      {'h', Trait::HISTO},
      {':', Trait::DOTTED}
  };

  pour(auto c: codest)
    si(present(format, c.code))
      res.impl->trait = c.trait;


  struct CodeMarqueur {char code; Marqueur marqueur;};
  static const CodeMarqueur codesm[] =
  {
      {'o', Marqueur::CERCLE},
      {'*', Marqueur::ETOILE},
      {'.', Marqueur::POINT},
      {'+', Marqueur::CROIX},
      {'s', Marqueur::CARRE},
      {'d', Marqueur::DIAMANT},
  };

  pour(auto c: codesm)
    si(present(format, c.code))
      res.impl->marqueur = c.marqueur;

  s->courbes.push_back(res);
  retourne res;
}

void Figure::def_pos_legende(cstring code)
{
  axes().get_config().legende.position = code;
}

void Figure::titre(cstring titre_global)
{
  impl->titre = titre_global;
}

void Figure::titres(cstring titre_global,
                    cstring axe_x,
                    cstring axe_y)
{
  soit &c = impl->axes.get_config();
  c.axe_horizontal.label  = axe_x;
  c.axe_vertical.label    = axe_y;
  impl->titre = titre_global;
}


void Figure::attente_ihm()
{
  stdo.fin();
}

tuple<entier, entier> Figures::Impl::get_nm() const
{
  entier nsubs = subplots.size();

  si(nsubs == 0)
    retourne {0, 0};

  si(m == -1)
  {
    si(nsubs == 1)
      retourne {1, 1};
    sinon si(nsubs == 2)
      retourne {2, 1};
    sinon si(nsubs == 3)
      retourne {3, 1};
    sinon
    {
      soit m = (entier) floor(sqrt(nsubs));
      soit n = (nsubs + m - 1) / m;
      retourne {n, m};
    }
  }
  retourne {m, n};
}



void Figure::def_nom(cstring s)
{
  impl->nom = s;
}

string Figure::lis_nom() const
{
  retourne impl->nom;
}

Figure Figure::clone() const
{
  Figure res;
  res.impl = shared_ptr<Impl>(new Figure::Impl(*impl));
  res.impl->axes = impl->axes.clone();
  retourne res;
}

Canva Figure::canva_pixel(const Dim &allocation)
{
  soit sz = axes().get_dim_graphique(allocation);
  soit a = canva();
  Rectf r = get_rdi();
  retourne impl->canva_utilisateur.clip(r, Rect{0, sz.h, sz.l, -sz.h});
}

Canva Figure::canva()
{
  retourne impl->canva_utilisateur;
}

Canva Figure::canva_pre()
{
  retourne impl->canva_utilisateur_pre;
}

void Figure::def_echelle(bouléen log_x, bouléen log_y)
{
  impl->log_x = log_x;
  impl->log_y = log_y;
}


sptr<const Rendable> Figure::rendable() const
{
  retourne impl;
}

void Figure::rendre1(Canva canva_, Canva canva) const
{
  retourne impl->rendre1(canva_, canva);
}

Canva Figure::rendre0(Canva canva_) const
{
  retourne impl->rendre0(canva_);
}

entier Figure::rendre00() const
{
  retourne impl->rendre00();
}


}

