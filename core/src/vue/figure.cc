#include "tsd/tsd.hpp"
#include "tsd/figure.hpp"
#include "tsd/fourier.hpp"
#include <cstdarg>
#include <deque>


#define DBG(AA)

namespace tsd::vue {

extern void stdo_attente_affichage_fait();

static bool mode_impression = false;

void set_mode_impression()
{
  mode_impression = true;
}

void Rendable::afficher(const std::string &titre, const Dim &dim) const
{
  stdo.affiche(shared_from_this(), titre, dim);
}

Image Rendable::genere_image(const Dim &dim_, const Couleur &arp) const
{
  Dim dim = dim_;
  if((dim.l <= 0) || (dim.h <= 0))
  {
    if(mode_impression)
      dim = {1500, 400};
    else
      dim = {1500, 600};
  }
  Canva c;
  c.set_allocation(dim);

  DBG(msg("Rendable/genere_image : dim={}, rendu -> canva...", dim);)

  rendre(c);

  //auto rdi = c.calc_rdi_englobante();

  //msg("Rendable/genere_image : dim={}, rdi_englobante={}.", dim, rdi);

  DBG(msg("Rendable/genere_image : rendu -> image...");)

  return c.rendre(dim, Rect{0, 0, dim.l, dim.h}, arp);
  //return c.rendre(dim, rdi, arp);
}

/** @brief Enregistrement sous la forme d'un fichier image */
void Rendable::enregistrer(const std::string &chemin_fichier, const Dim &dim) const
{
  genere_image(dim).enregister(chemin_fichier);
}


void PlotEvt::maj(ArrayXi &evt, const EvtPlotConfig &cfg)
{
  this->evt = evt;
  this->cfg = cfg;
}

void PlotEvt::rendre(Canva canva) const
{
  int nevt = cfg.evts.size();
  int nt = evt.rows();

  if(nt == 0)
    return;

  Axes axes;

  auto cfa = axes.get_config();
  cfa.axe_vertical.afficher = false;
  cfa.grille_majeure.afficher = false;
  cfa.grille_mineure.afficher = false;
  axes.configure(cfa);

  auto dim = canva.get_allocation();

  DBG(msg("plot EVT2 : nevt={}, nt={}, allocation={}", nevt, nt, dim);)

  // On garde 100 pixels à gauche pour les titres
  // Canva en pixels
  Canva cp = canva.vue(Rect{0, 0, dim.l, dim.h});

  int marge = 200;
  Canva c2 = cp.clip(Rect{marge, 0, dim.l - marge, dim.h}, Rect{0, 0, nt, nevt});

  // TODO : problématique
  auto c3 = axes.rendre(c2, Rect{0, nevt-1, nt, -nevt});

  cp.set_align(Align::DEBUT, Align::CENTRE);
  for(auto i = 0; i < nevt; i++)
  {
    if(i != cfg.idle_index)
    {
      Pointf p = c3.v2c({0.0f, i-0.5f});//nevt-i-1+0.5f});
      cp.texte(0, p.y, cfg.evts[i].nom, marge, (dim.h / nevt) - 2);
    }
    Pointf p = c3.v2c(Point{0, nevt-i-1});
    cp.ligne(0, p.y, dim.l, p.y);
  }

  for(auto i = 0; i < nt; i++)
  {
    int e = evt(i);
    //msg("e = {}", e);
    if(e != cfg.idle_index)
    {
      auto coul = cfg.evts[e].coul;
      c3.set_remplissage(true, coul);
      c3.rectangle(i, e-1, i+1, e);
      c3.set_remplissage(false);
    }
  }
}



struct Figure::Courbe::Impl
{
  ArrayXf x, y;

  ArrayXf ymin, ymax; // pour fill ymin ymax

  ArrayXXf Z; // pour dessin surface
  Rectf rdi_z;

  std::string nom, format;
  Couleur couleur = {0,0,0,180};
  int epaisseur = 1;
  bool remplissage = false;
  bool remplissage_vers_ymin = true;
  float remplissage_vmin = 0;
  bool accu = false;

  ArrayXXf couleurs_points_rvb;
  ArrayXf σ;

  Marqueur marqueur = Marqueur::AUCUN;
  int dim_marqueur = 5;
  Trait trait = Trait::AUCUN;
};



struct Figure::Impl: Rendable
{
  bool log_x = false, log_y = false;
  std::string nom;

  std::string titre;

  // TODO : enlever ce mutable
  mutable Axes axes;
  bool a_rdi_min = false, a_rdi = false;
  Rectf rm_rdi;
  std::vector<Figure::Courbe> courbes;

  // Canva utilisateur
  Canva canva_utilisateur;

  std::string pos_cartouche = "ne";

  Impl()
  {
    if(mode_impression)
      axes.get_config().legende.dim = 0.8;
    canva_utilisateur.set_allocation({1,1});
  }


  int calc_minmax(const ConfigAxes &config, float &xmin, float &ymin, float &xmax, float &ymax) const
  {
    DBG(msg("calcul min/max ({} courbes)...", courbes.size()));

    for(auto &c: courbes)
    {
      auto &ci = *(c.impl);

      if(ci.Z.rows() > 0)
      {
        xmin = ci.rdi_z.x;
        ymin = ci.rdi_z.y;
        xmax = ci.rdi_z.x + ci.rdi_z.l - 1;
        ymax = ci.rdi_z.y + ci.rdi_z.h - 1;
      }

      if(ci.x.rows() > 0)
      {
        //ArrayXf yn  = ci.y.isNaN().select(0,ci.y);
        //ArrayXf yn2 = yn.isInf().select(0,yn);
        //for(auto)
        float ci_ymin = ci.y(0), ci_ymax = ci.y(0);

        if(std::isinf(ci_ymin))
          ci_ymin = 0;

        if(std::isinf(ci_ymax))
          ci_ymax = 0;

        bool first = true;

        for(auto i = 0; i < ci.y.rows(); i++)
        {
          if(!std::isnan(ci.y(i)) && !(std::isinf(ci.y(i))))
          {
            if(!config.axe_vertical.echelle_logarithmique || (ci.y(i) > 0))
            {
              float s = (ci.σ.rows() > 0) ? ci.σ(i) : 0;
              float y1 = ci.y(i) + s;
              float y2 = ci.y(i) - s;

              if(config.axe_vertical.echelle_logarithmique)
              {
                // attention, y2 peut être négatif !!!
                y2 = std::pow(10.0f, 2 * log10(ci.y(i)) - log10(y1));
                if(std::isnan(y2) || std::isinf(y2))
                  y2 = 0;
              }

              if(!config.axe_vertical.echelle_logarithmique || (y2 > 0))
                if(first || (y2 < ci_ymin))
                  ci_ymin = y2;
              if(!config.axe_vertical.echelle_logarithmique || (y1 > 0))
                if(first || (y1 > ci_ymax))
                  ci_ymax = y1;
              first = false;
            }
          }
        }

        if((ci.trait == Trait::HISTO) || (ci.trait == Trait::BATON))
        {
          ci_ymin = std::min(ci_ymin, 0.0f);
        }

        //auto ci_ymin = yn2.minCoeff();
        //auto ci_ymax = yn2.maxCoeff();

        //infos("ci_ymin = %f, ci_ymax = %f", ci_ymin, ci_ymax);

        DBG(msg("ci_ymin = {}, ci_ymax = {}", ci_ymin, ci_ymax));

        if(&c == &(courbes[0]))
        {
          xmin = ci.x.minCoeff();
          xmax = ci.x.maxCoeff();
          ymin = ci_ymin;
          ymax = ci_ymax;
        }
        else
        {
          xmin = std::min(xmin, ci.x.minCoeff());
          xmax = std::max(xmax, ci.x.maxCoeff());
          ymin = std::min(ymin, ci_ymin);
          ymax = std::max(ymax, ci_ymax);
        }


        DBG(msg("xmin={},ymin={}, xmax={},ymax={}", xmin, ymin, xmax, ymax));

        /*if(ci.sigma.rows() > 0)
        {
          ymin = std::min(ymin, (yn2 - ci.sigma).minCoeff());
          ymax = std::max(ymax, (yn2 + ci.sigma).maxCoeff());
        }*/
        DBG(msg("ci ymin = {}, ci ymax = {}, ymin = {}, ymax = {}", ci_ymin, ci_ymax, ymin, ymax);)
      }
      if(ci.trait == Trait::BATON)
      {
        ymin = std::min(ymin, ymax / 10.0f);
        ymax = std::max(ymax, ymin / 10.0f);
      }

    }

    if(courbes.size() == 0)
    {
      xmin = rm_rdi.x;
      xmax = rm_rdi.x + rm_rdi.l;
      ymin = rm_rdi.y;
      ymax = rm_rdi.y + rm_rdi.h;
    }

    if(a_rdi_min)
    {
      //infos("a rdi min : xmin = %f", s.rm_xmin);
      xmin = std::min(xmin, rm_rdi.x);
      xmax = std::max(xmax, rm_rdi.x + rm_rdi.l);
      ymin = std::min(ymin, rm_rdi.y);
      ymax = std::max(ymax, rm_rdi.y + rm_rdi.h);
    }
    else if(a_rdi)
    {
      xmin = rm_rdi.x;
      xmax = rm_rdi.x + rm_rdi.l;
      ymin = rm_rdi.y;
      ymax = rm_rdi.y + rm_rdi.h;
    }

    if(   std::isnan(xmin) || std::isnan(xmax) || std::isnan(ymin) || std::isnan(ymax)
       || std::isinf(xmin) || std::isinf(xmax) || std::isinf(ymin) || std::isinf(ymax))
    {
      msg_avert("xmin={},xmax={},ymin={},ymax={}",xmin,xmax,ymin,ymax);
      return -1;
    }

    if(xmax == xmin)
    {
      if(xmax == 0)
      {
        xmin -= 1;
        xmax += 1;
      }
      else
      {
        xmin = xmin - 0.1 * xmin;
        xmax = xmax + 0.1 * xmax;
      }
    }
    if(ymax == ymin)
    {
      if(ymax == 0)
      {
        ymin -= 1;
        ymax += 1;
      }
      else
      {
        ymin = ymin - 0.1 * std::abs(ymin);
        ymax = ymax + 0.1 * std::abs(ymax);
      }
    }
    //infos("xmin=%f,xmax=%f,ymin=%f,ymax=%f",xmin,xmax,ymin,ymax);
    return 0;
  }

  void extension_rdi(float &xmin, float &ymin, float &xmax, float &ymax) const
  {
    const auto &config = axes.get_config();
    float dx = xmax - xmin, dy = ymax - ymin;
    if(dx <= 0)
      dx = 1;
    if(dy <= 0)
      dy = 1;

    float delta_x = config.axe_horizontal.afficher ? 0.1 : 0.02;
    float delta_y = config.axe_vertical.afficher ? 0.1 : 0.02;

    if(config.axe_horizontal.echelle_logarithmique)
    {
      float r = xmax / xmin;
      xmin = xmin / std::pow(r, 0.1);
      xmax = xmax * std::pow(r, 0.1);
    }
    else
    {
      xmin -= dx * delta_x;
      xmax += dx * delta_x;
    }


    if(config.axe_vertical.echelle_logarithmique)
    {
      // PB ICI !!!
      double r = ymax / ymin;
      double r2 = std::pow(r, 0.1);

      if(!std::isinf(r2))
      {
        DBG(msg("Extension RDI log (1) : ymin={},ymax={},r2={}", ymin, ymax, r2));

        ymin /= r2;

        //if(r2 > 10)
        //  r2 = 10;

        ymax *= r2;

        DBG(msg("Extension RDI log (2) : ymin={},ymax={}", ymin, ymax));
      }
    }
    else
    {
      ymin -= dy * delta_y;
      if(calc_titre().empty())
        ymax += dy * delta_y;
      else
        ymax += dy * (delta_y + 0.15);
    }
  }

  mutable float xmin = 0, xmax = 0, ymin = 0, ymax = 0;

  int rendre00() const
  {
    auto config = axes.get_config();
    if(calc_minmax(config, xmin, ymin, xmax, ymax))
      return -1;

    if(!a_rdi)
      extension_rdi(xmin, ymin, xmax, ymax);
    return 0;
  }

  std::string calc_titre() const
  {
    auto config = axes.get_config();

    std::string res = /*config.*/titre;
    if((courbes.size() == 1) && (res.empty()))
      res = courbes[0].impl->nom;

    return res;
  }

  // TODO : gestion cas d'échec calc_minmax
  Canva rendre0(Canva canva_) const
  {
    auto config = axes.get_config();

    config.legende.position = pos_cartouche;

    config.series.clear();

    if((!courbes.empty()) && (!courbes[0].impl->nom.empty()))
    {
      for(auto &c: courbes)
      {
        ConfigSerie cs;
        cs.nom          = c.impl->nom;
        cs.couleur      = c.impl->couleur;
        cs.remplissage  = c.impl->remplissage;
        cs.trait        = c.impl->trait;
        cs.marqueur     = c.impl->marqueur;
        if(!cs.nom.empty())
          config.series.push_back(cs);
      }
    }

    config.legende.afficher = true;
    if((courbes.size() == 1) &&
        ((config.titre.empty())
          || (config.titre == courbes[0].impl->nom)))
    {
      //infos("Pas de titre et une seule courbe -> pas de légende.");
      //config.titre = courbes[0].impl->nom;
      config.legende.afficher = false;
    }

    config.titre = calc_titre();

    //if(config.legende_dim)
    //config.legende_dim  = 6;
    axes.configure(config);

    DBG(msg("rendu axes, ymin={},ymax={},diff={}...", ymin, ymax, ymin-ymax));

    // Pb en mode logarithmique, si ymin = 1e-30, ymax = 2, on ne peut pas coder correctement ymin-ymax
    // Dans axes, le calcul suivant est fait :
    // ymax + (ymin-ymax) = 0....

    //Canva canva = axes.rendre(canva_, {xmin, ymax, xmax-xmin, ymin-ymax});

    // TODO : création d'un PTI ici

    Canva canva = axes.rendre(canva_, xmin, xmax, ymin, ymax);
    DBG(msg("ok.");)

    return canva;
  }

  void rendre1(Canva canva_, Canva canva) const
  {
    DBG(msg("dessin courbes..."));
    for(auto &c: courbes)
    {
      auto &ci = *(c.impl);
      auto &x = ci.x;
      auto &y = ci.y;

      auto npts = x.rows();
      tsd_assert(npts == y.rows());


      DBG(msg("courbe : {} points.", npts));

      canva.set_couleur(ci.couleur);
      canva.set_epaisseur(ci.epaisseur);

      if(ci.Z.rows() > 0)
      {
        unsigned int nr = ci.Z.rows(), nc = ci.Z.cols();

        float xpas = ci.rdi_z.l  * 1.0f / nr;
        float ypas = ci.rdi_z.h * 1.0f / nc;

        //infos("xpas = %f, width = %d", xpas, ci.rdi_z.l);

        ArrayXXf Zn = ci.Z - ci.Z.minCoeff();
        Zn /= Zn.maxCoeff();

        auto x0 = ci.rdi_z.x, y0 = ci.rdi_z.y;


        auto cmap = tsd::vue::cmap_parse(ci.format);//"jet");

        for(auto i = 0u; i < nr; i++)
        {
          for(auto j = 0u; j < nc; j++)
          {
            //float R = 0, V = 0, B = 0;
            //cmap_mono(Zn(i,j), R, V, B);
            //cmap_jet(Zn(i,j), R, V, B);
            //Couleur coul{255*R, 255*V, 255*B};
            auto coul = cmap->couleur(Zn(i,j));
            canva.set_remplissage(true, coul);
            canva.set_couleur(coul);
            canva.rectangle(x0 + i * xpas, y0 + j * ypas,
                x0 + (i+1) * xpas, y0 + (j+1) * ypas);
          }
        }
        canva.set_remplissage(false);

      }

      if(c.impl->σ.rows() > 0)
      {
        if(c.impl->σ.rows() != npts)
        {
          msg_erreur("Affichage courbe avec σ : le nombre de points doit être identique pour σ et y.");
          continue;
        }
        auto cr = ci.couleur.eclaircir(0.3);
        canva.set_remplissage(true, cr);
        canva.set_couleur(cr);
        for(auto i = 0; i + 1 < npts; i++)
        {
          // A clarifier

          auto cvt = [&](float y, float s) -> std::tuple<float,float>
          {
            float y1 = y + s;
            // y2 = y^2 / (y + s)
            float y2 = std::pow(10.0f, 2 * log10(y) - log10(y1));
            if(std::isnan(y2) || std::isinf(y2))
              y2 = 0;
            return {y2, y1};
          };

          // attention, y2 peut être négatif !!!

          auto [y0m,y0p] = cvt(y(i), c.impl->σ(i));
          auto [y1m,y1p] = cvt(y(i+1), c.impl->σ(i+1));

          if((y0p != 0) && (y1p != 0))
            canva.remplissage_vertical(x(i), y0m, x(i+1),  y1m, y0p, y1p);
        }
        canva.set_remplissage(false, cr);
        canva.set_couleur(ci.couleur);
      }

      if(ci.trait == Trait::MINMAX)
      {
        auto cr = ci.couleur.eclaircir(0.3);
        canva.set_couleur(ci.couleur);
        canva.set_remplissage(true, cr);
        for(auto i = 0; i + 1 < npts; i++)
          canva.remplissage_vertical(x(i), ci.ymin(i), x(i+1), ci.ymin(i+1), ci.ymax(i), ci.ymax(i+1));
        canva.set_remplissage(false, cr);
      }
      else if((ci.trait == Trait::LIGNE)
              || (ci.trait == Trait::DOTTED))
      {

        //////////////////////////////////////////////////////////////////////////////////
        // 1er pré-traitement :
        // on regarde les points qui sont vraiment visibles, on peut ignorer les autres
        //
        // 2eme pré-traitement :
        // Si il reste beaucoup de points, on fait une décimation spéciale
        // (calcul min / max sur des intervalle donnés)

        ArrayXf xr1, yr1;

        if(npts > 2000)
        {
          bool monotone = true;
          int idmin = -1, idmax = -1;

          for(auto i = 0; i + 1 < npts; i++)
          {
            // Vérification suite monotone, sinon on fait pas ça !
            if(x(i) > x(i+1))
            {
              monotone = false;
              break;
            }
            if((x(i) >= xmin) && (idmin == -1))
              idmin = i;
            if((x(i) > xmax) && (idmax == -1))
              idmax = i;
          }
          if((idmin >= 0) && (idmax >= 0))
          {
            int npts3 = idmax - idmin + 1;
            // Si monotone, et si la réduction en vaut la peine
            if(monotone && (npts3 < 0.5 * npts)) // Un peu arbitraire
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
        ArrayXf xr, yrmin, yrmax, yrmoy;
        bool dessine_sigma = false;
        if(npts >= 2000)
        {
          dessine_sigma = true;
          int R = floor(npts / 1000);
          int npts2 = npts / R;

          //msg("plot : npts={} -> R={}, npts2={}", npts, R, npts2);

          xr.resize(npts2);
          yrmin.resize(npts2);
          yrmax.resize(npts2);
          yrmoy.resize(npts2);

          for(auto i = 0; i < npts2; i++)
          {
            xr(i) = x(i * R);
            yrmin(i) = yrmax(i) = y(i * R);
            yrmoy(i) = y.segment(i * R, R).mean();
            for(auto j = 1; j < R; j++)
            {
              yrmax(i) = std::max(yrmax(i), y(i*R+j));
              yrmin(i) = std::min(yrmin(i), y(i*R+j));
            }
          }
          npts = npts2;
          x = xr;
          y = yrmoy;

          // Fait après remplissage éventuel
          /*Couleur cr1 = axes.get_config().est_mat_sur_clair() ? ci.couleur.eclaircir(0.7) : ci.couleur.assombrir(0.7);
          canva.set_couleur(cr1);
          for(auto i = 0; i + 1 < npts; i++)
          {
            if((x(i+1) <= xmax) && (x(i) >= xmin))
              canva.remplissage_vertical(x(i), yrmax(i), x(i+1),  yrmax(i+1), yrmin(i), yrmin(i+1));
          }
          canva.set_remplissage(false);*/
        }


        if(ci.remplissage)
        {
          Couleur cr1 = axes.get_config().est_mat_sur_clair() ?
              ci.couleur.eclaircir(0.1) : ci.couleur.assombrir(0.1);
          canva.set_couleur(cr1);
          //auto cr = ci.couleur.eclaircir(0.3);

          float vmin = ci.remplissage_vers_ymin ? ymin + (ymax-ymin)/200 : ci.remplissage_vmin;

          for(auto i = 0; i + 1 < npts; i++)
          {
            if((x(i+1) <= xmax) && (x(i) >= xmin) && (y(i) >= ymin) && (y(i+1) >= ymin))
            {
              canva.remplissage_vertical(x(i), y(i), x(i+1),  y(i+1), vmin, vmin);
            }
          }
          canva.set_remplissage(false);
          canva.set_couleur(ci.couleur);
        }

        if(dessine_sigma)
        {
          Couleur cr1 = axes.get_config().est_mat_sur_clair() ?
              ci.couleur.eclaircir(0.2) : ci.couleur.assombrir(0.2);
          canva.set_couleur(cr1);
          for(auto i = 0; i + 1 < npts; i++)
          {
            if((x(i+1) <= xmax) && (x(i) >= xmin))
              canva.remplissage_vertical(x(i), yrmax(i), x(i+1),  yrmax(i+1), yrmin(i), yrmin(i+1));
          }
          canva.set_remplissage(false);
        }

        bool dotted = ci.trait == Trait::DOTTED;
        canva.set_dotted(dotted);
        canva.set_couleur(ci.couleur);
        for(auto i = 0; i + 1 < npts; i++)
        {
          if((x(i+1) <= xmax) && (x(i) >= xmin) && (y(i) >= ymin) && (y(i+1) >= ymin))
            canva.ligne(x(i), y(i), x(i+1), y(i+1));
        }
        canva.set_dotted(false);
      }
      else if(c.impl->trait == Trait::BATON)
      {
        for(auto i = 0; i < npts; i++)
        {
          if(!std::isinf(y(i)))
          {
            canva.ligne(x(i), y(i), x(i), 0.0f);
            //canva.marqueur({x(i), y(i)}, Marqueur::CERCLE, 5);
          }
        }
      }
      else if(c.impl->trait == Trait::HISTO)
      {
        if(x.size() > 1)
        {
          canva.set_remplissage(true, ci.couleur);
          auto dx = x(1) - ci.x(0);
          for(auto i = 0; i < npts; i++)
          {
            if(!std::isinf(y(i)))
              canva.rectangle(x(i), y(i), x(i) + dx, 0);
          }
          canva.set_remplissage(false, ci.couleur);
        }
      }




      // Comment faire pour accumuler ?
      //  Idée : adapter par rapport à la densité maximale
      //         -> il faudrait pouvoir dessiner en flottant et normaliser
      // Ou alors, on ajoute 1 unité à chaque itération
      //  -> puis on normalise

      if(ci.accu)
      {
        //msg("figure : rendu accu ({} points)", npts);
        ArrayXcf pts(npts);
        canva.set_couleur(ci.couleur);
        for(auto i = 0; i < npts; i++)
        {
          pts(i).real(ci.x(i));
          pts(i).imag(ci.y(i));
        }
        canva.dessine_accu(pts);
      }

      if(ci.marqueur == Marqueur::AUCUN)
        continue;


      for(auto i = 0; i < npts; i++)
      {
        if(ci.accu)
          continue;

        Couleur couleur = ci.couleur;

        if((ci.couleurs_points_rvb.rows() > 0) && (i < ci.couleurs_points_rvb.cols()))
          couleur = Couleur{
          ci.couleurs_points_rvb(2, i),
          ci.couleurs_points_rvb(1, i),
          ci.couleurs_points_rvb(0, i)};

        canva.set_couleur(couleur);

        auto x = ci.x(i), y = ci.y(i);
        canva.set_remplissage(true, couleur);
        canva.marqueur({x, y}, ci.marqueur, ci.dim_marqueur);
        canva.set_remplissage(false);
      }
    }
    axes.post_rendu(canva_);
    canva_utilisateur.forward(canva);
  }

  void rendre(Canva canva_) const
  {
    // Détection ici ????
    if(rendre00())
      return;
    auto canva = rendre0(canva_);
    rendre1(canva_, canva);
  }

  Figure::Courbe plot(const ArrayXf &x, const ArrayXf &y, const std::string &format,
      const std::string &titre);
};




void Figure::Courbe::def_remplissage(bool actif, bool vers_vmin, float vmin)
{
  impl->remplissage           = actif;
  impl->remplissage_vers_ymin = !vers_vmin;
  impl->remplissage_vmin      = vmin;
}

void Figure::Courbe::def_epaisseur(int ep)
{
  impl->epaisseur = ep;
}

void Figure::Courbe::def_couleur(const tsd::vue::Couleur &coul)
{
  impl->couleur = coul;
}

void Figure::Courbe::def_σ(IArrayXf σ)
{
  impl->σ = σ;
}

void Figure::Courbe::def_légende(const std::string &titre)
{
  impl->nom = titre;
}

void Figure::Courbe::def_dim_marqueur(int dim)
{
  impl->dim_marqueur = dim;
}

// Définit la couleur de chacun des points de la courbe
void Figure::Courbe::def_couleurs(IArrayXf c, const std::string cmap_nom)
{
  ArrayXXf rvb(3, c.rows());

  auto cmap = tsd::vue::cmap_parse(cmap_nom);

  for(auto i = 0; i < c.rows(); i++)
  {
    float R = 0, V = 0, B = 0;
    cmap->calc(c(i), R, V, B);
    rvb(0, i) = R;
    rvb(1, i) = V;
    rvb(2, i) = B;
  }
  impl->couleurs_points_rvb = 255 * rvb;
}

Figure::Figure(const std::string &nom)
{
  impl = std::make_shared<Impl>();
  impl->nom = nom;
}

Axes Figure::axes()
{
  return impl->axes;
}

void Figure::clear()
{
  impl->courbes.clear();
  impl->canva_utilisateur.clear();
}

struct Figures::Impl: Rendable
{
  int n = -1, m = -1, pos = 1;
  std::deque<Figure> subplots;

  std::tuple<int, int> get_nm() const;


  Impl(int n = 1, int m = 1)
  {
    this->n = n;
    this->m = m;
    pos = 1;
  }

  Figure *gcs()
  {
    if(subplots.empty())
      subplots.push_back(Figure());

    if((pos-1) >= (int) subplots.size())
    {
      msg_erreur("gcs() : pos = {}, subplots.size() = {}", pos, subplots.size());
      return &(subplots[0]);
    }
    return &(subplots[pos - 1]);
  }

  void rendre(Canva canva) const
  {
    //infos("Rendu figure en cours...");
    int nsubs = subplots.size();

    auto [nrows, ncols] = get_nm();

    tsd_assert((ncols >= 1) && (nrows >= 1));
    tsd_assert_msg(nsubs <= (ncols * nrows), "nsubs = {}, m = {}, n = {}", nsubs, ncols, nrows);

    //std::vector<Image> sO;
    auto rdi = canva.get_rdi();
    int mx = rdi.l / ncols,
        my = rdi.h / nrows;

    DBG(msg("parcours subplots..."));
    int i = 0;
    for(auto &s: subplots)
    {
      int px = (i % ncols) * mx;
      int py = (i / ncols) * my;
      Canva sub = canva.clip(Rect{px, py, mx, my}, Rect{0, 0, mx, my});
      s.rendre(sub);
      i++;
    }
      //rendre_subplot(s, sO, mx, my);
    //DBG(msg("fin parcours subplots."));

    DBG(msg("ok.");)

    //if((ncols == 1) && (nrows == 1)) //&& (sO.size() == 1))
      //return;

    /*Image O(sx, sy);
    O.remplir(Couleur::Blanc);

    for(auto i = 0u; i < sO.size(); i++)
    {
      int px = (i % ncols) * mx;
      int py = (i / ncols) * my;
      O.puti({px, py, mx, my}, sO[i]);
    }

    DBG(msg("Fait.");)
    return O;*/
  }


};

sptr<const Rendable> Figures::rendable() const
{
  return impl;
}


void Figures::afficher(const std::string &titre, const Dim &dim) const
{
  Dim dim2 = dim;
  if(dim.l < 0)
  {
    dim2.l = impl->m * 800;
    dim2.h = impl->n * 500;
  }
  ARendable::afficher(titre, dim2);
}

Figures::Figures(int n, int m)
{
  impl = std::make_shared<Impl>(n, m);
}

void Figures::clear()
{
  *impl = Impl();
}

Figure Figures::gf(int sel)
{
  if(sel >= (int) impl->subplots.size())
    echec("Figures::gcf({}) : index invalide (n figures = {})", sel, impl->subplots.size());
  return impl->subplots[sel];
}

Figure Figures::subplot(int n, int m, int pos_)
{
  int pos = pos_;
  if(((impl->n != n) || (impl->m != m)) && (impl->subplots.size() > 0))
  {
    impl->subplots.clear();
    //clf();
  }
  impl->n   = n;
  impl->m   = m;

  if(pos <= 0)
  {
    impl->subplots.push_back(Figure());
    pos = impl->subplots.size();
  }

  if(pos > n * m)
  {
    msg_erreur("subplot : n = {}, m = {}, pos = {} (pos prm = {}), nbsubplots={}.", n, m, pos, pos_, impl->subplots.size());
  }

  //infos("subplot : %d, %d, %d", n, m, pos);

  while(pos > (int) impl->subplots.size())
    impl->subplots.push_back(Figure());

  /*if(pos > impl->pos)
  {
    impl->subplots.push_back(SubPlot());
  }*/
  impl->pos = pos;
  return impl->subplots[pos-1];
}

Figure Figures::gcf()
{
  return impl->subplots.back();
}

Figure Figures::subplot(int i)
{
  if(i < 0)
  {
    impl->n = impl->m  = -1;
    impl->subplots.push_back(Figure());
    impl->pos = impl->subplots.size();
    return impl->subplots.back();
  }

  int n = i / 100;
  int m = (i - 100 * n) / 10;
  int pos = i % 10;

  return subplot(n, m, pos);
}


std::vector<Figure::Courbe> &Figure::courbes()
{
  return impl->courbes;
}

void Figure::def_rdi_min(const Rectf &rdi)
{
  tsd_assert(impl);
  impl->a_rdi_min = true;
  impl->rm_rdi    = rdi;

}

// TODO : ordre pas cohérent
void Figure::def_rdi(const Rectf &rdi)
{
  impl->a_rdi   = true;
  impl->rm_rdi  = rdi;
}

Rectf Figure::get_rdi() const
{
  return impl->rm_rdi;//{impl->rm_xmin, impl->rm_ymin, impl->rm_xmax - impl->rm_xmin, impl->rm_ymax - impl->rm_ymin};
}

Figure::Courbe Figure::plot_minmax(const ArrayXf &x, const ArrayXf &y1, const ArrayXf &y2)
{
  Figure::Courbe res;
  res.impl = std::make_shared<Figure::Courbe::Impl>();
  res.impl->x    = x;
  res.impl->ymin = y1;
  res.impl->ymax = y2;
  res.impl->y = (y1 + y2) / 2;
  res.impl->trait = Trait::MINMAX;

  impl->courbes.push_back(res);
  return res;
}

Figure::Courbe Figure::plot_iq_int(const ArrayXcf &z, const std::string &format, const std::string &titre)
{
  axes().set_isoview(true);
  auto x = z.real();
  auto y = z.imag();
  auto res = impl->plot(x, y, format, titre);
  return res;
}


Figure::Courbe Figure::plot_img(float xmin, float xmax, float ymin, float ymax, IArrayXXf &Z, const std::string &format)
{
  Figure::Courbe res;
  res.impl = std::make_shared<Figure::Courbe::Impl>();

  res.impl->Z = Z;
  res.impl->format = format;
  res.impl->rdi_z.x = xmin;
  res.impl->rdi_z.y = ymin;
  res.impl->rdi_z.l = xmax - xmin;// + 1;
  res.impl->rdi_z.h = ymax - ymin;// + 1;

  impl->courbes.push_back(res);

  return res;
}

Figure::Courbe Figure::plot_img(IArrayXXf &Z, const std::string &format)
{
  return plot_img(0, Z.rows() - 1, 0, Z.cols() - 1, Z, format);
}



void Figure::plot(const float &x, const float &y, const std::string &format)
{
  ArrayXf vx(1), vy(1);
  vx(0) = x;
  vy(0) = y;
  plot(vx, vy, format);
}

Figure::Courbe Figure::plot_int(const ArrayXf &y, const std::string &format, const std::string &titre)
{
  return impl->plot(ArrayXf(), y, format, titre);
}




Figure::Courbe Figure::plot_psd_int(IArrayXcf y, float fe, const std::string &format, const std::string &titre)
{
  if(y.rows() == 0)
    return Courbe();

  auto max_coef_imag = y.imag().abs().maxCoeff();
  bool est_reel = max_coef_imag == 0;


  ArrayXcf yc = y;
  ArrayXf fen = tsd::filtrage::fenetre("hn", y.rows(), false);
  ArrayXcf ycf = yc * fen;

  ArrayXf Y;

  Y = 10.0f * ((tsd::fourier::fft(ycf)).abs2().log10());
  ArrayXf freq;

  if(est_reel)
  {
    Y = Y.head(Y.rows() / 2).eval();
    freq = linspace(0, fe/2, Y.rows());
  }
  else
  {
    Y = (tsd::fourier::fftshift(Y)).eval();
    freq = linspace(-fe/2, fe/2, Y.rows());
  }

  impl->axes.get_config().axe_horizontal.label = "Fréquence";
  //if(fe != 1)
    //impl->axes.get_config().axe_horizontal.label += " (Hz)";

  return plot(freq, Y, format, titre);
}

Figure::Courbe Figure::plot_int(const ArrayXf &x, const ArrayXf &y, const std::string &format, const std::string &titre)
{
  auto res = impl->plot(x, y, format, titre);
  return res;
}

Figure::Courbe Figure::Impl::plot(const ArrayXf &x, const ArrayXf &y_, const std::string &format_, const std::string &titre)
{
  Figure::Courbe res;
  res.impl = std::make_shared<Figure::Courbe::Impl>();
  res.impl->nom = titre;

  std::string format = format_;

  if(format.empty())
  {
    int nc = courbes.size();
    std::vector<std::string> fdef = {"b-", "g-", "r-", "y-", "c-", "m-", "k-"};
    format = fdef[nc % 7];
  }

  auto n = y_.rows();

  if(n > 1e6)
  {
    msg_avert("plot : {} échantillons.", n);
  }

  ArrayXf y = y_;

  if(n == 0)
  {
    DBG(msg_avert("plot : y.rows() == 0. Titre = [{} / {}]", nom, res.impl->nom);)
    return res;
  }

  auto s = this;

  ArrayXf X;
  if(x.rows() == 0)
    X = linspace(0, n - 1, n);
  else
    X = x;

  if(X.rows() != n)
  {
    msg_erreur("Figure::plot(x,y,...) : x.rows() = {}, y.rows() = {}, y_.rows() = {}", X.rows(), y.rows(), y_.rows());
    return res;
  }
  tsd_assert(X.rows() == y.rows());

  res.impl->x       = X;
  res.impl->y       = y;
  res.impl->format  = format;

  auto present = [](const std::string &s, char c)
  {
    return s.find(c) != std::string::npos;
  };

  static const Couleur
    BleuSombre{0,0,180},
    VertSombre{0,100,0},
    RougeSombre{210,0,0},
    CyanSombre{0,192,192},
    VioletSombre{192,0,192},
    JauneSombre{192,192,0},
    MarronSombre{128,70,0},
    OrangeSombre{250,140,0};


  struct CodeCouleur {char code; Couleur couleur;};
  CodeCouleur codes[] =
  {
      {'b', BleuSombre},
      {'g', VertSombre},
      {'r', RougeSombre},

      {'y', JauneSombre},
      {'c', CyanSombre},
      {'m', VioletSombre},
      {'a', OrangeSombre},
      {'k', Couleur::Noir},
      {'w', Couleur::Blanc}

  };

  int nc = courbes.size();
  res.impl->couleur = codes[nc % 7].couleur;//fdef[nc % 6];

  for(auto c: codes)
    if(present(format, c.code))
      res.impl->couleur = c.couleur;

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

  for(auto c: codest)
    if(present(format, c.code))
      res.impl->trait = c.trait;


  struct CodeMarqueur {char code; Marqueur marqueur;};
  CodeMarqueur codesm[] =
  {
      //{'+', Marqueur::PLUS},
      {'o', Marqueur::CERCLE},
      {'*', Marqueur::ETOILE},
      {'.', Marqueur::POINT},
      //{'x', Marqueur::CROIX},
      {'s', Marqueur::CARRE},
      {'d', Marqueur::DIAMANT},
  };

  for(auto c: codesm)
    if(present(format, c.code))
      res.impl->marqueur = c.marqueur;

  s->courbes.push_back(res);
  return res;
}

void Figure::def_pos_legende(const std::string &code)
{
  impl->pos_cartouche = code;
}

void Figure::titre(const std::string &titre_global)
{
  impl->/*axes.get_config().*/titre = titre_global;
}

void Figure::titres(const std::string &titre_global,
                    const std::string &axe_x,
                    const std::string &axe_y)
{
  auto &c = impl->axes.get_config();
  c.axe_horizontal.label  = axe_x;
  c.axe_vertical.label    = axe_y;
  impl->titre = titre_global;
  //c.titre                 = titre_global;
}


void Figure::attente_ihm()
{
  stdo.fin();
}

std::tuple<int, int> Figures::Impl::get_nm() const
{
  int nsubs = subplots.size();
  if(m == -1)
  {
    if(nsubs == 1)
      return {1, 1};
    else if(nsubs == 2)
      return {2, 1};
    else if(nsubs == 3)
      return {3, 1};
    else
    {
      int m = (int) floor(std::sqrt(nsubs));
      int n = (nsubs + m - 1) / m;
      return {n, m};
    }
  }
  return {m, n};
}



void Figure::def_nom(const std::string &s)
{
  impl->nom = s;
}

std::string Figure::lis_nom() const
{
  return impl->nom;
}

Figure Figure::clone() const
{
  Figure res;
  res.impl = std::shared_ptr<Impl>(new Figure::Impl(*impl));
  res.impl->axes = impl->axes.clone();
  return res;
}

Canva Figure::canva_pixel(const Dim &allocation)
{
  auto sz = axes().get_dim_graphique(allocation);
  auto a = canva();
  Rectf r = get_rdi();
  return impl->canva_utilisateur.clip(r, Rect{0, sz.h, sz.l, -sz.h});
}

Canva Figure::canva()
{
  return impl->canva_utilisateur;
}

void Figure::def_echelle(bool log_x, bool log_y)
{
  impl->log_x = log_x;
  impl->log_y = log_y;
}

/*void Figure::rendre(Canva canva) const
{
  impl->rendre(canva);
}*/

sptr<const Rendable> Figure::rendable() const
{
  return impl;
}

void Figure::rendre1(Canva canva_, Canva canva) const
{
  return impl->rendre1(canva_, canva);
}

Canva Figure::rendre0(Canva canva_) const
{
  return impl->rendre0(canva_);
}

int Figure::rendre00() const
{
  return impl->rendre00();
}


}

