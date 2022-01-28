#include "tsd/tsd.hpp"
#include "tsd/figure.hpp"
#include <deque>

#define DBG(AA)


using namespace std;

namespace tsd::vue
{

struct Canva::PointIntermediaire
{
  Image img;
};

struct Canva::Impl
{

  struct Pinceau
  {
    Align alignement_vertical   = Align::DEBUT,
          alignement_horizontal = Align::DEBUT;
    Couleur couleur_trait{0, 0, 0},
            couleur_remplissage{255,255,255};
    int epaisseur = 1;
    bool remplir = false;
    Orientation orient = Orientation::HORIZONTALE;
    // Une chaine de caractère : fonte précisée, ou bien,
    // un rectangle englobant
    float dim_fonte         = -1;
    bool dotted             = false;
    float alpha             = 1;
    //bool coord_mode_pixels  = false;
  };


  // Permet de passer d'un système de coordonnées à un autre
  struct ClipGen
  {
    // Position dans le canva parent
    Rectf target{0,0,0,0};
    // RDI = système de coordonnée virtuel
    Rectf rdi{0,0,1,1};

    bool log_x = false, log_y = false;

    inline float v2c_x(float x) const
    {
      if(log_x)
        x = log10(x);

      if(rdi.l == 0)
        echec("ClipGen::v2c_x : rdi.l = 0.");

      return target.x + target.l * (x - rdi.x) / rdi.l;
    }

    inline float v2c_y(float y) const
    {
      if(log_y)
        y = log10(y);

      if(rdi.h == 0)
        echec("ClipGen::v2c_x : rdi.h = 0.");

      return target.y + target.h * (y - rdi.y) / rdi.h;
    }

    Pointf c2v(const Pointf &p) const
    {
      if((target.l == 0) || (target.h == 0))
        echec("ClipGen::c2v : target = {}.", target);


      return {(p.x - target.x) * rdi.l / target.l + rdi.x,
              (p.y - target.y) * rdi.h / target.h + rdi.y};
    }

    // Système interne ("virtuel") vers externe ("concret")
    Pointf v2c(const Pointf &p) const
    {
      return {v2c_x(p.x), v2c_y(p.y)};
    }

    Dimf v2c_l(const Dimf &d) const
    {
      return {v2c_lx(d.l), v2c_ly(d.h)};
    }

    // Attention : n'a pas de sens en échelle logarithmique
    float v2c_lx(float largeur) const
    {
      return largeur * target.l / rdi.l;
    }

    float v2c_ly(float hauteur) const
    {
      return hauteur * target.h / rdi.h;
    }

    ClipGen vue(const Rectf &rdi_) const
    {
      ClipGen res;
      res.target = this->rdi;
      res.rdi    = rdi_;
      return res;
    }

    ClipGen clip(const Rectf &rdi_, const Rectf &vue) const
    {
      ClipGen res;
      res.target = rdi_;
      res.rdi    = vue;
      return res;
    }
  };


  struct CanvaElement
  {
    enum Type
    {
      RECT, LIGNE, FLECHE, CHAINE, CERCLE, MARQUEUR, RVEC, ACCU, IMAGE, PTI, ELLIPSE
    } type = MARQUEUR;
    std::string chaine;
    Pointf p0, p1;
    Dimf dim;
    float r = 3, y2 = 0, y3 = 0;
    float a, b;
    //int dim_pixels = 3;
    int dim_pixels = 5;


    Marqueur marqueur = Marqueur::CARRE;
    Pinceau pinceau;
    ArrayXcf accu_points;
    Image img;
    sptr<PointIntermediaire> pti;


    Rectf get_rdi() const
    {
      switch(type)
      {
      case RECT:
      case LIGNE:
      case FLECHE:
        return Rectf::englobant(p0, p1);
      case IMAGE:
      {
        // Pb si coordonnées pixel !!!!
        auto res = Rectf::englobant(p0, p1);// + Pointf(1,1));

        //msg("get_rdi(image) : res = {}", res);
        // Bricolage !
        res.l++;
        res.h++;

        return res;
      }
      case CHAINE:
      {
        return Rectf(p0, Dimf{1.0f,1.0f}); // TODO
      }
      case MARQUEUR:
        return Rectf(p0, Dimf{0.0f,0.0f}); // TODO
      case CERCLE:
        return Rectf(p0 - Pointf(r/2,r/2), Dimf{r/2,r/2});
      default:
        return Rectf(p0, Dimf{1.0f,1.0f}); // TODO
      }
    }

  };


  ///sptr<Impl> parent;
  std::weak_ptr<Impl> parent;
  ClipGen clip;

  Pinceau pinceau;
  std::deque<CanvaElement> elements;
  Image img_fond;


  Dim get_allocation() const
  {
    auto p = parent.lock();
    if(p)
    {
      auto dp = p->get_allocation();
      Dim res;
      Rectf rect_parent = p->clip.rdi;
      res.l = dp.l * (clip.target.l / rect_parent.l);
      res.h = dp.h * (clip.target.h / rect_parent.h);
      return res;
    }
    return {(int) clip.target.l, (int) clip.target.h};
  }

  Rectf get_rdi() const
  {
    return clip.rdi;
  }

  float pixels_par_lsb_x() const
  {
    return clip.v2c_lx(1.0f);
  }
  float pixels_par_lsb_y() const
  {
    return clip.v2c_ly(1.0f);
  }

  // pixel -> pos
  Pointf pixel2pos(int x, int y)
  {
    return clip.c2v({(float)x,(float)y});
  }

  // pos -> pixel
  Point pos(const Pointf &p) const
  {
    auto res = clip.v2c(p);
    return {(int) res.x, (int) res.y};
  }

  Pointf posf(const Pointf &p) const
  {
    return clip.v2c(p);
  }

  void marqueur(Image O, CanvaElement &d)
  {
    auto ep = d.pinceau.epaisseur;
    auto pos0 = pos(d.p0);

    O.def_epaisseur(ep);
    O.def_couleur_dessin(d.pinceau.couleur_trait);
    O.def_couleur_remplissage(d.pinceau.couleur_remplissage);


    //DBG_AXES(msg("marqueur r={}...", d.r);)
    int r2 = (int) round(d.dim_pixels / sqrt(2.0f));
    if(d.marqueur == Marqueur::CARRE)
    {
      //msg("Marqueur carré : r = {}, r2 = {}", d.r, r2);
      Point p0{(int) (pos0.x-r2),(int) (pos0.y-r2)},
            p1{(int) (pos0.x+r2),(int) (pos0.y+r2)};
      O.rectangle_plein(p0, p1);
      O.rectangle(p0, p1);
    }
    else if(d.marqueur == Marqueur::ETOILE)
    {
      O.ligne_aa(pos0.x-r2, pos0.y-r2, pos0.x+r2, pos0.y+r2);
      O.ligne_aa(pos0.x-r2, pos0.y, pos0.x+r2, pos0.y);
      O.ligne_aa(pos0.x, pos0.y-r2, pos0.x, pos0.y+r2);
      O.ligne_aa(pos0.x-r2, pos0.y+r2, pos0.x+r2, pos0.y-r2);
    }
    else if(d.marqueur == Marqueur::CROIX)
    {
      O.ligne_aa(pos0.x-r2, pos0.y, pos0.x+r2, pos0.y);
      O.ligne_aa(pos0.x, pos0.y-r2, pos0.x, pos0.y+r2);
    }
    else if(d.marqueur == Marqueur::DIAMANT)
    {
      O.ligne_aa(pos0.x-r2, pos0.y, pos0.x, pos0.y-r2);
      O.ligne_aa(pos0.x, pos0.y-r2, pos0.x+r2, pos0.y);
      O.ligne_aa(pos0.x+r2, pos0.y, pos0.x, pos0.y+r2);
      O.ligne_aa(pos0.x, pos0.y+r2, pos0.x-r2, pos0.y);
    }
    else if(d.marqueur == Marqueur::POINT)
    {
      O.point(pos0);
    }
    else if(d.marqueur == Marqueur::CERCLE)
    {
      //msg("Marqueur cercle : r = {}", d.r);
      O.cercle_plein(pos0, d.dim_pixels);
      O.cercle(pos0, d.dim_pixels);
    }
    else
    {
      msg_erreur("Marqueur de type inconnu ({})", (int) d.type);
    }
  }

  static Rectf union_rdi(const Rectf &r1, const Rectf &r2)
  {
    Rectf r;
    r.x = min(r1.x, r2.x);
    r.y = min(r1.y, r2.y);
    r.l = max(r1.x + r1.l, r2.x + r2.l) - r.x;
    r.h = max(r1.y + r1.h, r2.y + r2.h) - r.y;
    return r;
  }

  Rectf calc_rdi_englobante() const
  {
    Rectf rdi;
    for(auto &d: elements)
    {
      Rectf r = d.get_rdi();


      if(std::isinf(r.l) || (std::isinf(r.h)))
      {
        echec("calc_rdi_englobante(): r={}, type={}, p0={}, p1={}", r, (int) d.type, d.p0, d.p1);
      }

      if(&d == &(elements.front()))
        rdi = r;
      else
        rdi = union_rdi(rdi, r);

      //msg(" -> r={}, rdi={}", r, rdi);
    }
    //msg("Calc RDI englobante : {}", rdi);
    return rdi;
  }


  static CanvaElement clip_elmt(const CanvaElement &elt, const ClipGen &clip)
  {
    CanvaElement elt2 = elt;
    elt2.p0   = clip.v2c(elt.p0);

    elt2.p1   = clip.v2c(elt.p1);
    elt2.y2   = clip.v2c_y(elt.y2);
    elt2.y3   = clip.v2c_y(elt.y3);
    if(elt.dim.l > 0)
      elt2.dim  = clip.v2c_l(elt.dim);
    elt2.a    = clip.v2c_lx(elt.a);
    elt2.b    = clip.v2c_ly(elt.b);

    // Si commenté : pb cercles
    // Si décommenté : pb marqueurs
    elt2.r    = clip.v2c_lx(elt.r);

    for(auto i = 0u; i < elt.accu_points.rows(); i++)
    {
      elt2.accu_points(i).real(clip.v2c_x(elt.accu_points(i).real()));
      elt2.accu_points(i).imag(clip.v2c_y(elt.accu_points(i).imag()));
    }
    return elt2;
  }


  void ajoute_element(const CanvaElement &elt)
  {
    auto p = parent.lock();
    if(p)
    {
      auto elt2 = clip_elmt(elt, clip);
      p->ajoute_element(elt2);
      return;
    }
    elements.push_back(elt);
  }

  void ellipse(const Pointf &p, float a, float b, float theta, const Pinceau &pinceau)
  {
    Impl::CanvaElement d;
    d.p0 = p;
    d.a  = a;
    d.b  = b;
    d.type = Impl::CanvaElement::ELLIPSE;
    d.pinceau = pinceau;
    ajoute_element(d);
  }

  void cercle(const Pointf &p, float r, const Pinceau &pinceau)
  {
    Impl::CanvaElement d;
    d.p0 = p;
    d.r  = r;
    //d.type = Impl::CanvaElement::CERCLE;
    d.a  = d.b = r;
    d.type = Impl::CanvaElement::ELLIPSE;
    d.pinceau = pinceau;
    ajoute_element(d);
  }

  void ligne(const Pointf &p0, const Pointf &p1, const Pinceau &pinceau)
  {
    if(std::isinf(p0.y) || std::isinf(p1.y))
      return;
    Impl::CanvaElement d;
    d.p0      = p0;
    d.p1      = p1;
    d.type    = Impl::CanvaElement::LIGNE;
    d.pinceau = pinceau;
    ajoute_element(d);
  }

  void fleche(const Pointf &p0, const Pointf &p1, const Pinceau &pinceau)
  {
    if(std::isinf(p0.y) || std::isinf(p1.y))
      return;
    Impl::CanvaElement d;
    d.p0      = p0;
    d.p1      = p1;
    d.type    = Impl::CanvaElement::FLECHE;
    d.pinceau = pinceau;
    ajoute_element(d);
  }

  void dessine_accu(const ArrayXcf &x, const Pinceau &pinceau)
  {
    Impl::CanvaElement d;
    d.type        = Impl::CanvaElement::ACCU;
    d.accu_points = x;
    d.pinceau     = pinceau;
    ajoute_element(d);
  }

  void dessine_img(const Pointf &p0, const Pointf &p1, Image img, const Pinceau &pinceau)
  {
    Impl::CanvaElement d;
    d.p0 = p0;
    d.p1 = p1;
    d.type = Impl::CanvaElement::IMAGE;
    d.img = img;
    d.pinceau = pinceau;
    ajoute_element(d);
  }

  void rectangle(const Pointf &p0, const Pointf &p1, const Pinceau &pinceau)
  {
    Impl::CanvaElement d;
    d.p0 = p0;
    d.p1 = p1;
    d.type = Impl::CanvaElement::RECT;
    d.pinceau = pinceau;
    ajoute_element(d);
  }

  void remplissage_vertical(const Pointf &p0, const Pointf &p1, float y2, float y3, const Pinceau &pinceau)
  {
    if(std::isinf(p1.y) || std::isinf(y2) || std::isinf(y3))
      return;

    Impl::CanvaElement d;
    d.p0 = p0;
    d.p1 = p1;
    d.y2 = y2;
    d.y3 = y3;
    d.type = Impl::CanvaElement::RVEC;
    d.pinceau = pinceau;
    ajoute_element(d);
  }


  void marqueur(const Pointf &p, Marqueur m, float dim_pixels, const Pinceau &pinceau)
  {
    Impl::CanvaElement d;
    d.p0          = p;
    d.dim_pixels  = dim_pixels;
    d.marqueur    = m;
    d.type        = Impl::CanvaElement::MARQUEUR;
    d.pinceau     = pinceau;
    ajoute_element(d);
  }



  void texte(const Pointf &p, const std::string &texte, const Dimf &dim_max, const Pinceau &pinceau)
  {
    Impl::CanvaElement d;
    d.p0      = p;
    d.dim     = dim_max;
    d.type    = Impl::CanvaElement::CHAINE;
    d.pinceau = pinceau;
    d.chaine  = texte;
    ajoute_element(d);
  }

  Image rendre(const Dim &dim_, const Rectf &rdi_, const Couleur &arp_)
  {
    Dim dim = dim_;

    if((dim.l < 0) && (clip.target.l > 0))
    {
      dim.l = clip.target.l;
      dim.h = clip.target.h;
    }


    if(dim.l < 0)
    {
      dim.l = 800;
      dim.h = 600;
    }

    Rectf rdi = rdi_;

    /*if((rdi.l == 0) && (rdi.h == 0))
    {
      rdi = clip.rdi;
    }*/


    if((rdi.l == 0) && (rdi.h == 0))
    {
      //rdi = calc_rdi_englobante();
      //msg("Canva/rendre() : rdi=0, allocation={}, rdi englobante={}", dim, rdi);
      rdi.x = 0;
      rdi.y = 0;
      rdi.l = dim.l;
      rdi.h = dim.h;
    }

    clip.target.x   = clip.target.y = 0;
    clip.target.l   = dim.l;
    clip.target.h   = dim.h;
    clip.rdi        = rdi;


    DBG(msg("canva rendre(): dim = {}, rdi = {}.", dim, rdi));
    Image O1;
    Couleur arp = arp_;
    arp.alpha = 0;

    if(!img_fond.empty())
      O1 = img_fond.redim(dim);
    else
    {
      O1.resize(dim.l, dim.h);
      O1.remplir(arp);
    }

    DBG(msg("Dessins ({})...", elements.size());)

    int cnt_rvec = 0, cnt_rec = 0, cnt_ligne = 0, cnt_chaine = 0, cnt_cercle = 0, cnt_marqueur = 0;

    for(auto &d: elements)
    {

      auto ep = d.pinceau.epaisseur;

      auto pos0 = pos(d.p0);
      auto pos1 = pos(d.p1);

      O1.def_epaisseur(ep);
      O1.def_couleur_dessin(d.pinceau.couleur_trait);
      O1.def_couleur_remplissage(d.pinceau.couleur_remplissage);

      // Sera possible à partir de GCC 11 🙂
      // using enum Impl::CanvaElement::Type;

      if(d.type == Impl::CanvaElement::PTI)
      {
        // Point intermédiaire
        // il faut sauvegarder l'état du dessin à ce moment là
        tsd_assert(d.pti);
        d.pti->img = O1.clone();
      }
      else if(d.type == Impl::CanvaElement::ACCU)
      {
        //msg("Rendu accu : (1)");
        DBG(msg("accu...");)
        auto c = d.pinceau.couleur_trait;
        auto n = d.accu_points.rows();
        auto dim = get_allocation();
        DBG(msg("axes2: rendu accu, allocation = {}, npts = {}", dim, n));
        ArrayXXf a = ArrayXXf::Zero(dim.l, dim.h);
        for(auto i = 0; i < n; i++)
        {
          Point p = pos({d.accu_points(i).real(), d.accu_points(i).imag()});
          int r = 1;
          if((p.x >= r) && (p.y >= r) && (p.x < dim.l-r) && (p.y < dim.h-r))
            a.block(p.x-r, p.y-r, 2*r+1, 2*r+1) += 1;
        }

        //msg("Rendu accu : (2)");

        auto mc = a.maxCoeff();
        DBG(msg("axes2: rendu accu, max = {}", mc));
        a *= 1.0f / (mc + 1e-20);
        for(auto y = 0; y < dim.h; y++)
        {
          for(auto x = 0; x < dim.l; x++)
          {
            if(a(x,y) > 0)
              O1.point({x,y}, c, a(x,y));
          }
        }
        //msg("Rendu accu : (fin)");
      }
      else if(d.type == Impl::CanvaElement::IMAGE)
      {
        auto r = Rect::englobant(pos0, pos1);

        // r = rectangle en pixel où doit doit être posé l'image
        // Calcul centre

        double cx = O1.sx() / 2, cy = O1.sy() / 2;
        cx = cx - r.x;
        cy = cy - r.y;
        cx *= d.img.sx();
        cy *= d.img.sy();
        cx /= r.l;
        cy /= r.h;

        //msg("Canva / dessin image rect = {}, dim img = {}, dim O1={}", r, d.img.get_dim(), O1.get_dim());
        //msg("Centre affichage dans coord image : {} x {}", cx, cy);



        auto i = d.img;

        // r : unités pixels

        // Ne prends en compte que ce qui est à l'intérieur de l'image
        if((r.x < 0) || (r.y < 0) || (r.x + r.l > O1.sx()) || (r.y + r.h > O1.sy()))
        {
          msg("CUT! (r = {}, pos0 = {}, pos1 = {}, d.p0={}, d.p1={})", r, pos0, pos1, d.p0, d.p1);
          if(r.x < 0)
          {
            // Coupe la partie gauche de l'image
            int cut = (-i.sx() * r.x) / r.l;
            i = i.sous_image(cut, 0, i.sx() - cut, i.sy());
            r.l += r.x;
            r.x = 0;
          }
          if(r.x + r.l > O1.sx())
          {
            // Coupe la partie droite de l'image
            int cut = (i.sx() * (r.x+r.l-O1.sx())) / r.l;
            i = i.sous_image(0, 0, i.sx() - cut, i.sy());
            r.l = O1.sx()-1-r.x;
          }
          if(r.y < 0)
          {
            // Coupe la partie haute de l'image
            int cut = (-i.sy() * r.y) / r.h;
            i = i.sous_image(0, cut, i.sx(), i.sy() - cut);
            r.h += r.y;
            r.y = 0;
          }
          if(r.y + r.h > O1.sy())
          {
            // Coupe la partie basse de l'image
            int cut = (i.sy() * (r.y+r.h-O1.sy())) / r.h;
            i = i.sous_image(0, 0, i.sx(), i.sy() - cut);
            r.h = O1.sy()-1-r.y;
          }

          msg("Après cut: rect = {}, dim img = {}, dim O1={}",
              r, i.get_dim(), O1.get_dim());
        }

        //msg("p0={}, p1={}", d.p0, d.p1);
        O1.puti(r, i);
      }
      else if(d.type == Impl::CanvaElement::RECT)
      {
        cnt_rec++;
        DBG(msg("rect...");)
        auto r = Rect::englobant(pos0, pos1);
        auto p0 = r.tl();
        auto p1 = r.br();

        for(auto k = 0; k < ep; k++)
          O1.rectangle({p0.x+k,p0.y+k}, {p1.x-k,p1.y-k});
        if(d.pinceau.remplir && (p0.x != p1.x) && (p0.y != p1.y))
          O1.rectangle_plein({p0.x+ep, p0.y+ep},{p1.x-ep,p1.y-ep});
      }
      else if(d.type == Impl::CanvaElement::RVEC)
      {
        cnt_rvec++;
        auto pos2 = pos({d.p0.x, d.y2});
        auto pos3 = pos({d.p1.x, d.y3});

        if((pos2.y != -1) && (pos3.y != -1))
        {
          for(auto x = pos0.x; x <= pos1.x; x++)
          {
            // Interpolation linéaire entre les points
            int yh = pos0.y + ((pos1.y - pos0.y) * (x - pos0.x) * 1.0f) / (pos1.x - pos0.x);
            int yb = pos2.y + ((pos3.y - pos2.y) * (x - pos0.x) * 1.0f) / (pos1.x - pos0.x);
            O1.ligne({x, yh}, {x, yb});
          }
        }
      }
      else if(d.type == Impl::CanvaElement::LIGNE)
      {
        cnt_ligne++;
        auto [x0,y0] = posf(d.p0);
        auto [x1,y1] = posf(d.p1);
        O1.ligne_aa(x0, y0, x1, y1, d.pinceau.dotted ? Image::POINTILLEE : Image::PLEINE);
      }
      else if(d.type == Impl::CanvaElement::FLECHE)
      {
        O1.fleche(pos0, pos1, d.pinceau.dotted ? Image::POINTILLEE : Image::PLEINE);
      }
      else if(d.type == Impl::CanvaElement::CHAINE)
      {
        cnt_chaine++;
        DBG(msg("chaine '{}'...", d.chaine););

        float fonte = 1;

        if(d.pinceau.dim_fonte != -1)
          fonte = d.pinceau.dim_fonte;

        tsd::vue::TexteConfiguration config;

        config.couleur    = d.pinceau.couleur_trait;
        config.scale      = fonte;
        config.thickness  = d.pinceau.epaisseur;
        config.dim_max    = Dim{O1.sx() - pos0.x, O1.sy() - pos0.y};


        if(d.dim.l > 0)
          config.dim_max.l = min((int) clip.v2c_lx(d.dim.l), config.dim_max.l);
        if(d.dim.h > 0)
          config.dim_max.h = min((int) clip.v2c_ly(d.dim.h), config.dim_max.h);

        if(d.pinceau.orient == Orientation::VERTICALE)
        {
          std::swap(config.dim_max.l, config.dim_max.h);
        }

        DBG(msg("  creation image, dmax = {}, dimax_in = {}...", config.dim_max, d.dim););

        auto i = tsd::vue::texte_creation_image(d.chaine, config);

        DBG(msg("  ok : {}...", i.get_dim()););
        if(i.empty())
          continue;

        if((pos0.x < 0) || (pos0.y < 0))
        {
          //msg_avert("Axes : affichage texte [{}], position invalide {}, position flt = {}.", d.chaine, pos0, d.p0);
          //msg("mode pixels = {}", d.pinceau.coord_mode_pixels);
        }

        if(d.pinceau.orient == Orientation::VERTICALE)
          i = i.rotation_90();

        //O1.puti_avec_gamma(pos0, i);

        int px, py;

        if(d.pinceau.alignement_horizontal == Align::DEBUT)
          px = pos0.x;
        else if(d.pinceau.alignement_horizontal == Align::FIN)
          px = pos0.x + d.p1.x - i.sx();
        else // Centre
          px = pos0.x - i.sx() / 2;


        if(d.pinceau.alignement_vertical == Align::DEBUT)
          py = pos0.y;
        else if(d.pinceau.alignement_vertical == Align::FIN)
          py = pos0.y + d.p1.y - i.sy();
        else // Centre
          py = pos0.y - i.sy() / 2 - 1;

        DBG(msg("  puti : @{}x{}, dim = {}.", px, py, i.get_dim()));


        /*if(d.chaine.substr(0, 3) == "0 dBm")
        {
          msg("  puti : @{}x{}, dim = {}.", px, py, i.get_dim());
        }*/

        //msg("  puti : '{}' @{}x{}, dim = {}, pos0.y={}, sy={}, d.p1.y={}, d.p0={}.", d.chaine, px, py, i.get_dim(), pos0.y, i.sy(), d.p1.y, d.p0);

        //O1.puti(Point{px, py}, i);
        try
        {
          O1.puti_avec_gamma(Point{px, py}, i);
        }
        catch(...)
        {
          msg_avert("Canva / dessin texte : échec puti avec gamma, chaine = '{}'.", d.chaine);
        }
        DBG(msg("  puti ok."));
      }
      else if(d.type == Impl::CanvaElement::CERCLE)
      {
        cnt_cercle++;
        DBG(msg("cercle r={}...", d.r);)
        int dx = d.r * pixels_par_lsb_x(),
            dy = d.r * pixels_par_lsb_y();
        DBG(msg("  -> dx={} pixels, dy = {} pixels", dx, dy);)

        Point tl{pos0.x-dx,pos0.y-dy}, br{pos0.x+dx,pos0.y+dy};

        //msg("Cercle")
        O1.ellipse(tl, br);

        if(d.pinceau.remplir)
          O1.ellipse_pleine(tl + Point{ep,ep}, br - Point{ep,ep});
          //O1.cercle_plein(pos0, d.r - ep); // rayon ????
      }
      else if(d.type == Impl::CanvaElement::ELLIPSE)
      {
        cnt_cercle++;
        DBG(msg("ellipse a={}, b={}...", d.a, d.b);)
        int dx = d.a * pixels_par_lsb_x(),
            dy = d.b * pixels_par_lsb_y();
        DBG(msg("  -> dx={} pixels, dy = {} pixels", dx, dy);)

        Point tl{pos0.x-dx,pos0.y-dy}, br{pos0.x+dx,pos0.y+dy};

        O1.ellipse(tl, br);

        if(d.pinceau.remplir)
          O1.ellipse_pleine(tl + Point{ep,ep}, br - Point{ep,ep});
      }
      else if((d.type == Impl::CanvaElement::MARQUEUR) && (pos0.y != -1))
      {
        cnt_marqueur++;
        marqueur(O1, d);
      }
    }

    DBG(msg("rvec:{}, rect:{}, lignes:{}, chaines:{}, cercles:{}, marqueurs:{}.", cnt_rvec, cnt_rec, cnt_ligne, cnt_chaine, cnt_cercle, cnt_marqueur));

    return O1;
  }
};



Rectf Canva::calc_rdi_englobante() const
{
  return impl->calc_rdi_englobante();
}


sptr<Canva::PointIntermediaire> Canva::enregistre_pti()
{
  auto res = std::make_shared<Canva::PointIntermediaire>();
  Impl::CanvaElement dess;
  dess.type = Impl::CanvaElement::PTI;
  dess.pti  = res;
  impl->elements.push_back(dess);
  return res;
}

void Canva::restaure_pti(sptr<Canva::PointIntermediaire> pti)
{
  tsd_assert(pti);
  // Efface tous les éléments actuels
  clear();
  // TODO : c'est du bricolage
  auto dim = get_allocation();
  this->dessine_img(0, 0, dim.l-1, dim.h-1, pti->img);
}

Image Canva::rendre(const Dim &dim, const Rectf &rdi, const Couleur &arp)
{
  return impl->rendre(dim, rdi, arp);
}

void Canva::set_dim_fonte(float echelle)
{
  impl->pinceau.dim_fonte = echelle;
}

void Canva::set_rdi(const Rectf &rdi)
{
  impl->clip.rdi    = rdi;
}

void Canva::set_allocation(const Dim &dim)
{
  impl->clip.target = Rect{0, 0, dim.l, dim.h};
  impl->clip.rdi    = Rect{0, 0, dim.l, dim.h};
}

Rectf Canva::get_rdi() const
{
  return impl->clip.rdi;
}

Dim Canva::get_allocation() const
{
  return impl->get_allocation();
}

void Canva::set_orientation(Orientation orient)
{
  impl->pinceau.orient = orient;
}

void Canva::set_epaisseur(int ep)
{
  impl->pinceau.epaisseur = ep;
}

void Canva::set_dotted(bool dotted)
{
  impl->pinceau.dotted = dotted;
}

void Canva::set_remplissage(bool remplir, const Couleur &coul)
{
  impl->pinceau.remplir = remplir;
  impl->pinceau.couleur_remplissage = coul;
}

Canva::Canva()
{
  impl = std::make_shared<Impl>();
}

void Canva::forward(Canva dest) const
{
  for(auto d: impl->elements)
    dest.impl->ajoute_element(d);
}

void Canva::clear()
{
  DBG(msg("axes : clear elements ({})...", impl->elements.size());)
  impl->elements.clear();
  DBG(msg("ok.");)
}

void Canva::set_align(Align hor, Align vert)
{
  impl->pinceau.alignement_horizontal = hor;
  impl->pinceau.alignement_vertical   = vert;
}


void Canva::set_couleur(const Couleur &coul)
{
  impl->pinceau.couleur_trait = coul;
}

void Canva::fleche(const Pointf &p0, const Pointf &p1)
{
  impl->fleche(p0,p1, impl->pinceau);
}

void Canva::ligne(const Pointf &p0, const Pointf &p1)
{
  impl->ligne(p0,p1, impl->pinceau);
}

void Canva::ligne(float x0, float y0, float x1, float y1)
{
  impl->ligne({x0,y0},{x1,y1}, impl->pinceau);
}

void Canva::marqueur(const Pointf &p, Marqueur m, float dim_pixels)
{
  impl->marqueur(p, m, dim_pixels, impl->pinceau);
}

void Canva::rectangle(float x0, float y0, float x1, float y1)
{
  impl->rectangle({x0, y0}, {x1, y1}, impl->pinceau);
}

void Canva::rectangle(const Pointf &p0, const Pointf &p1)
{
  rectangle(Rectf::englobant(p0, p1));
}

void Canva::rectangle(const Rectf &r)
{
  impl->rectangle(r.tl(), r.br(), impl->pinceau);
}

void Canva::dessine_img(float x0, float y0, float x1, float y1, Image img)
{
  impl->dessine_img({x0, y0}, {x1, y1}, img, impl->pinceau);
}

void Canva::dessine_accu(const ArrayXcf &pts)
{
  impl->dessine_accu(pts, impl->pinceau);
}



void Canva::remplissage_vertical(float x0, float y0, float x1, float y1, float y2, float y3)
{
  impl->remplissage_vertical({x0, y0}, {x1, y1}, y2, y3, impl->pinceau);
}

void Canva::def_image_fond(Image img)
{
  impl->img_fond = img;
}

/*void Canva::set_coord_mode_pixels(bool pixels)
{
  impl->pinceau.coord_mode_pixels = pixels;
}*/

void Canva::texte(const Pointf &p,
    const std::string &texte,
    const Dimf &dim_max)
{
  impl->texte(p, texte, dim_max, impl->pinceau);
}

void Canva::texte(float x0, float y0, const std::string &texte, float dx_max, float dy_max)
{
  this->texte(Pointf{x0, y0}, texte, Dimf{dx_max, dy_max});
}


void Canva::ellipse(const Pointf &p, float a, float b, float theta)
{
  impl->ellipse(p, a, b, theta, impl->pinceau);
}

void Canva::cercle(float x0, float y0, float r)
{
  impl->cercle({x0, y0}, r, impl->pinceau);
}

Pointf Canva::v2c(const Pointf &p) const
{
  return impl->clip.v2c(p);
}

Canva Canva::clip_alt(const Rectf &rdi, float xmin, float xmax, float ymin, float ymax, bool log_x, bool log_y) const
{
  Canva res;
  Rectf vue{xmin, ymin, xmax - xmin, ymax - ymin};
  res.impl->clip = impl->clip.clip(rdi, vue);

  if(log_x)
  {
    res.impl->clip.rdi.x = log10(xmin);
    res.impl->clip.rdi.l = log10(xmax) - log10(xmin);
  }
  if(log_y)
  {
    res.impl->clip.rdi.y = log10(ymin);
    res.impl->clip.rdi.h = log10(ymax) - log10(ymin);
  }

  res.impl->clip.log_x = log_x;
  res.impl->clip.log_y = log_y;
  res.impl->parent          = impl;
  return res;
}

Canva Canva::clip(const Rectf &rdi, const Rectf &vue, bool log_x, bool log_y) const
{
  Canva res;
  res.impl->clip = impl->clip.clip(rdi, vue);

  if(log_x)
  {
    res.impl->clip.rdi.x = log10(vue.x);
    res.impl->clip.rdi.l = log10(vue.x + vue.l) - log10(vue.x);
  }
  if(log_y)
  {
    res.impl->clip.rdi.y = log10(vue.y);
    res.impl->clip.rdi.h = log10(vue.y + vue.h) - log10(vue.y);
  }

  res.impl->clip.log_x = log_x;
  res.impl->clip.log_y = log_y;
  res.impl->parent          = impl;
  return res;
}

Canva Canva::vue(const Rectf &rdi) const
{
  Canva res;
  res.impl->clip = impl->clip.vue(rdi);
  res.impl->parent    = impl;
  return res;
}

Pointf Canva::coord_pixel(const Point &p)
{
  return impl->clip.c2v(p);
}


}
