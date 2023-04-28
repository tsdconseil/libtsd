#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"
#include <deque>

#define DBG(AA)


using namespace std;

namespace tsd::vue
{




//struct Canva::Impl::CanvaElement;

struct Canva::GroupeTextes
{
  bool done = non;
  float echelle = 0;
  vector</*Canva::Impl::CanvaElement*/void *> elements;
};

struct Canva::PointIntermediaire
{
  Image img;
};

struct Canva::Impl
{

  vector<sptr<GroupeTextes>> grp;

  struct Pinceau
  {
    Align alignement_vertical   = Align::DEBUT,
          alignement_horizontal = Align::DEBUT;
    Couleur couleur_trait{0, 0, 0},
            couleur_remplissage{255,255,255};
    entier epaisseur = 1;
    bouléen remplir  = non;
    bouléen dessiner_contours = oui;
    Orientation orient = Orientation::HORIZONTALE;
    // Une chaine de caractère : fonte précisée, ou bien,
    // un rectangle englobant
    float dim_fonte         = -1;
    bouléen dotted             = non;
    float alpha             = 1;


    bouléen texte_arrière_plan_actif = non;
    Couleur texte_arrière_couleur = Couleur::Blanc;
  };


  // Permet de passer d'un système de coordonnées à un autre
  struct ClipGen
  {
    // Position dans le canva parent
    Rectf target{0,0,0,0};
    // RDI = système de coordonnée virtuel
    Rectf rdi{0,0,1,1};

    bouléen log_x = non, log_y = non;

    inline float v2c_x(float x) const
    {
      si(log_x)
        x = log10(x);

      si(rdi.l == 0)
      {
        msg_avert("ClipGen::v2c_x : rdi.l = 0.");
        retourne target.x;
      }

      retourne target.x + target.l * (x - rdi.x) / rdi.l;
    }

    inline float v2c_y(float y) const
    {
      si(log_y)
        y = log10(y);

      si(rdi.h == 0)
      {
        msg_avert("ClipGen::v2c_x : rdi.h = 0.");
        retourne target.y;
      }

      retourne target.y + target.h * (y - rdi.y) / rdi.h;
    }

    Pointf c2v(const Pointf &p) const
    {
      si((target.l == 0) || (target.h == 0))
      {
        msg_avert("ClipGen::c2v : target = {}.", target);
        retourne {rdi.x,rdi.y};
      }


      retourne {(p.x - target.x) * rdi.l / target.l + rdi.x,
              (p.y - target.y) * rdi.h / target.h + rdi.y};
    }

    // Système interne ("virtuel") vers externe ("concret")
    Pointf v2c(const Pointf &p) const
    {
      retourne {v2c_x(p.x), v2c_y(p.y)};
    }

    Dimf v2c_l(const Dimf &d) const
    {
      retourne {v2c_lx(d.l), v2c_ly(d.h)};
    }

    Rectf v2c(const Rectf &r) const
    {
      retourne {v2c_x(r.x), v2c_y(r.y), v2c_lx(r.l), v2c_ly(r.h)};
    }

    // Attention : n'a pas de sens en échelle logarithmique
    float v2c_lx(float largeur) const
    {
      retourne largeur * target.l / rdi.l;
    }

    float v2c_ly(float hauteur) const
    {
      retourne hauteur * target.h / rdi.h;
    }

    ClipGen vue(const Rectf &rdi_) const
    {
      ClipGen res;
      res.target = this->rdi;
      res.rdi    = rdi_;
      retourne res;
    }

    ClipGen clip(const Rectf &rdi_, const Rectf &vue) const
    {
      ClipGen res;
      res.target = rdi_;
      res.rdi    = vue;
      retourne res;
    }
  };


  struct CanvaElement
  {
    enum Type
    {
      RECT, LIGNE, FLECHE, CHAINE, CERCLE, MARQUEUR, RVEC, ACCU, IMAGE, PTI, ELLIPSE
    } type = MARQUEUR;
    string chaine;
    sptr<GroupeTextes> grp;
    Pointf p0, p1;
    Dimf dim;
    float r = 3, y2 = 0, y3 = 0;
    float a, b;
    //entier dim_pixels = 3;
    entier dim_pixels = 5;


    Marqueur marqueur = Marqueur::CARRE;
    Pinceau pinceau;
    Veccf accu_points;
    Image img;
    sptr<PointIntermediaire> pti;


    Rectf get_rdi() const
    {
      switch(type)
      {
      case RECT:
      case LIGNE:
      case FLECHE:
        retourne Rectf::englobant(p0, p1);
      case IMAGE:
      {
        // Pb si coordonnées pixel !!!!
        soit res = Rectf::englobant(p0, p1);// + Pointf(1,1));

        //msg("get_rdi(image) : res = {}", res);
        // Bricolage !
        res.l++;
        res.h++;

        retourne res;
      }
      case CHAINE:
      {
        retourne Rectf(p0, Dimf{1.0f,1.0f}); // TODO
      }
      case MARQUEUR:
        retourne Rectf(p0, Dimf{0.0f,0.0f}); // TODO
      case CERCLE:
        retourne Rectf(p0 - Pointf(r/2,r/2), Dimf{r/2,r/2});
      default:
        retourne Rectf(p0, Dimf{1.0f,1.0f}); // TODO
      }
    }

  };


  std::weak_ptr<Impl> parent;
  ClipGen clip;

  Pinceau pinceau;
  std::deque<CanvaElement> elements;
  Image img_fond;


  Dim get_allocation() const
  {
    soit p = parent.lock();
    si(p)
    {
      soit dp = p->get_allocation();
      Dim res;
      Rectf rect_parent = p->clip.rdi;
      res.l = dp.l * (clip.target.l / rect_parent.l);
      res.h = dp.h * (clip.target.h / rect_parent.h);
      retourne res;
    }
    retourne {(entier) clip.target.l, (entier) clip.target.h};
  }

  Rectf get_rdi() const
  {
    retourne clip.rdi;
  }

  float pixels_par_lsb_x() const
  {
    retourne clip.v2c_lx(1.0f);
  }
  float pixels_par_lsb_y() const
  {
    retourne clip.v2c_ly(1.0f);
  }

  // pixel -> pos
  Pointf pixel2pos(entier x, entier y)
  {
    retourne clip.c2v(Pointf{x,y});
  }

  // pos -> pixel
  Point pos(const Pointf &p) const
  {
    soit res = clip.v2c(p);
    retourne {(entier) res.x, (entier) res.y};
  }

  Pointf posf(const Pointf &p) const
  {
    retourne clip.v2c(p);
  }

  void marqueur(Image O, CanvaElement &d)
  {
    soit ep = d.pinceau.epaisseur;
    soit pos0 = pos(d.p0);

    O.def_epaisseur(ep);
    O.def_couleur_dessin(d.pinceau.couleur_trait);
    O.def_couleur_remplissage(d.pinceau.couleur_remplissage);


    //DBG_AXES(msg("marqueur r={}...", d.r);)
    entier r2 = (entier) round(d.dim_pixels / sqrt(2.0f));
    si(d.marqueur == Marqueur::CARRE)
    {
      //msg("Marqueur carré : r = {}, r2 = {}", d.r, r2);
      Point p0{(entier) (pos0.x-r2),(entier) (pos0.y-r2)},
            p1{(entier) (pos0.x+r2),(entier) (pos0.y+r2)};
      O.rectangle_plein(p0, p1);
      O.rectangle(p0, p1);
    }
    sinon si(d.marqueur == Marqueur::ETOILE)
    {
      O.ligne_aa(pos0.x-r2, pos0.y-r2, pos0.x+r2, pos0.y+r2);
      O.ligne_aa(pos0.x-r2, pos0.y, pos0.x+r2, pos0.y);
      O.ligne_aa(pos0.x, pos0.y-r2, pos0.x, pos0.y+r2);
      O.ligne_aa(pos0.x-r2, pos0.y+r2, pos0.x+r2, pos0.y-r2);
    }
    sinon si(d.marqueur == Marqueur::CROIX)
    {
      O.ligne_aa(pos0.x-r2, pos0.y, pos0.x+r2, pos0.y);
      O.ligne_aa(pos0.x, pos0.y-r2, pos0.x, pos0.y+r2);
    }
    sinon si(d.marqueur == Marqueur::DIAMANT)
    {
      O.ligne_aa(pos0.x-r2, pos0.y, pos0.x, pos0.y-r2);
      O.ligne_aa(pos0.x, pos0.y-r2, pos0.x+r2, pos0.y);
      O.ligne_aa(pos0.x+r2, pos0.y, pos0.x, pos0.y+r2);
      O.ligne_aa(pos0.x, pos0.y+r2, pos0.x-r2, pos0.y);
    }
    sinon si(d.marqueur == Marqueur::POINT)
    {
      O.point(pos0);
    }
    sinon si(d.marqueur == Marqueur::CERCLE)
    {
      //msg("Marqueur cercle : r = {}", d.r);
      O.cercle_plein(pos0, d.dim_pixels);
      O.cercle(pos0, d.dim_pixels);
    }
    sinon
    {
      msg_erreur("Marqueur de type inconnu ({})", (entier) d.type);
    }
  }

  static Rectf union_rdi(const Rectf &r1, const Rectf &r2)
  {
    Rectf r;
    r.x = min(r1.x, r2.x);
    r.y = min(r1.y, r2.y);
    r.l = max(r1.x + r1.l, r2.x + r2.l) - r.x;
    r.h = max(r1.y + r1.h, r2.y + r2.h) - r.y;
    retourne r;
  }

  Rectf calc_rdi_englobante() const
  {
    Rectf rdi;
    pour(auto &d: elements)
    {
      Rectf r = d.get_rdi();


      si(std::isinf(r.l) || (std::isinf(r.h)))
      {
        echec("calc_rdi_englobante(): r={}, type={}, p0={}, p1={}", r, (entier) d.type, d.p0, d.p1);
      }

      si(&d == &(elements.front()))
        rdi = r;
      sinon
        rdi = union_rdi(rdi, r);

      //msg(" -> r={}, rdi={}", r, rdi);
    }
    //msg("Calc RDI englobante : {}", rdi);
    retourne rdi;
  }


  static CanvaElement clip_elmt(const CanvaElement &elt, const ClipGen &clip)
  {
    CanvaElement elt2 = elt;
    elt2.p0   = clip.v2c(elt.p0);

    elt2.p1   = clip.v2c(elt.p1);
    elt2.y2   = clip.v2c_y(elt.y2);
    elt2.y3   = clip.v2c_y(elt.y3);
    si(elt.dim.l > 0)
      elt2.dim  = clip.v2c_l(elt.dim);
    elt2.a    = clip.v2c_lx(elt.a);
    elt2.b    = clip.v2c_ly(elt.b);

    // si commenté : pb cercles
    // si décommenté : pb marqueurs
    elt2.r    = clip.v2c_lx(elt.r);

    pour(auto i = 0; i < elt.accu_points.rows(); i++)
    {
      elt2.accu_points(i).real(clip.v2c_x(elt.accu_points(i).real()));
      elt2.accu_points(i).imag(clip.v2c_y(elt.accu_points(i).imag()));
    }
    retourne elt2;
  }


  void ajoute_element(const CanvaElement &elt)
  {
    soit p = parent.lock();
    si(p)
    {
      soit elt2 = clip_elmt(elt, clip);
      p->ajoute_element(elt2);
      retourne;
    }
    elements.push_back(elt);
    if(elt.grp)
    {
      elt.grp->elements.push_back(&(elements.back()));
      grp.push_back(elt.grp);
    }
  }

  void ellipse(const Pointf &p, float a, float b, float theta, const Pinceau &pinceau)
  {
    Impl::CanvaElement d;
    d.p0      = p;
    d.a       = a;
    d.b       = b;
    d.type    = Impl::CanvaElement::ELLIPSE;
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
    si(std::isinf(p0.y) || std::isinf(p1.y))
      retourne;
    Impl::CanvaElement d;
    d.p0      = p0;
    d.p1      = p1;
    d.type    = Impl::CanvaElement::LIGNE;
    d.pinceau = pinceau;
    ajoute_element(d);
  }

  void fleche(const Pointf &p0, const Pointf &p1, const Pinceau &pinceau, float dim_pixel = 5)
  {
    si(std::isinf(p0.y) || std::isinf(p1.y))
      retourne;
    Impl::CanvaElement d;
    d.p0          = p0;
    d.p1          = p1;
    d.type        = Impl::CanvaElement::FLECHE;
    d.pinceau     = pinceau;
    d.dim_pixels  = dim_pixel;
    ajoute_element(d);
  }

  void dessine_accu(const Veccf &x, const Pinceau &pinceau)
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
    si(std::isinf(p1.y) || std::isinf(y2) || std::isinf(y3))
      retourne;

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



  void texte(const Pointf &p, cstring texte, const Dimf &dim_max, const Pinceau &pinceau, sptr<GroupeTextes> grp)
  {
    Impl::CanvaElement d;
    d.p0      = p;
    d.dim     = dim_max;
    d.type    = Impl::CanvaElement::CHAINE;
    d.pinceau = pinceau;
    d.chaine  = texte;
    d.grp     = grp;
    ajoute_element(d);
  }

  Image rendre(const Dim &dim_, const Rectf &rdi_, const Couleur &arp_)
  {
    for(auto &g: grp)
      g->done = non;


    Dim dim = dim_;

    si((dim.l < 0) && (clip.target.l > 0))
    {
      dim.l = clip.target.l;
      dim.h = clip.target.h;
    }


    si(dim.l < 0)
    {
      dim.l = 800;
      dim.h = 600;
    }

    Rectf rdi = rdi_;

    /*si((rdi.l == 0) && (rdi.h == 0))
    {
      rdi = clip.rdi;
    }*/


    si((rdi.l == 0) && (rdi.h == 0))
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

    si(!img_fond.empty())
      O1 = img_fond.redim(dim);
    sinon
    {
      O1.resize(dim.l, dim.h);
      O1.remplir(arp);
    }

    DBG(msg("Dessins ({})...", elements.size());)

    entier cnt_rvec = 0, cnt_rec = 0, cnt_ligne = 0, cnt_chaine = 0, cnt_cercle = 0, cnt_marqueur = 0;

    pour(auto &d: elements)
    {

      soit ep = d.pinceau.epaisseur;

      soit pos0 = pos(d.p0);
      soit pos1 = pos(d.p1);

      O1.def_epaisseur(ep);
      O1.def_couleur_dessin(d.pinceau.couleur_trait);
      O1.def_couleur_remplissage(d.pinceau.couleur_remplissage);

      // Sera possible à partir de GCC 11 🙂
      // using enum Impl::CanvaElement::Type;

      si(d.type == Impl::CanvaElement::PTI)
      {
        // Point intermédiaire
        // il faut sauvegarder l'état du dessin à ce moment là
        tsd_assert(d.pti);
        d.pti->img = O1.clone();
      }
      sinon si(d.type == Impl::CanvaElement::ACCU)
      {
        //msg("Rendu accu : (1)");
        DBG(msg("accu...");)
        soit c = d.pinceau.couleur_trait;
        soit n = d.accu_points.rows();
        soit dim = get_allocation();
        DBG(msg("axes2: rendu accu, allocation = {}, npts = {}", dim, n));
        soit a = Tabf::zeros(dim.l, dim.h);
        pour(auto i = 0; i < n; i++)
        {
          Point p = pos({d.accu_points(i).real(), d.accu_points(i).imag()});
          entier r = 1;
          si((p.x >= r) && (p.y >= r) && (p.x < dim.l-r) && (p.y < dim.h-r))
            a.block(p.x-r, p.y-r, 2*r+1, 2*r+1) += 1;
        }

        //msg("Rendu accu : (2)");

        soit mc = a.valeur_max();
        DBG(msg("axes2: rendu accu, max = {}", mc));
        a *= 1.0f / (mc + 1e-20);
        pour(auto y = 0; y < dim.h; y++)
        {
          pour(auto x = 0; x < dim.l; x++)
          {
            si(a(x,y) > 0)
              O1.point({x,y}, c, a(x,y));
          }
        }
        //msg("Rendu accu : (fin)");
      }
      sinon si(d.type == Impl::CanvaElement::IMAGE)
      {
        DBG(msg("image...");)
        soit r = Rect::englobant(pos0, pos1);

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
        soit i = d.img;

        // r : unités pixels

        // Ne prends en compte que ce qui est à l'intérieur de l'image
        si((r.x < 0) || (r.y < 0) || (r.x + r.l > O1.sx()) || (r.y + r.h > O1.sy()))
        {
          msg("CUT! (r = {}, pos0 = {}, pos1 = {}, d.p0={}, d.p1={})", r, pos0, pos1, d.p0, d.p1);
          si(r.x < 0)
          {
            // Coupe la partie gauche de l'image
            entier cut = (-i.sx() * r.x) / r.l;
            i = i.sous_image(cut, 0, i.sx() - cut, i.sy());
            r.l += r.x;
            r.x = 0;
          }
          si(r.x + r.l > O1.sx())
          {
            // Coupe la partie droite de l'image
            entier cut = (i.sx() * (r.x+r.l-O1.sx())) / r.l;
            i = i.sous_image(0, 0, i.sx() - cut, i.sy());
            r.l = O1.sx()-1-r.x;
          }
          si(r.y < 0)
          {
            // Coupe la partie haute de l'image
            entier cut = (-i.sy() * r.y) / r.h;
            i = i.sous_image(0, cut, i.sx(), i.sy() - cut);
            r.h += r.y;
            r.y = 0;
          }
          si(r.y + r.h > O1.sy())
          {
            // Coupe la partie basse de l'image
            entier cut = (i.sy() * (r.y+r.h-O1.sy())) / r.h;
            i = i.sous_image(0, 0, i.sx(), i.sy() - cut);
            r.h = O1.sy()-1-r.y;
          }

          msg("Après cut: rect = {}, dim img = {}, dim O1={}",
              r, i.get_dim(), O1.get_dim());
        }

        //msg("p0={}, p1={}", d.p0, d.p1);
        O1.puti(r, i);

        DBG(msg("image ok.");)
      }
      sinon si(d.type == Impl::CanvaElement::RECT)
      {
        cnt_rec++;
        DBG(msg("rect...");)
        soit r = Rect::englobant(pos0, pos1);
        soit p0 = r.tl();
        soit p1 = r.br();

        /*static int cnt = 0;
        if(cnt++ < 100)
        {
          msg("p0 = {}, p1 = {}", p0, p1);
        }*/

        O1.def_epaisseur(1);

        //O1.rectangle({p0.x,p0.y}, {p1.x,p1.y});
        //O1.rectangle_plein({p0.x, p0.y},{p1.x,p1.y});

        si(d.pinceau.dessiner_contours)
          pour(auto k = 0; k < ep; k++)
            O1.rectangle({p0.x+k,p0.y+k}, {p1.x-k,p1.y-k});

        si(d.pinceau.remplir && (p0.x != p1.x) && (p0.y != p1.y))
        {
          soit e = d.pinceau.dessiner_contours ? ep : 0;
          O1.rectangle_plein({p0.x+e, p0.y+e},{p1.x-e,p1.y-e});
        }

        O1.def_epaisseur(ep);

        DBG(msg("rect ok.");)
      }
      sinon si(d.type == Impl::CanvaElement::RVEC)
      {
        DBG(msg("rvec...");)
        cnt_rvec++;
        soit pos2 = pos({d.p0.x, d.y2});
        soit pos3 = pos({d.p1.x, d.y3});

        si((pos2.y != -1) && (pos3.y != -1))
        {
          pour(auto x = pos0.x; x <= pos1.x; x++)
          {
            // Interpolation linéaire entre les points
            entier yh = pos0.y + ((pos1.y - pos0.y) * (x - pos0.x) * 1.0f) / (pos1.x - pos0.x);
            entier yb = pos2.y + ((pos3.y - pos2.y) * (x - pos0.x) * 1.0f) / (pos1.x - pos0.x);
            O1.ligne({x, yh}, {x, yb});
          }
        }
        DBG(msg("rvec ok.");)
      }
      sinon si(d.type == Impl::CanvaElement::LIGNE)
      {

        cnt_ligne++;
        soit [x0,y0] = posf(d.p0);
        soit [x1,y1] = posf(d.p1);
        DBG(msg("ligne ({},{}) - ({},{})...", x0, y0, x1, y1);)
        O1.ligne_aa(x0, y0, x1, y1, d.pinceau.dotted ? Image::POINTILLEE : Image::PLEINE);
        DBG(msg("ligne ok.");)
      }
      sinon si(d.type == Impl::CanvaElement::FLECHE)
      {
        O1.fleche(pos0, pos1, d.pinceau.dotted ? Image::POINTILLEE : Image::PLEINE, d.dim_pixels);
      }
      sinon si(d.type == Impl::CanvaElement::CHAINE)
      {
        if(d.chaine.empty())
          continue;

        cnt_chaine++;
        DBG(msg("chaine '{}'...", d.chaine););

        float fonte = 1;

        si(d.pinceau.dim_fonte != -1)
          fonte = d.pinceau.dim_fonte;

        tsd::vue::TexteConfiguration config;

        si(d.pinceau.texte_arrière_plan_actif)
          config.couleur_fond = d.pinceau.texte_arrière_couleur;
        config.couleur    = d.pinceau.couleur_trait;
        config.scale      = fonte;
        config.thickness  = d.pinceau.epaisseur;


        entier m = 0;

        si(d.pinceau.texte_arrière_plan_actif)
          m = 2;

        config.dim_max    = Dim{O1.sx() - pos0.x - 2*m, O1.sy() - pos0.y - 2*m};


        si(d.dim.l > 0)
          config.dim_max.l = min((entier) clip.v2c_lx(d.dim.l), config.dim_max.l);
        si(d.dim.h > 0)
          config.dim_max.h = min((entier) clip.v2c_ly(d.dim.h), config.dim_max.h);

        si(d.pinceau.orient == Orientation::VERTICALE)
        {
          std::swap(config.dim_max.l, config.dim_max.h);
        }

        DBG(msg("  creation image, echelle = {}, dmax = {}, dimax_in = {}...", fonte, config.dim_max, d.dim););

        // Si appartient à un groupe
        si(d.grp)
        {
          si(!d.grp->done)
          {
            DBG(msg("Canva: groupe de texte commun détecté.");)
            float s = 0;
            for(auto elptr: d.grp->elements)
            {
              auto el = (CanvaElement *) elptr;
              TexteProps props;
              soit i = tsd::vue::texte_creation_image(el->chaine, config, &props);
              if(s == 0)
                s = props.scale_out;
              else
                s = min(s, props.scale_out);
            }
            d.grp->done = oui;
            d.grp->echelle = s;
            DBG(msg("Canva/groupe de texte: échelle : {} -> {}", config.scale, s);)
          }

          config.scale = d.grp->echelle;
        }


        soit i = tsd::vue::texte_creation_image(d.chaine, config);

        DBG(msg("  ok : {}...", i.get_dim()););
        si(i.empty())
        {
          DBG(msg("  (empty)"));
          continue;
        }



        si(d.pinceau.texte_arrière_plan_actif)
        {
          // Ajoute une marge
          Image pg(i.sx() + 2 * m, i.sy() + 2 * m);
          pg.remplir(d.pinceau.texte_arrière_couleur);
          pg.puti({m,  m}, i);
          i = pg;
        }

        si((pos0.x < 0) || (pos0.y < 0))
        {
          //msg_avert("Axes : affichage texte [{}], position invalide {}, position flt = {}.", d.chaine, pos0, d.p0);
          //msg("mode pixels = {}", d.pinceau.coord_mode_pixels);
        }

        si(d.pinceau.orient == Orientation::VERTICALE)
          i = i.rotation_90();

        //O1.puti_avec_gamma(pos0, i);

        entier px, py;

        si(d.pinceau.alignement_horizontal == Align::DEBUT)
          px = pos0.x;
        sinon si(d.pinceau.alignement_horizontal == Align::FIN)
          px = pos0.x + d.p1.x - i.sx();
        sinon // Centre
          px = pos0.x - i.sx() / 2;


        si(d.pinceau.alignement_vertical == Align::DEBUT)
          py = pos0.y;
        sinon si(d.pinceau.alignement_vertical == Align::FIN)
          py = pos0.y + d.p1.y - i.sy();
        sinon // Centre
          py = pos0.y - i.sy() / 2 - 1;

        DBG(msg("  puti : @{}x{}, dim = {}.", px, py, i.get_dim()));

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
        DBG(msg("  chaine: ok."));
      }
      sinon si(d.type == Impl::CanvaElement::CERCLE)
      {
        cnt_cercle++;
        DBG(msg("cercle r={}...", d.r);)
        entier dx = d.r * pixels_par_lsb_x(),
            dy = d.r * pixels_par_lsb_y();
        DBG(msg("  -> dx={} pixels, dy = {} pixels", dx, dy);)

        Point tl{pos0.x-dx,pos0.y-dy}, br{pos0.x+dx,pos0.y+dy};

        si(d.pinceau.dessiner_contours)
          O1.ellipse(tl, br);

        si(d.pinceau.remplir)
        {
          soit e = d.pinceau.dessiner_contours ? ep : 0;
          Point δ{e, e};
          O1.ellipse_pleine(tl + δ, br - δ);
        }
        DBG(msg("cercle ok.");)
      }
      sinon si(d.type == Impl::CanvaElement::ELLIPSE)
      {
        cnt_cercle++;
        DBG(msg("ellipse a={}, b={}...", d.a, d.b);)
        entier dx = d.a * pixels_par_lsb_x(),
               dy = d.b * pixels_par_lsb_y();
        DBG(msg("  -> dx={} pixels, dy = {} pixels", dx, dy);)

        Point tl{pos0.x-dx,pos0.y-dy}, br{pos0.x+dx,pos0.y+dy};



        si(d.pinceau.remplir)
        {
          soit e = d.pinceau.dessiner_contours ? ep : 0;
          Point δ{e, e};
          si(abs(dx) == abs(dy))
            O1.cercle_plein(pos0, dx);
          sinon
            O1.ellipse_pleine(tl + δ, br - δ);
        }

        si(d.pinceau.dessiner_contours)
        {
          si(abs(dx) == abs(dy))
            O1.cercle(pos0, dx);
          sinon
            O1.ellipse(tl, br);
        }

        DBG(msg("ellipse ok.");)
      }
      sinon si((d.type == Impl::CanvaElement::MARQUEUR) && (pos0.y != -1))
      {
        DBG(msg("marqueur...");)
        cnt_marqueur++;
        marqueur(O1, d);
        DBG(msg("marqueur ok.");)
      }
    }

    DBG(msg("rvec:{}, rect:{}, lignes:{}, chaines:{}, cercles:{}, marqueurs:{}.", cnt_rvec, cnt_rec, cnt_ligne, cnt_chaine, cnt_cercle, cnt_marqueur));

    retourne O1;
  }
};



Rectf Canva::calc_rdi_englobante() const
{
  retourne impl->calc_rdi_englobante();
}


sptr<Canva::PointIntermediaire> Canva::enregistre_pti()
{
  soit res = std::make_shared<Canva::PointIntermediaire>();
  Impl::CanvaElement dess;
  dess.type = Impl::CanvaElement::PTI;
  dess.pti  = res;
  impl->elements.push_back(dess);
  retourne res;
}

sptr<Canva::GroupeTextes> Canva::groupe_textes()
{
  soit res = std::make_shared<Canva::GroupeTextes>();
  retourne res;
}

void Canva::restaure_pti(sptr<Canva::PointIntermediaire> pti)
{
  tsd_assert(pti);
  // Efface tous les éléments actuels
  clear();
  // TODO : c'est du bricolage
  soit dim = get_allocation();
  this->dessine_img(0, 0, dim.l-1, dim.h-1, pti->img);
}

Image Canva::rendre(const Dim &dim, const Rectf &rdi, const Couleur &arp)
{
  retourne impl->rendre(dim, rdi, arp);
}

void Canva::plot_cmap(const Tabf &Z, const Rectf &rdi, sptr<CMap> cmap)
{
  soit nr = Z.rows(), nc = Z.cols();
  ParamGrille pg(nr, nc, rdi.tl(), rdi.br());

  //DBG(msg("plot img : tl = {}, br = {}", ci.rdi_z.tl(), ci.rdi_z.br());)

  //canva.ligne(ci.rdi_z.tl(), ci.rdi_z.br());

  active_contours(non);
  active_remplissage(oui);

  pour(auto i = 0; i < nr; i++)
  {
    pour(auto j = 0; j < nc; j++)
    {
      soit coul = cmap->couleur(Z(i,j));

      //si(((i % 5) == 0) && ((j % 5) == 0))
      ///{
      //coul = Couleur::rand();
      //coul.alpha = 128;
      set_couleur_remplissage(coul);
      rectangle(pg.tl(i, j), pg.br(i, j));
      //}
    }
  }
  active_contours(oui);
  active_remplissage(non);
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
  retourne impl->clip.rdi;
}

Dim Canva::get_allocation() const
{
  retourne impl->get_allocation();
}

void Canva::set_orientation(Orientation orient)
{
  impl->pinceau.orient = orient;
}

void Canva::set_epaisseur(entier ep)
{
  impl->pinceau.epaisseur = ep;
}

void Canva::set_dotted(bouléen dotted)
{
  impl->pinceau.dotted = dotted;
}

void Canva::set_texte_arrière_plan(bouléen actif, const Couleur &coul)
{
  impl->pinceau.texte_arrière_plan_actif = actif;
  impl->pinceau.texte_arrière_couleur = coul;
}

void Canva::set_couleur_remplissage(const Couleur &coul)
{
  impl->pinceau.couleur_remplissage = coul;
}

void Canva::set_remplissage(bouléen remplir, const Couleur &coul)
{
  impl->pinceau.remplir = remplir;
  impl->pinceau.couleur_remplissage = coul;
}

void Canva::active_contours(bouléen contour_actif)
{
  impl->pinceau.dessiner_contours = contour_actif;
}

void Canva::active_remplissage(bouléen remplissage_actif)
{
  impl->pinceau.remplir = remplissage_actif;
}

Canva::Canva()
{
  impl = std::make_shared<Impl>();
}

void Canva::forward(Canva dest) const
{
  pour(auto d: impl->elements)
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

void Canva::fleche(const Pointf &p0, const Pointf &p1, float dim)
{
  impl->fleche(p0,p1, impl->pinceau, dim);
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

void Canva::dessine_accu(const Veccf &pts)
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

/*void Canva::set_coord_mode_pixels(bouléen pixels)
{
  impl->pinceau.coord_mode_pixels = pixels;
}*/

void Canva::texte(const Pointf &p,
    cstring texte,
    const Dimf &dim_max, sptr<GroupeTextes> grp)
{
  impl->texte(p, texte, dim_max, impl->pinceau, grp);
}

void Canva::texte(float x0, float y0, cstring texte, float dx_max, float dy_max, sptr<GroupeTextes> grp)
{
  this->texte(Pointf{x0, y0}, texte, Dimf{dx_max, dy_max}, grp);
}


void Canva::ellipse(const Pointf &p, float a, float b, float theta)
{
  impl->ellipse(p, a, b, theta, impl->pinceau);
}

void Canva::cercle(const Pointf &p, float r)
{
  impl->cercle(p, r, impl->pinceau);
}

void Canva::cercle(float x0, float y0, float r)
{
  impl->cercle({x0, y0}, r, impl->pinceau);
}

Pointf Canva::v2c(const Pointf &p) const
{
  retourne impl->clip.v2c(p);
}

Rectf Canva::v2c(const Rectf &r) const
{
  retourne impl->clip.v2c(r);
}

Canva Canva::clip_alt(const Rectf &rdi, float xmin, float xmax, float ymin, float ymax, bouléen log_x, bouléen log_y) const
{
  Canva res;
  Rectf vue{xmin, ymin, xmax - xmin, ymax - ymin};
  res.impl->clip = impl->clip.clip(rdi, vue);

  si(log_x)
  {
    res.impl->clip.rdi.x = log10(xmin);
    res.impl->clip.rdi.l = log10(xmax) - log10(xmin);
  }
  si(log_y)
  {
    res.impl->clip.rdi.y = log10(ymin);
    res.impl->clip.rdi.h = log10(ymax) - log10(ymin);
  }

  res.impl->clip.log_x = log_x;
  res.impl->clip.log_y = log_y;
  res.impl->parent          = impl;
  retourne res;
}

Canva Canva::clip(const Rectf &rdi, const Rectf &vue, bouléen log_x, bouléen log_y) const
{
  Canva res;
  res.impl->clip = impl->clip.clip(rdi, vue);

  si(log_x)
  {
    res.impl->clip.rdi.x = log10(vue.x);
    res.impl->clip.rdi.l = log10(vue.x + vue.l) - log10(vue.x);
  }
  si(log_y)
  {
    res.impl->clip.rdi.y = log10(vue.y);
    res.impl->clip.rdi.h = log10(vue.y + vue.h) - log10(vue.y);
  }

  res.impl->clip.log_x = log_x;
  res.impl->clip.log_y = log_y;
  res.impl->parent          = impl;
  retourne res;
}

Canva Canva::vue(const Rectf &rdi) const
{
  Canva res;
  res.impl->clip = impl->clip.vue(rdi);
  res.impl->parent    = impl;
  retourne res;
}

Pointf Canva::coord_pixel(const Point &p)
{
  retourne impl->clip.c2v(p);
}


}
