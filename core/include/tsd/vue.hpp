

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#ifndef TSD_VUE_HPP
#define TSD_VUE_HPP

#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"


namespace tsd::vue::unites {
extern tuple<entier,entier> calc_expo_nb_chiffres(double t, cstring unite);
extern string valeur_vers_chaine(double t, cstring unite, entier expo, entier nb_chiffres);
extern tuple<entier, entier> calc_expo_nb_chiffres_commun(const vector<double> &tics, cstring unite);
extern string valeur_vers_chaine(double t, cstring unite);
extern string valeur_vers_chaine(double t, cstring unite, entier nb_chiffres);
}



namespace tsd::vue
{

/** @addtogroup vue
 *  @{
 */

extern Couleur eclaircir_couleur(const Couleur &src, float ratio = 0.5);



struct ConfigGrille
{
  entier afficher = -1;
  entier lg_tiret = -1, lg_trou = -1;
  Couleur couleur{0,0,0,0};
  entier epaisseur = 1;

  // Par défaut : 50 pour majeur, 4 pour mineur
  entier pixels_par_tic_min_x = -1;
  entier pixels_par_tic_min_y = -1;

  entier sous_graduation = 5;

  void from_style(cstring style);
  string vers_style() const;

  std::strong_ordering operator <=>(const ConfigGrille &c) const = default;
};

/** @brief Configuration pour l'affichage d'un axe (horizontal ou vertical) */
struct ConfigAxe
{
  string id, label, unité;

  /** @brief Pour zoom manuel. */
  float vmin = -1.0, vmax = 1.0;

  bouléen afficher                   = oui;
  entier ecart_pixel                 = -1;
  float ecart_valeur              = -1;

  bouléen echelle_logarithmique      = non;

  /** @brief Si vrai, les tics affichés sont ceux précisés dans le vecteur "tics". */
  bouléen tic_manuels = non;

  vector<std::pair<float, string>> tics;

  /** @brief Dim fonte valeurs des tics. */
  float valeurs_echelle = 0.6;

  /** @brief Hauteur max. pour les valeurs des tics, en pixels. */
  entier valeurs_hauteur_max = 18;

  /** @brief Décalage avec l'axe, en pixels. */
  entier valeurs_décalage = 4;


  bouléen operator <=>(const ConfigAxe &c) const = default;

  operator string() const;
};

enum class Marqueur
{
  AUCUN, CARRE, DIAMANT, ETOILE, POINT, CERCLE, CROIX
};

enum class Trait
{
  AUCUN, LIGNE, DOTTED, BATON, HISTO, MINMAX
};

/** @brief Channel display configuration */
struct ConfigSerie
{
  string nom;
  Couleur couleur{0,0,0};
  entier largeur_de_trait  = 1;
  string unite;
  bouléen visible          = oui;
  bouléen remplissage      = non;
  string style;

  Marqueur marqueur = Marqueur::AUCUN;
  Trait trait = Trait::LIGNE;

  bouléen operator <=>(const ConfigSerie &c) const = default;
};



struct ConfigAxes
{
  void ajoute_serie(cstring nom);

  void def_rdi_visible(float xmin, float xmax, float ymin, float ymax);
  void def_rdi_visible_abs(float xmin, float xmax, float ymin, float ymax);

  /** @brief List of channels */
  vector<ConfigSerie> series;

  ConfigGrille grille_majeure, grille_mineure;
  ConfigAxe axe_vertical, axe_horizontal;

  string titre;

  float titre_echelle = 1.0;

  /** @brief Background color */
  Couleur couleur_arriere_plan = Couleur{255,255,255};
  Couleur couleur_avant_plan = Couleur{0,0,0};

  //bouléen affiche_cartouche_legende  = oui;

  bouléen est_mat_sur_clair() const;

  struct Legende
  {
    bouléen afficher = oui;

    struct Position
    {
      enum Type
      {
        NORD_EST = 0,
        SUD_EST,
        SUD_OUEST,
        NORD_OUEST,
        AUCUNE
      } type = NORD_EST;

      operator string() const;
      Position(cstring id = "ne");
      Position(const char *id): Position(string(id)){}

      std::strong_ordering operator <=>(const Position &c) const = default;
    } position;

    float dim = 0.5f;

    std::strong_ordering operator <=>(const Legende &c) const = default;

  } legende;



  //entier legende_dim = 6;

  enum PositionAxeVertical
  {
    AXE_VERTICAL_GAUCHE = 0,
    AXE_VERTICAL_DROITE
  } position_axe_vertical = AXE_VERTICAL_DROITE;

  bouléen mode_isoview = non;

  float aspect_ratio = 1.0f;

  bouléen operator <=>(const ConfigAxes &c) const = default;

  operator string() const;
};

struct Axe
{
  ConfigAxe config;

  entier dim; // nb pixels
  vector<double> tics_majeurs_pos;
  vector<double> tics_mineurs_pos;

  double ecart_tics_majeurs_valeur;
  entier   ecart_tics_majeurs_pixels;

  double ecart_tics_mineurs_valeur;
  entier   ecart_tics_mineurs_pixels;

  bouléen tics_mineurs_presents;
  /** Center of X / Y axis, in lsb units */
  float   decalage_zero;

  /** This defines the display gain, for Y axis */
  float   pixels_par_lsb;

  bouléen vertical = non;

  void configure(const ConfigAxe &config, entier dim);

  bouléen tic_ok(float v);
  bouléen tic_ok_souple(float v);

  float pixel_vers_valeur(float p) const;
  entier valeur_vers_pixel(float v) const;
  float valeur_vers_pixelf(float v) const;
  float vmin() const;
  float vmax() const;
};

enum class Align
{
  DEBUT, CENTRE, FIN
};


struct Canva
{
  Canva();

  // Dessine dans un sous-rectangle du Canva passé en paramètre
  // rdi est spécifiée en coordonnées virtuelle de l'ancien canva.
  // Exemples pour une grille :
  // Canva c;
  // c.clip({0,0,0.5,0.5});
  // figure1.rendre(c);
  //   -> Dans figure::rendre
  //       auto c1 = c.vue({xmin,ymin,xmax-xmin,ymax-ymin})
  //       c1.ligne(x(i), y(i), x(i+1), y(i+1))

  // clip permet de modifier la cible des dessins
  // vue  permet de modifier la source des dessins (système de coordonnées)

  // x = x_clip + rdi_clip.l * (xv - rdi_view.x) / rdi_view.l

  // soit :
  // x = alpha + beta * xv


  // Donc dans l'ordre :
  // - canva initial, de dimension e.g. 640 x 480
  // - vue initiale = aussi 640 x 480
  // - C'est à dire qu'initialement, x virtuel = x concret
  // - On pose alpha = 1, beta = 0 (coefs de conversion virtuel -> concret)
  // - Effet de clip :
  //    beta_x += rdi.x, alpha_x = inchangé
  // - Effet de view :
  //    beta_x += -rdi.x / rdi.l, alpha_x *= 1 / rdi.l

  // Dessine dans un sous-rectangle de la vue en cours
  Canva clip(const Rectf &rdi, const Rectf &vue, bouléen log_x = non, bouléen log_y = non) const;
  Canva clip_alt(const Rectf &rdi, float xmin, float xmax, float ymin, float ymax, bouléen log_x, bouléen log_y) const;

  // Change l'échelle de la vue en cours
  Canva vue(const Rectf &vue) const;

  Pointf v2c(const Pointf &p) const;
  Rectf v2c(const Rectf &r) const;

  void forward(Canva dest) const;


  void clear();

  // A VOIR
  void set_allocation(const Dim &dim);
  void set_rdi(const Rectf &rdi);
  Dim get_allocation() const;

  Rectf get_rdi() const;

  void set_orientation(Orientation orient);
  void set_dim_fonte(float echelle);
  void set_couleur(const Couleur &coul);
  void set_couleur_remplissage(const Couleur &coul);
  void set_epaisseur(entier ep);
  void set_remplissage(bouléen remplir, const Couleur &coul = Couleur{255,255,255});

  void active_contours(bouléen contour_actif);
  void active_remplissage(bouléen remplissage_actif);


  void set_texte_arrière_plan(bouléen actif, const Couleur &coul = Couleur::Blanc);

  void set_dotted(bouléen dotted);
  void set_align(Align hor, Align vert);
  void rectangle(float x0, float y0, float x1, float y1);
  void rectangle(const Rectf &r);
  void rectangle(const Pointf &p0, const Pointf &p1);
  void ligne(float x0, float y0, float x1, float y1);
  void ligne(const Pointf &p0, const Pointf &p1);
  void fleche(const Pointf &p0, const Pointf &p1, float dim = 5);

  struct GroupeTextes;
  void texte(float x0, float y0,
      cstring texte,
      float dx_max = -1, float dy_max = -1, sptr<GroupeTextes> grp = sptr<GroupeTextes>());

  void texte(const Pointf &p,
      cstring texte,
      const Dimf &dim_max = {-1.0f, -1.0f}
      , sptr<GroupeTextes> grp = sptr<GroupeTextes>());



  void cercle(float x0, float y0, float r);
  void cercle(const Pointf &p, float r);
  void ellipse(const Pointf &p, float a, float b, float theta = 0);
  void remplissage_vertical(float x0, float y0, float x1, float y1, float y2, float y3);
  void dessine_accu(const Veccf &pts);
  void def_image_fond(Image img);
  // TODO : décrire différence entre dessine_img et plot_cmap...
  void dessine_img(float x0, float y0, float y1, float y2, Image img);
  void plot_cmap(const Tabf &Z, const Rectf &rdi, sptr<CMap> cmap);

  void marqueur(const Pointf &p, Marqueur m, float dim_pixels);
  void get_pixel_size(float &width, float &heigth);
  /** @brief Coordonnées pixel vers unités */
  Pointf coord_pixel(const Point &p);


  // TODO : rdi nécessaire ? n'est-il pas déjà précisé dans le clip interne ?
  Image rendre(const Dim &dim = {-1, -1}, const Rectf &rdi = {0,0,0,0}, const Couleur &arp = Couleur::Blanc);
  struct PointIntermediaire;

  sptr<PointIntermediaire> enregistre_pti();
  void restaure_pti(sptr<Canva::PointIntermediaire> pti);


  /** Création d'un groupe de texte dont la fonte sera forcée à être identique */
  sptr<GroupeTextes> groupe_textes();

  Rectf calc_rdi_englobante() const;

  _PIMPL_
};


struct Rendable: std::enable_shared_from_this<Rendable>
{
  virtual void rendre(Canva canva) const = 0;

  /** @brief Enregistrement sous la forme d'un fichier image */
  void enregistrer(cstring chemin_fichier, const Dim &dim = {-1, -1}) const;

  virtual void afficher(cstring titre = "", const Dim &dim = {-1, -1}) const;

  Image genere_image(const Dim &dim = Dim{-1,-1}, const Couleur &arp = Couleur::Blanc) const;


  virtual Dim dimensions_idéales() const{retourne {-1,-1};}
};


struct ARendable
{
  virtual sptr<const Rendable> rendable() const = 0;

  void enregistrer(cstring chemin_fichier, const Dim &dim = {-1, -1}) const
  {
    rendable()->enregistrer(chemin_fichier, dim);
  }
  Image genere_image(const Dim &dim = Dim{-1,-1}, const Couleur &arp = Couleur::Blanc) const
  {
    return rendable()->genere_image(dim, arp);
  }
  void afficher(cstring titre = "", const Dim &dim = {-1, -1}) const
  {
    rendable()->afficher(titre, dim);
  }
  void rendre(Canva canva) const
  {
    rendable()->rendre(canva);
  }
};




struct Axes
{
  Axes();

  Rectf get_rdi() const;
  void def_rdi_visible(float xmin, float xmax, float ymin, float ymax);
  void def_rdi_visible_abs(float xmin, float xmax, float ymin, float ymax);

  void set_isoview(bouléen isoview);
  void fixe_aspect_ratio(float ysx);

  void active_grilles(bouléen actives);
  void supprime_decorations();

  /** @brief Définition de l'échelle d'affichage (lin ou log)
   *  @todo : n'a pas l'air de fonctionner ? */
  void def_echelle(cstring x, cstring y);


  Canva rendre(Canva canva_, float xmin, float xmax, float ymin, float ymax) const;
  Canva rendre(Canva canva, const Rectf &rdi) const;
  void post_rendu(Canva canva_) const;



  Dim get_dim_graphique(const Dim &dim_totale) const;
  void configure(const ConfigAxes &config);
  ConfigAxes &get_config();

  Axes clone();

  _PIMPL_
};




/** @brief Système de coordonnées pour l'affichage d'une surface 2d. */
struct ParamGrille
{
  entier nc = 0, nr = 0;
  float x0 = 0, y0 = 0;
  float ẟx = 0, ẟy = 0;

  ParamGrille(entier nc, entier nr, const Pointf &p0, const Pointf &p1);
  Dimf dim() const;
  Pointf ctr(entier ix, entier iy) const;
  Pointf tl(entier ix, entier iy) const;
  Pointf br(entier ix, entier iy) const;
  Pointf bl(entier ix, entier iy) const;
  Pointf tr(entier ix, entier iy) const;
};




/** @brief %Figure avec API de type Matlab / Scilab.
 *
 *  <h3>%Figure avec API de type Matlab / Scilab</h3>
 *
 *  @par Exemple, avec différentes couleurs et différents marqueurs :
 *  @snippet exemples/src/ex-vue.cc ex_figure
 *  @image html figure.png width=800px
 */
struct Figure: ARendable
{
  /** @brief %Courbe (telle que renvoyée par la méthode plot()).
   *
   * Cette classe permet de modifier a posteriori les propriétés d'un courbe (épaisseur de trait, couleur, etc.)
   *
   */
  struct Courbe
  {
    /** @brief Changement de couleur */
    void def_couleur(const tsd::vue::Couleur &coul);

    /** @brief Définit une couleur différente pour chaque point de la courbe */
    void def_couleurs(const Vecf &c, cstring cmap = "");

    /** @brief Changement d'épaisseur du trait */
    void def_epaisseur(entier ep);

    /** @brief Remplissage sous la courbe */
    void def_remplissage(bouléen actif, bouléen vers_vmin = non, float vmin = 0);

    /** @brief Définit, pour chaque point de la courbe, l'écart-type */
    void def_σ(const Vecf &σ);

    /** @brief Définit la dimension (en pixels) des marqueurs. */
    void def_dim_marqueur(entier dim);

    void def_légende(cstring titre);

  //private:
    struct Impl;
    sptr<Impl> impl;
  };

  Figure(cstring nom = "");

  /** @brief Efface la figure */
  void clear();



  // void def_couleur(float r, float g, float b);



  Courbe plot(const float &x, const float &y, cstring format = "");

  Courbe plot_int(const Vecf &x, const Vecf &y, cstring format = "", cstring titre = "");

  Courbe plot_int(const Vecf &y, cstring format = "", cstring titre = "");

  template<typename T, typename ... Ts>
  Courbe plot(const VecT<T> &y, cstring format = "", cstring label = "", Ts &&... args)
  {
    if constexpr(est_complexe<T>())
    {
      if(abs2(imag(y)).somme() < abs2(real(y)).somme() * 1e-5)
      {
        return plot_int(real(y), format, sformat(FMT_RUNTIME(label), args...));
      }
      else
      {
        plot_int(real(y), format, sformat(FMT_RUNTIME(label), args...) + " (re)");
        return plot_int(imag(y), format, sformat(FMT_RUNTIME(label), args...) + " (im)");
      }
    }
    else
    {
      return plot_int(y, format, sformat(FMT_RUNTIME(label), args...));
    }
  }


  /** @brief Affiche une courbe
   *
   *  <h3>Affichage d'une courbe</h3>
   *  @param x      Vecteur des valeurs en abcisse.
   *  @param y      Vecteur des valeurs en ordonnée.
   *  @param format Chaine de caractère optionnelle pour modifier le format.
   *  @param label  Argument optionnel indiquant le nom de la courbe
   *  @param args   Argument optionnel pour générer le nom de la courbe (format librairie fmt).
   *
   *  Le format est une chaine de caractère qui permet de choisir la couleur, le type de trait,
   *  et le type de curseur :
   *    - Couleurs : "b" (black), "g" (green), "r" (red), "m" (magenta), "c" (cyan), etc.
   *    - Traits : "-" (lignes), "$|$" (barres), "" (aucun)
   *    - Curseurs : "s" (carré), "d" (diamant), "*" (étoile), "." (point), "o" (rond), "" (aucun)
   *
   *  @par Exemple, avec différentes couleurs et différents marqueurs :
   *  @snippet exemples/src/ex-vue.cc ex_figure
   *  @image html figure.png width=800px
   *
   *  @sa Figure::plot_psd(), Figure::plot_img()
   *
   */
  template<typename T, typename ... Ts>
  Courbe plot(const Vecf &x, const VecT<T> &y, cstring format = "", cstring label = "", Ts &&... args)
  {
    auto t = FMT_RUNTIME(label);
    if constexpr(est_complexe<T>())
    {
      if(abs2(imag(y)).somme() < abs2(real(y)).somme() * 1e-5)
      {
        return plot_int(x, real(y), format, sformat(t, args...));
      }
      else
      {
        plot_int(x, real(y), format, sformat(t, args...) + " (re)");
        return plot_int(x, imag(y), format, sformat(t, args...) + " (im)");
      }
    }
    else
    {
      return plot_int(x, y, format, sformat(t, args...));
    }
  }

  template<typename T, typename ... Ts>
  Courbe plot(const tuple<Vecf, VecT<T>> &xy, cstring format = "", cstring label = "", Ts &&... args)
  {
    soit [x, y] = xy;
    retourne plot(x, y, format, sformat(FMT_RUNTIME(label), args...));
  }

  Courbe plot_psd_int(const Veccf &y, float fe, cstring format, cstring label);

  /** @brief Affiche le spectre d'un signal
   *
   *  <h3>Affichage du spectre d'un signal</h3>
   *
   *  @param y      Signal à analyser (domaine temporel)
   *  @param fe     Fréquence d'échantillonnage
   *  @param format Spécification de format (optionnel)
   *  @param label  Nom de la courbe (optionnel)
   *  @param args   ...
   *
   *  @sa Figure::plot(), Figure::plot_img()
   */
  template<typename T, typename ... Ts>
  Courbe plot_psd(const VecT<T> &y, float fe = 1.0, cstring format = "", cstring label = "", Ts &&... args)
  {
    return plot_psd_int(y.as_complex(), fe, format, sformat(FMT_RUNTIME(label), args...));
  }



  /** @brief Dessine une surface 2d avec des niveaux de couleur.
   *
   *  <h3>Surface 2d</h3>
   *
   *  Cette méthode dessine une surface 2d grâce à des niveaux de couleurs (par défaut, carte "jet").
   *
   *  @param Z      Tableau 2d à afficher.
   *  @param format Spécification optionnelle de format.
   *  @param rdi    Rectangle où placer l'image.
   *
   *
   *  @image html plot-img.png width=600px
   *
   *
   *  @sa Figure::plot(), Figure::plot_psd()
   */
  Courbe plot_img(const Rectf &rdi, const Tabf &Z, cstring format = "jet");


  /** @brief Dessine une surface 2d avec des niveaux de couleur (dimensions par défaut). */
  Courbe plot_img(const Tabf &Z, cstring format = "jet");

  Courbe plot_minmax(const Vecf &x, const Vecf &y1, const Vecf &y2);

  Courbe plot_iq_int(const Veccf &z, cstring format, cstring label);

  template<typename ... Ts>
    Courbe plot_iq(const Veccf &z, cstring format = "", cstring label = "", Ts &&... args)
  {
    return plot_iq_int(z, format, sformat(FMT_RUNTIME(label), args...));
  }

  vector<Courbe> &courbes();

  void def_rdi_min(const Rectf &rdi);
  void def_rdi(const Rectf &rdi);

  Rectf get_rdi() const;

  Axes axes();

  /** @brief Définition du titre de la figure (ou de la sous-figure en cours).
   *
   * <h3>Définition du titre</h3>
   *
   * Cette fonction permet de changer le titre de la figure (ou de la sous-figure en cours).
   *
   * @param titre_global Texte affiché en haut de la figure.
   *
   * @sa Figure::titres()
   */
  void titre(cstring titre_global);

  /** @brief Définition des titres.
   *
   * <h3>Définition des titres</h3>
   *
   * @param titre_global Texte affiché en haut de la figure
   * @param axe_x Label de l'axe horizontal
   * @param axe_y Label de l'axe vertical
   *
   * @sa Figure::titre()
   */
  void titres(cstring titre_global, cstring axe_x = "", cstring axe_y = "");

  void def_pos_legende(cstring code);

  /** @brief  Affichage à l'écran dans une fenêtre à part.
   *
   * <h3>Affichage d'une figure</h3>
   *
   * @warning Cette fonction n'aura aucun effet si la dépendance optionnelle GTKMM n'est pas présente
   * au moment de la compilation de la librairie.
   *
   * @param nom_fenetre Titre de la fenêtre graphique.
   *
   * @par Exemple
   * @code
   * Figure f;
   * f.plot(x);
   * f.afficher();
   * @endcode
   *
   * @sa Figure::enregistrer()
   */
  //void afficher(cstring nom_fenetre = "");

  /** @brief Rendu vers une image bitmap */
  //Image rendre(uint32_t sx = 1200, uint32_t sy = 800, const Rectf &rdi, const Couleur &arp = Couleur::Blanc) const;

  Canva canva();
  Canva canva_pre();

  Canva canva_pixel(const Dim &allocation);


  void def_echelle(bouléen log_x, bouléen log_y);

  //void rendre(Canva canva) const;


  entier rendre00() const;
  Canva rendre0(Canva canva_) const;
  void rendre1(Canva canva_, Canva canva) const;


  /** @brief Enregistrement sous la forme d'un fichier image */
  //void enregistrer(cstring chemin_fichier, entier sx = -1, entier sy = -1) const;

  string lis_nom() const;
  void def_nom(cstring s);

  /** @brief Attends que l'utilisateur ferme toutes les fenêtres graphiques.
   */
  static void attente_ihm();

  Figure clone() const;

  sptr<const Rendable> rendable() const;


  _PIMPL_
  friend class FenetreFigure;
};


  struct Figures: ARendable
  {

    Figures(entier rows = -1, entier cols = -1);

    void clear();

    /** @brief Partitionne la figure en sous-plots.
      *
      *  <h3>Partition de la figure en sous-plots</h3>
      *  @param rows   Nombre de lignes
      *  @param cols   Nombre de colonnes
      *  @param sel    Sélection du subplot en cours (numéro entre 1 et rows*cols). Ou -1 pour automatique (nouveau subplot). */
     Figure subplot(entier rows, entier cols, entier sel = -1);

     /** @brief Partitionne la figure en sous-plots.
      *
      *  <h3>Partition de la figure en sous-plots</h3>
      *
      *  Cette surcharge est juste une manière plus compacte de spécifier un sous-plot.
      *
      *  @param p Entier entre 000 et 999, où les trois chiffres sont, respectivement,
      *  les nombre de lignes, de colonnes, et le numéro de subplot en cours (index commençant à 1).
      *  Si une valeur négative est passée en paramètre, un nouveau subplot est créé, et la position
      *  est déterminée de manière automatique.
      *  @par Exemple :
      *  @code
      *    Figure f;
      *    // 2 lignes, 1 colonne, première sous-figure
      *    f.subplot(211); // équivalent à f.subplot(2,1,1);
      *    f.plot(...);
      *    // 2 lignes, 1 colonne, deuxième sous-figure
      *    f.subplot(212); // équivalent à f.subplot(2,1,2);
      *    f.plot(...);
      *  @endcode
      *  @sa subplot() */
     Figure subplot(entier p = -1);

     vector<Figure> subplots(int n);

     Figure gcf();

     Figure gf(entier sel);


     /** @brief Enregistrement sous la forme d'un fichier image */
     //void enregistrer(cstring chemin_fichier, entier sx = -1, entier sy = -1) const;
     //void rendre(Canva canva) const;


     sptr<const Rendable> rendable() const;


     // Sur-charge
     void afficher(cstring titre = "", const Dim &dim = {-1, -1}) const;

     _PIMPL_
  };


  struct AFigures
  {
    Figures figures;
  };



  struct Stdo
  {
    void def_dossier_sortie(cstring chemin);

    void affiche(sptr<const Rendable> s, cstring titre = "", const Dim &dim = {-1,-1});

    void printf(cstring s);

    void flush();

    void fin();
  };

  extern Stdo stdo;

  /** @} */

  extern std::function<void ()> stdo_fin;
  extern std::function<void (sptr<const Rendable> fig, cstring titre)> stdo_ajoute_figure;

  // Active les réglages pour un mode impression par défaut
  extern void set_mode_impression(bool mode_impression = oui);

  extern bouléen get_mode_impression();

}

namespace tsd::ihm {
extern void init(bouléen create_gtk_thread = non);
}


#endif
