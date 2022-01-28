#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/vue/image.hpp"


namespace tsd::vue::unites {
extern std::tuple<int,int> calc_expo_nb_chiffres(double t, const std::string &unite);
extern std::string valeur_vers_chaine(double t, const std::string &unite, int expo, int nb_chiffres);
extern std::tuple<int, int> calc_expo_nb_chiffres_commun(const std::vector<double> &tics, const std::string &unite);
extern std::string valeur_vers_chaine(double t, const std::string &unite);
extern std::string valeur_vers_chaine(double t, const std::string &unite, int nb_chiffres);
}



namespace tsd::vue
{

/** @addtogroup vue
 *  @{
 */

extern Couleur eclaircir_couleur(const Couleur &src, float ratio = 0.5);



struct ConfigGrille
{
  int afficher = -1;
  int lg_tiret = -1, lg_trou = -1;
  Couleur couleur{0,0,0,0};
  int epaisseur = 1;

  // Par défaut : 50 pour majeur, 4 pour mineur
  int pixels_par_tic_min = -1;

  int sous_graduation = 5;

  void from_style(const std::string &style);
  std::string vers_style() const;

  std::strong_ordering operator <=>(const ConfigGrille &c) const = default;
};

struct ConfigAxe
{
  std::string id, label, unite;
  /** @brief For manual zoom selection. */
  float vmin = -1.0, vmax = 1.0;

  bool afficher                   = true;
  int ecart_pixel                 = -1;
  float ecart_valeur              = -1;

  bool echelle_logarithmique      = false;

  /** Si vrai, les tics affichés sont ceux précisés dans le vecteur "tics" */
  bool tic_manuels = false;
  std::vector<std::pair<float, std::string>> tics;


  float valeurs_echelle = 0.6;


  bool operator <=>(const ConfigAxe &c) const = default;

  operator std::string() const;
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
  std::string nom;
  Couleur couleur{0,0,0};
  int largeur_de_trait  = 1;
  std::string unite;
  bool visible          = true;
  bool remplissage      = false;
  std::string style;

  Marqueur marqueur = Marqueur::AUCUN;
  Trait trait = Trait::LIGNE;

  bool operator <=>(const ConfigSerie &c) const = default;
};



struct ConfigAxes
{
  void ajoute_serie(const std::string &nom);

  void def_rdi_visible(float xmin, float xmax, float ymin, float ymax);
  void def_rdi_visible_abs(float xmin, float xmax, float ymin, float ymax);

  /** @brief List of channels */
  std::vector<ConfigSerie> series;

  ConfigGrille grille_majeure, grille_mineure;
  ConfigAxe axe_vertical, axe_horizontal;

  std::string titre;

  int titre_dim = 10;

  /** @brief Background color */
  Couleur couleur_arriere_plan = Couleur{255,255,255};

  //bool affiche_cartouche_legende  = true;

  bool est_mat_sur_clair() const;

  struct Legende
  {
    bool afficher = true;

    struct Position
    {
      enum Type
      {
        NORD_EST = 0,
        SUD_EST,
        SUD_OUEST,
        NORD_OUEST
      } type = NORD_EST;

      operator std::string() const;
      Position(const std::string &id = "ne");

      std::strong_ordering operator <=>(const Position &c) const = default;
    } position;

    float dim = 0.5f;

    std::strong_ordering operator <=>(const Legende &c) const = default;

  } legende;



  //int legende_dim = 6;

  enum PositionAxeVertical
  {
    AXE_VERTICAL_GAUCHE = 0,
    AXE_VERTICAL_DROITE
  } position_axe_vertical = AXE_VERTICAL_DROITE;

  bool mode_isoview = false;

  float aspect_ratio = 1.0f;

  bool operator <=>(const ConfigAxes &c) const = default;

  operator std::string() const;
};

struct Axe
{
  ConfigAxe config;

  int dim; // nb pixels
  std::vector<double> tics_majeurs_pos;
  std::vector<double> tics_mineurs_pos;

  double ecart_tics_majeurs_valeur;
  int   ecart_tics_majeurs_pixels;

  bool tics_mineurs_presents;
  /** Center of X / Y axis, in lsb units */
  float   decalage_zero;

  /** This defines the display gain, for Y axis */
  float   pixels_par_lsb;

  bool vertical = false;

  void configure(const ConfigAxe &config, int dim);

  bool tic_ok(float v);
  bool tic_ok_souple(float v);

  float pixel_vers_valeur(float p) const;
  int valeur_vers_pixel(float v) const;
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

  // Soit :
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
  Canva clip(const Rectf &rdi, const Rectf &vue, bool log_x = false, bool log_y = false) const;
  Canva clip_alt(const Rectf &rdi, float xmin, float xmax, float ymin, float ymax, bool log_x, bool log_y) const;

  // Change l'échelle de la vue en cours
  Canva vue(const Rectf &vue) const;

  Pointf v2c(const Pointf &p) const;

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
  void set_epaisseur(int ep);
  void set_remplissage(bool remplir, const Couleur &coul = Couleur{255,255,255});
  void set_dotted(bool dotted);
  void set_align(Align hor, Align vert);
  void rectangle(float x0, float y0, float x1, float y1);
  void rectangle(const Rectf &r);
  void rectangle(const Pointf &p0, const Pointf &p1);
  void ligne(float x0, float y0, float x1, float y1);
  void ligne(const Pointf &p0, const Pointf &p1);
  void fleche(const Pointf &p0, const Pointf &p1);

  void texte(float x0, float y0,
      const std::string &texte,
      float dx_max = -1, float dy_max = -1);

  void texte(const Pointf &p,
      const std::string &texte,
      const Dimf &dim_max = {-1.0f, -1.0f});



  void cercle(float x0, float y0, float r);
  void ellipse(const Pointf &p, float a, float b, float theta = 0);
  void remplissage_vertical(float x0, float y0, float x1, float y1, float y2, float y3);
  void dessine_accu(const ArrayXcf &pts);
  void def_image_fond(Image img);
  void dessine_img(float x0, float y0, float y1, float y2, Image img);
  void marqueur(const Pointf &p, Marqueur m, float dim_pixels);
  void get_pixel_size(float &width, float &heigth);
  /** @brief Coordonnées pixel vers unités */
  Pointf coord_pixel(const Point &p);


  // TODO : rdi nécessaire ? n'est-il pas déjà précisé dans le clip interne ?
  Image rendre(const Dim &dim = {-1, -1}, const Rectf &rdi = {0,0,0,0}, const Couleur &arp = Couleur::Blanc);
  struct PointIntermediaire;

  sptr<PointIntermediaire> enregistre_pti();
  void restaure_pti(sptr<Canva::PointIntermediaire> pti);


  Rectf calc_rdi_englobante() const;

private:
  struct Impl;
  sptr<Impl> impl;
};


struct Rendable: std::enable_shared_from_this<Rendable>
{
  virtual void rendre(Canva canva) const = 0;

  /** @brief Enregistrement sous la forme d'un fichier image */
  void enregistrer(const std::string &chemin_fichier, const Dim &dim = {-1, -1}) const;

  virtual void afficher(const std::string &titre = "", const Dim &dim = {-1, -1}) const;

  Image genere_image(const Dim &dim = Dim{-1,-1}, const Couleur &arp = Couleur::Blanc) const;
  //Image rendre() const {return genere_image();}
};


struct ARendable
{
  virtual sptr<const Rendable> rendable() const = 0;

  void enregistrer(const std::string &chemin_fichier, const Dim &dim = {-1, -1}) const
  {
    rendable()->enregistrer(chemin_fichier, dim);
  }
  Image genere_image(const Dim &dim = Dim{-1,-1}, const Couleur &arp = Couleur::Blanc) const
  {
    return rendable()->genere_image(dim, arp);
  }
  void afficher(const std::string &titre = "", const Dim &dim = {-1, -1}) const
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

  void set_isoview(bool isoview);
  void fixe_aspect_ratio(float ysx);

  void active_grilles(bool actives);
  void supprime_decorations();

  /** @brief Définition de l'échelle d'affichage (lin ou log)
   *  @todo : n'a pas l'air de fonctionner ? */
  void def_echelle(const std::string &x, const std::string &y);


  Canva rendre(Canva canva_, float xmin, float xmax, float ymin, float ymax) const;
  Canva rendre(Canva canva, const Rectf &rdi) const;
  void post_rendu(Canva canva_) const;



  Dim get_dim_graphique(const Dim &dim_totale) const;
  void configure(const ConfigAxes &config);
  ConfigAxes &get_config();

  Axes clone();

private:
  struct Impl;
  sptr<Impl> impl;
};








  /** @brief %Figure avec API de type Matlab / Scilab.
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
    void def_couleurs(IArrayXf c, const std::string cmap = "");

    /** @brief Changement d'épaisseur du trait */
    void def_epaisseur(int ep);

    /** @brief Remplissage sous la courbe */
    void def_remplissage(bool actif, bool vers_vmin = false, float vmin = 0);

    /** @brief Définit, pour chaque point de la courbe, l'écart-type */
    void def_σ(IArrayXf σ);

    /** @brief Définit la dimension (en pixels) des marqueurs. */
    void def_dim_marqueur(int dim);
  //private:
    struct Impl;
    sptr<Impl> impl;
  };

  Figure(const std::string &nom = "");

  /** @brief Efface la figure */
  void clear();



  // void def_couleur(float r, float g, float b);



  void plot(const float &x, const float &y, const std::string &format = "");
  Courbe plot_int(const ArrayXf &x, const ArrayXf &y, const std::string &format = "", const std::string &titre = "");
  Courbe plot_int(const ArrayXf &y, const std::string &format = "", const std::string &titre = "");

  template<typename derived, typename ... Ts>
  Courbe plot(const Eigen::ArrayBase<derived> &y, const std::string &format = "", const std::string &titre = "", Ts &&... args)
  {
    if constexpr(est_complexe<typename derived::Scalar>())
    {
      if(y.imag().abs2().sum() < y.real().abs2().sum() * 1e-5)
      {
        return plot_int(y.real(), format, fmt::format(FMT_RUNTIME(titre), args...));
      }
      else
      {
        plot_int(y.real(), format, fmt::format(FMT_RUNTIME(titre), args...) + " (re)");
        return plot_int(y.imag(), format, fmt::format(FMT_RUNTIME(titre), args...) + " (im)");
      }
    }
    else
    {
      return plot_int(y, format, fmt::format(FMT_RUNTIME(titre), args...));
    }
  }


  /** @brief Affiche une courbe
   *
   *  <h3>Affichage d'une courbe</h3>
   *  @param x      Vecteur des valeurs en abcisse.
   *  @param y      Vecteur des valeurs en ordonnée.
   *  @param format Chaine de caractère optionnelle pour modifier le format.
   *  @param titre  Argument optionnel indiquant le nom de la courbe
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
  template<typename derived, typename ... Ts>
  Courbe plot(const ArrayXf &x, const Eigen::ArrayBase<derived> &y, const std::string &format = "", const std::string &titre = "", Ts &&... args)
  {
    auto t = FMT_RUNTIME(titre);
    if constexpr(est_complexe<typename derived::Scalar>())
    {
      if(y.imag().abs2().sum() < y.real().abs2().sum() * 1e-5)
      {
        return plot_int(x, y.real(), format, fmt::format(t, args...));
      }
      else
      {
        plot_int(x, y.real(), format, fmt::format(t, args...) + " (re)");
        return plot_int(x, y.imag(), format, fmt::format(t, args...) + " (im)");
      }
    }
    else
    {
      return plot_int(x, y, format, fmt::format(t, args...));
    }
  }

  Courbe plot_psd_int(IArrayXcf y, float fe, const std::string &format, const std::string &titre);

  /** @brief Affiche le spectre d'un signal
   *
   *  <h3>Affichage du spectre d'un signal</h3>
   *
   *  @param y      Signal à analyser (domaine temporel)
   *  @param fe     Fréquence d'échantillonnage
   *  @param format Spécification de format (optionnel)
   *  @param titre  Nom de la courbe (optionnel)
   *  @param args   ...
   *
   *  @sa Figure::plot(), Figure::plot_img()
   */
  template<typename Derived, typename ... Ts>
  Courbe plot_psd(const Eigen::ArrayBase<Derived> &y, float fe = 1.0, const std::string &format = "", const std::string &titre = "", Ts &&... args)
  {
    ArrayXcf y2 = y;
    return plot_psd_int(y2, fe, format, fmt::format(FMT_RUNTIME(titre), args...));
  }

  /** @brief Dessine une surface 2d avec des niveaux de couleur.
   *
   *  <h3>Surface 2d</h3>
   *
   *  Cette méthode dessine une surface 2d grâce à des niveaux de couleurs (par défaut, carte "jet").
   *
   *  @param Z Tableau 2d à afficher
   *  @param format Spécification optionnelle de format
   *
   *
   *  @sa Figure::plot(), Figure::plot_psd()
   */
  Courbe plot_img(IArrayXXf &Z, const std::string &format = "jet");

  Courbe plot_img(float xmin, float xmax, float ymin, float ymax, IArrayXXf &Z, const std::string &format = "jet");

  Courbe plot_minmax(const ArrayXf &x, const ArrayXf &y1, const ArrayXf &y2);

  Courbe plot_iq_int(const ArrayXcf &z, const std::string &format, const std::string &titre);

  template<typename ... Ts>
    Courbe plot_iq(const ArrayXcf &z, const std::string &format = "", const std::string &titre = "", Ts &&... args)
  {
    return plot_iq_int(z, format, fmt::format(FMT_RUNTIME(titre), args...));
  }



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
  void titre(const std::string &titre_global);

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
  void titres(const std::string &titre_global, const std::string &axe_x = "", const std::string &axe_y = "");

  void def_pos_legende(const std::string &code);

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
  //void afficher(const std::string &nom_fenetre = "");

  /** @brief Rendu vers une image bitmap */
  //Image rendre(uint32_t sx = 1200, uint32_t sy = 800, const Rectf &rdi, const Couleur &arp = Couleur::Blanc) const;

  Canva canva();

  Canva canva_pixel(const Dim &allocation);


  void def_echelle(bool log_x, bool log_y);

  //void rendre(Canva canva) const;


  int rendre00() const;
  Canva rendre0(Canva canva_) const;
  void rendre1(Canva canva_, Canva canva) const;


  /** @brief Enregistrement sous la forme d'un fichier image */
  //void enregistrer(const std::string &chemin_fichier, int sx = -1, int sy = -1) const;

  std::string lis_nom() const;
  void def_nom(const std::string &s);

  /** @brief Attends que l'utilisateur ferme toutes les fenêtres graphiques.
   */
  static void attente_ihm();

  Figure clone() const;

  sptr<const Rendable> rendable() const;

private:
  friend class FenetreFigure;
  struct Impl;
  sptr<Impl> impl;
};


  struct Figures: ARendable
  {

    Figures(int rows = -1, int cols = -1);

    void clear();

    /** @brief Partitionne la figure en sous-plots.
      *
      *  <h3>Partition de la figure en sous-plots</h3>
      *  @param rows   Nombre de lignes
      *  @param cols   Nombre de colonnes
      *  @param sel    Sélection du subplot en cours (numéro entre 1 et rows*cols). Ou -1 pour automatique (nouveau subplot). */
     Figure subplot(int rows, int cols, int sel = -1);

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
     Figure subplot(int p = -1);

     Figure gcf();


     /** @brief Enregistrement sous la forme d'un fichier image */
     //void enregistrer(const std::string &chemin_fichier, int sx = -1, int sy = -1) const;
     //void rendre(Canva canva) const;


     sptr<const Rendable> rendable() const;


     // Sur-charge
     void afficher(const std::string &titre = "", const Dim &dim = {-1, -1}) const;

  private:
    struct Impl;
    sptr<Impl> impl;
  };



  struct EvtPlotConfig
  {
    struct Evt
    {
      std::string nom;
      tsd::vue::Couleur coul = Couleur::Blanc;
    };
    std::vector<Evt> evts;
    int idle_index = -1;
  };



  struct PlotEvt: Rendable
  {
    ArrayXi evt;
    EvtPlotConfig cfg;

    void maj(ArrayXi &evt, const EvtPlotConfig &cfg);
    void rendre(Canva canva) const;
  };


  struct Stdo
  {
    void def_dossier_sortie(const std::string &chemin);

    void affiche(sptr<const Rendable> s, const std::string &titre = "", const Dim &dim = {-1,-1});

    void printf(const std::string &s);

    void flush();

    void fin();
  };

  extern Stdo stdo;

  /** @} */

  extern std::function<void ()> stdo_fin;
  extern std::function<void (sptr<const Rendable> fig, const std::string &titre)> stdo_ajoute_figure;

  // Active les réglages pour un mode impression par défaut
  extern void set_mode_impression();

}


