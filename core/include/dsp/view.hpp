#pragma once

#include "dsp/dsp.hpp"
#include "tsd/vue/image.hpp"
#include "tsd/vue.hpp"
#include "dsp/time.hpp"

namespace dsp::view
{

/** @addtogroup vue
 *  @{
 */

  using Dim    = tsd::vue::Dim;
  using Rectf  = tsd::vue::Rectf;
  using Rect   = tsd::vue::Rect;
  using Pointf = tsd::vue::Pointf;
  using Point  = tsd::vue::Point;
  using Color  = tsd::vue::Couleur;

  /** @brief %Figure with Matlab / Scilab type API.
   *
   *  @par Example:
   *  @snippet exemples/src/ex-vue.cc ex_figure
   *  @image html figure.png
   */
  struct Figure: tsd::vue::ARendable
  {
    tsd::vue::Figure f;

    Figure(tsd::vue::Figure f_): f(f_){}

    /** @brief %Courbe (telle que renvoyée par la méthode plot()).
     *
     * Cette classe permet de modifier a posteriori les propriétés d'un courbe (épaisseur de trait, couleur, etc.)
     *
     */
    struct Curve
    {
      tsd::vue::Figure::Courbe courbe;

      Curve(tsd::vue::Figure::Courbe c){courbe = c;}

      /** @brief Changement de couleur */
      void set_color(const tsd::vue::Couleur &color)
      {
        courbe.def_couleur(color);
      }

      /** @brief Définit une couleur différente pour chaque point de la courbe */
      void set_colors(const Vecf &c, const std::string cmap = "")
      {
        courbe.def_couleurs(c.fr, cmap);
      }

      /** @brief Changement d'épaisseur du trait */
      void set_width(int width)
      {
        courbe.def_epaisseur(width);
      }

      /** @brief Remplissage sous la courbe */
      void set_fill(bool active)
      {
        courbe.def_remplissage(active);
      }

      /** @brief Définit, pour chaque point de la courbe, l'écart-type */
      void set_σ(const Vecf &σ)
      {
        courbe.def_σ(σ.fr);
      }

      /** @brief Définit la dimension (en pixels) des marqueurs. */
      void set_marker_width(int width)
      {
        courbe.def_dim_marqueur(width);
      }
    };

    Figure(const std::string &name = ""):f(name){}

    /** @brief Efface la figure */
    void clear(){f.clear();}


    Curve plot(const float &x, const float &y, const std::string &format = ""){return f.plot(x,y,format);}

    template<typename T, typename ... Ts>
    Curve plot(const Vector<T> &y, const std::string &format = "", const std::string &label = "", Ts &&... args)
    {
      return f.plot(y.fr, format, label, args...);
    }



    /** @brief Affiche une courbe
     *
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
     *  @image html figure.png
     *
     *  @sa Figure::plot_psd(), Figure::plot_img()
     *
     */
    template<typename ... Ts>
    Curve plot(const Vecf &x, const Vecf &y, const std::string &format = "", const std::string &label = "", Ts &&... args)
    {
      return f.plot(x.fr, y.fr, format, fmt::format(FMT_RUNTIME(label), args...));
    }

    /** @brief Affiche le spectre d'un signal
     *
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
    Curve plot_psd(const Vector<T> &y, float fe = 1.0, const std::string &format = "", const std::string &label = "", Ts &&... args)
    {
      return f.plot_psd(y.fr, fe, format, label, args...);
    }

    /** @brief Dessine une surface 2d avec des niveaux de couleur.
     *
     *
     *  Cette méthode dessine une surface 2d grâce à des niveaux de couleurs (par défaut, carte "jet").
     *
     *  @param Z Tableau 2d à afficher
     *  @param format Spécification optionnelle de format
     *
     *
     *  @sa Figure::plot(), Figure::plot_psd()
     */
    Curve plot_img(const Tabf &Z, const std::string &format = "jet")
    {
      return f.plot_img(Z.fr, format);
    }

    Curve plot_img(const Rectf &rdi, const Tabf &Z, const std::string &format = "jet")
    {
      return f.plot_img(rdi, Z, format);
    }

    Curve plot_minmax(const Vecf &x, const Vecf &y1, const Vecf &y2)
    {
      return f.plot_minmax(x.fr, y1.fr, y2.fr);
    }

    template<typename ... Ts>
    Curve plot_iq(const Veccf &z, const std::string &format = "", const std::string &label = "", Ts &&... args)
    {
      return f.plot_iq(z.fr, format, label, args...);
    }



    void set_rdi_min(const Rectf &rdi)
    {
      f.def_rdi_min(rdi);
    }
    void set_rdi(const Rectf &rdi)
    {
      f.def_rdi(rdi);
    }

    tsd::vue::Axes axes()
    {
      return f.axes();
    }

    /** @brief Définition du titre de la figure (ou de la sous-figure en cours).
     *
     *
     * Cette fonction permet de changer le titre de la figure (ou de la sous-figure en cours).
     *
     * @param titre_global Texte affiché en haut de la figure.
     *
     * @sa Figure::titres()
     */
    void title(const std::string &titre_global)
    {
      f.titre(titre_global);
    }

    /** @brief Définition des titres.
     *
     * @param titre_global Texte affiché en haut de la figure
     * @param axe_x Label de l'axe horizontal
     * @param axe_y Label de l'axe vertical
     *
     * @sa Figure::titre()
     */
    void titles(const std::string &titre_global, const std::string &axe_x = "", const std::string &axe_y = "")
    {
      f.titres(titre_global, axe_x, axe_y);
    }

    void set_legend_position(const std::string &code)
    {
      f.def_pos_legende(code);
    }

    void enable_time_axis(const dsp::time::DateTime &t0)
    {
      f.active_axe_temporel(t0.dt());
    }


    tsd::vue::Canva canva()
    {
      return f.canva();
    }

    tsd::vue::Canva canva_pre()
    {
      return f.canva_pre();
    }

    tsd::vue::Canva canva_pixel(const Dim &allocation)
    {
      return f.canva_pixel(allocation);
    }

    /** @brief  Affichage à l'écran dans une fenêtre à part.
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
    void show(const std::string &nom_fenetre = "")
    {
      f.afficher(nom_fenetre);
    }

    /** @brief Rendu vers une image bitmap */
    /*tsd::vue::Image render(uint32_t sx = 1200, uint32_t sy = 800) const
    {
      return f.rendre(sx, sy);
    }*/

    /** @brief Enregistrement sous la forme d'un fichier image */
    void save(const std::string &chemin_fichier, const Dim &dim = {-1, -1}) const
    {
      f.enregistrer(chemin_fichier, dim);
    }

    std::string get_name() const
    {
      return f.lis_nom();
    }
    void set_name(const std::string &s)
    {
      f.def_nom(s);
    }

    Figure clone() const
    {
      return f.clone();
    }

    sptr<const tsd::vue::Rendable> rendable() const
    {
      return f.rendable();
    }
  };


  struct Figures
  {
    tsd::vue::Figures f;

    Figures(int rows = -1, int cols = -1): f(rows, cols){}
    Figures(const tsd::vue::Figures &f)
    {
      this->f = f;
    }

    void clear()
    {
      f.clear();
    }

    /** @brief Partitionne la figure en sous-plots.
      *
      *  @param rows   Nombre de lignes
      *  @param cols   Nombre de colonnes
      *  @param sel    Sélection du subplot en cours (numéro entre 1 et rows*cols). Ou -1 pour automatique (nouveau subplot). */
     Figure subplot(int rows, int cols, int sel = -1)
     {
       return f.subplot(rows, cols, sel);
     }

     /** @brief Partitionne la figure en sous-plots.
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
     Figure subplot(int p = -1)
     {
       return f.subplot(p);
     }

     vector<Figure> subplots(int n)
     {
       vector<Figure> res;
       let l = f.subplots(n);
       for(auto &s: l)
         res.push_back(s);
       return res;
     }

     Figure gcf()
     {
       return f.gcf();
     }

     Figure gf(int sel)
     {
       return f.gf(sel);
     }


     sptr<const tsd::vue::Rendable> rendable() const
     {
       return f.rendable();
     }


     // Sur-charge
     void show(const std::string &title = "", const Dim &dim = {-1, -1}) const
     {
       f.afficher(title, dim);
     }

     void save(const std::string &chemin_fichier, const Dim &dim = {-1, -1}) const
     {
       f.enregistrer(chemin_fichier, dim);
     }
  };




  /** @} */

}


