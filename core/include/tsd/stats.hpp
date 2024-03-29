#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"

namespace tsd::stats {


  template<typename T>
    T variance(const Vecteur<T> &vec)
  {
    return carré(vec).moyenne() - carré(vec.moyenne());
  }

  template<typename T>
    T stdev(const Vecteur<T> &vec)
  {
    return sqrt(variance(vec));
  }


  /** @addtogroup stats
   *  @{
   */

  /** @brief Résolution de l'équation @f$Ra = -r@f$ par l'algorithme de Levinson - Durbin
   *
   * Résolution du système de Toeplitz
   * @f[
\left(\begin{array}{cccc}
r_0 &r_1 &\cdots &r_{n-1}\\
r_1 &r_0 &\cdots& r_{n-2}\\
\vdots &\ddots &\ddots &\vdots\\
r_{n-1} &r_{n-2}& \cdots &r_0\\
\end{array}\right)
\cdot
\left(\begin{array}{c}
a_1\\
\vdots\\
a_n\\
\end{array}\right)
=
-\left(\begin{array}{c}
r_1\\
r_2\\
\vdots\\
r_n\\
\end{array}\right)
   * @f]
   * par l'algorithme de Levinson-Durbin.
   *
   * @param r Vecteur d'auto-corrélation du signal
   * @returns Vecteur @f$a@f$ (@f$n+1@f$ coefficients), avec @f$a_0=1@f$.
   *
   * Cette méthode est typiquement utilisable pour un problème de prédication linéaire, où l'on cherche à modéliser un
   * processus @f$x@f$ suivant un filtre AR :
   * @f[
   *  x_n = \sum_{k=1}^{n} a_k x_{n-k} + e_n
   * @f] */
  Vecf levinson_reel(const Vecf &r);


  /** @brief Récursion de Levinson - Durbin (cas général)
   *
   *  Résolution du système de Toeplitz général
   * @f[
      \left(\begin{array}{cccc}
      r_0 &r_{-1} &\cdots &r_{-(n-1)}\\
      r_1 &r_0 &\cdots& r_{-(n-2)}\\
      \vdots &\ddots &\ddots &\vdots\\
      r_{n-1} &r_{n-2}& \cdots &r_0\\
      \end{array}\right)
      \cdot
      \left(\begin{array}{c}
      x_1\\
      \vdots\\
      x_n\\
      \end{array}\right)
      =
      \left(\begin{array}{c}
      y_1\\
      y_2\\
      \vdots\\
      y_n\\
      \end{array}\right)
   * @f]
   * par l'algorithme de Levinson-Durbin (https://en.wikipedia.org/wiki/Levinson_recursion).
   *
   *
   *  @param l1   Première ligne de la matrice de Toeplitz
   *  @param c1   Première colonne de la matrice
   *  @param y    Second membre
   *  @return     Solution @f$x@f$
   *
   *  */
  extern Vecf levinson(const Vecf &l1, const Vecf &c1, const Vecf &y);

  /** @brief Calcul d'une matrice d'auto-corrélation à partir d'un vecteur de corrélations
   *
   * Cette fonction construit une matrice d'auto-corrélation (matrice de Toeplitz pour un signal stationnaire) à partir du vecteur
   * d'auto-corrélation :
   *
   * @f[
      R = \left(\begin{array}{cccc}
      r_0 &r_1^\star &\cdots &r_{n-1}^\star\\
      r_1 &r_0 &\cdots& r_{n-2}^\star\\
      \vdots &\ddots &\ddots &\vdots\\
      r_{n-1} &r_{n-2}& \cdots &r_0\\
      \end{array}\right)
   * @f]
   *
   * @param r Auto-corrélation du signal
   * @return R
   */
  template<typename T>
  auto r_vers_R(const Vecteur<T> &r)
  {
    soit n = r.rows();
    TabT<T,2> R(n, n);
    Pour(auto i = 0; i < n; i++)
    {
      Pour(auto j = 0; j < n; j++)
      {
        R(i,j) = r(abs(i-j));
        Si constexpr(est_complexe<T>())
          Si(j > i)
            R(i,j) = conj(R(i,j));
      }
    }
    return R;
  }

  /** @brief Calcul de la matrice de covariance pour un signal supposé <b>stationnaire</b> et <b>centré</b>.
   *
   *  Calcul de la matrice de covariance pour un signal supposé <b>stationnaire</b> et <b>centré</b>.
   *  @param x Signal à analyser
   *  @param m Dimension de la matrice (nombre de délais examinés)
   *  @returns Matrice de covariance :
   *  @f[
   *  R_{ij} = \mathbb{E}\left[x_{.+i} \cdot x_{.+j}^\star\right]
   *  @f]
   **/
  template<typename T>
  auto covmtx(const Vecteur<T> &x, entier m)
  {
    // Si le signal est stationnaire :
    // R_ij = cov(x_(k+i) , x_(k+j)) = cov(x_0, x_(j-i))
    // Donc ligne 0 = cov(x0,x0), cov(x0,x1), etc.
    //      ligne 1 = cov(x1,x0), cov(x1,x1), etc.
    //              = cov(x0,x1)*, cov(x0,x0) etc.
    // Donc matrice auto-ajointe ou symétrique ?
    soit [lags, cr] = tsd::fourier::xcorr(x, x, m);
    return r_vers_R(cr.tail(m));
  }


  /** @brief Analyse LPC (prédiction linéaire)
   *
   * Cette fonction calcul les coefficients @f$a_k@f$ du modèle AR :
   * @f[
   *  x_n = \sum_{k=1}^{n} a_k x_{n-k} + e_n
   * @f]
   *
   * @param x Signal à analyser
   * @param p Nombre de coefficients du modèle
   * @returns Les coefficients du filtre @f$a(z)@f$ et le vecteur d'erreur.
   */
  extern tuple<Vecf, Vecf> lpc(const Vecf &x, entier p);

  /** @brief %Filtre de Wiener (RIF)
   *
   *  Modèle :
   *  @f[
   *  y = g \star x + b
   *  @f]
   *  et recherche du filtre RIF optimal @f$h@f$ (au sens des moindres carrés) tel que
   *  @f[
   *  h \star y \sim x
   *  @f]
   *
   *  @param Rxy Matrice de corrélation croisée entre les signaux @f$x@f$ (entrées) et @f$y@f$ (observations)
   *  @param rx  Vecteur d'auto-corrélation du signal @f$x@f$
   *  @param p   Ordre du modèle (nombre de coefficients du filtre @f$h@f$)
   *  @return Coefficients du filtre d'égalisation optimal (filtre @f$h@f$)
   *
   */
  extern Vecf wiener_rif(const Tabf &Rxy, const Vecf &rx, entier p);


  /** @brief Paramètrage analyse sous-espace */
  struct SubSpaceSpectrumConfig
  {
    /** @brief Méthode de calcul */
    enum Method
    {
      /** @brief MUltiple SIgnal Classification */
      MUSIC,
    } methode = MUSIC;

    /** @brief Si activé, trace des figures */
    bouléen debug_actif = non;

    /** @brief Nombre d'exponentielles complexes / de signaux sources (ou -1 pour détermination automatique) */
    entier Ns = 0;

    /** @brief Nombre de points d'analyse */
    entier Nf = 1000;

    /** @brief Si vrai, la recherche est effectuée parmi les exponentielles complexes,
     *         par une méthode rapide type FFT, et la fonction de balayage n'est pas utilisée. */
    bouléen est_sigexp = oui;

    /** @brief Fonction de balayage
     *
     *  Fonction de balayage (par défaut simple exponentielle complexe, de fréquence dans l'intervalle [-0.5,0.5[,
     *  c'est à dire adaptée à la recherche de signaux périodiques)
     *  - Paramètre i : index
     *  - Paramètre n : nombre de points
     *  - Paramètre m : dimension du vecteur à retourner (nombre de récepteurs en DOA) */
    fonction<tuple<Vecf, Veccf> (entier i, entier n, entier m)> balayage = [](entier i, entier n, entier m) -> tuple<Vecf, Veccf>
    {
      // Calcul du "vecteur de steering"
      Vecf f(1);
      f(0) = (1.0f*i)/n-0.5f;
      return {f, sigexp(f(0), m)};
    };
  };

  /** @brief Résultat d'une analyse de sous-espace */
  struct SubSpaceSpectrum
  {
    /** @brief Valeurs des variables pour chaque point de spectre.
     *
     * (index ligne = valeurs de spectre, index colonne = variable)
     * Par exemple, pour une simple estimation de fréquence, c'est un simple vecteur colonne avec les fréquences normalisée. */
    Tabf var;

    /** @brief Valeurs du spectre */
    Vecf spectrum;

    /** @brief Nombre de sources détectées */
    entier Ns = 0;
  };


  /** @brief Calcul d'une réponse générale par la méthode des sous-espaces
   *
   *  Calcul d'une réponse suivant l'algorithme MUSIC :
   *  @f[
   *  S_i = \frac{1}{\sum_{k=p+1}^N \left|v_k^H \cdot e_i \right|}
   *  @f]
   *
   *  les @f$v_k@f$ étant les vecteurs propres associés au bruit
   *  (en supposant que @f$p@f$ premiers vecteurs propres soient associés au signal),
   *  et les @f$e_i@f$ étant les signaux tests (typiquement, des exponentielles complexes, pour la recherche
   *  de signaux périodiques).
   *
   *  @param R  Matrice d'auto-corrélation
   *  @param config Configuration (voir @ref SubSpaceSpectrumConfig)
   *  @returns La réponse générale (voir structure @ref SubSpaceSpectrum)
   *
   *  @par Bibliographie
   *  <i>Statistical signal processing and modelling</i>, M.H. Hayes, 1996
   */
  extern SubSpaceSpectrum subspace_spectrum(const Tabcf &R, const SubSpaceSpectrumConfig &config);



  /** @} */


}
