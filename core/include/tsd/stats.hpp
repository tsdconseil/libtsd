#pragma once

/** (C) 2022 J. Arzi / GPL V3 - voir fichier LICENSE. */

#include "tsd/tsd.hpp"
#include "tsd/fourier.hpp"

namespace tsd::stats {


  /** @addtogroup stats
   *  @{
   */

  /** @brief Résolution de l'équation @f$Ra = -r@f$
   *
   * <h3>Algorithme de Levinson - Durbin</h3>
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
  ArrayXf levinson_reel(const ArrayXf &r);


  /** @brief Récursion de Levinson - Durbin (cas général)
   *
   *  <h3>Récursion de Levinson - Durbin</h3>
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
  extern ArrayXf levinson(const ArrayXf &l1, const ArrayXf &c1, const ArrayXf &y);

  /** @brief Calcul d'une matrice d'auto-corrélation à partir d'un vecteur de corrélations
   *
   * <h3>Calcul de la matrice d'auto-corrélation</h3>
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
  template<typename derived>
  auto r_vers_R(const Eigen::ArrayBase<derived> &r)
  {
    using T = typename derived::Scalar;
    using rtype = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    int n = r.rows();
    rtype R(n, n);
    for(auto i = 0; i < n; i++)
    {
      for(auto j = 0; j < n; j++)
      {
        R(i,j) = r(std::abs(i-j));
        if constexpr(est_complexe<T>())
          if(j > i)
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
  auto covmtx(const Vecteur<T> &x, int m)
  {
    // Si le signal est stationnaire :
    // R_ij = cov(x_(k+i) , x_(k+j)) = cov(x_0, x_(j-i))
    // Donc ligne 0 = cov(x0,x0), cov(x0,x1), etc.
    //      ligne 1 = cov(x1,x0), cov(x1,x1), etc.
    //              = cov(x0,x1)*, cov(x0,x0) etc.
    // Donc matrice auto-ajointe ou symétrique ?
    auto [lags, cr] = tsd::fourier::xcorr(x, x, m);
    return r_vers_R(cr.tail(m));
  }


  /** @brief Analyse LPC (prédiction linéaire)
   *
   * <h3>Analyse LPC (prédiction linéaire)</h3>
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
  extern std::tuple<ArrayXf, ArrayXf> lpc(const ArrayXf &x, int p);

  /** @brief %Filtre de Wiener (RIF)
   *
   *  <h3>%Filtre de Wiener (RIF)</h3>
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
  extern ArrayXf wiener_rif(const MatrixXf &Rxy, const ArrayXf &rx, int p);


  /** @brief Paramètrage analyse sous-espace */
  struct SubSpaceSpectrumConfig
  {
    /** @brief Méthode de calcul */
    enum
    {
      /** @brief MUltiple SIgnal Classification */
      MUSIC,
    } methode = MUSIC;

    /** @brief Si activé, trace des figures */
    bool debug_actif = false;

    /** @brief Nombre d'exponentielles complexes / de signaux sources (ou -1 pour détermination automatique) */
    int Ns = 0;

    /** @brief Nombre de points d'analyse */
    int Nf = 1000;

    /** @brief Si vrai, la recherche est effectuée parmi les exponentielles complexes,
     *         par une méthode rapide type FFT, et la fonction de balayage n'est pas utilisée. */
    bool est_sigexp = true;

    /** @brief Fonction de balayage
     *
     *  Fonction de balayage (par défaut simple exponentielle complexe, de fréquence dans l'intervalle [-0.5,0.5[,
     *  c'est à dire adaptée à la recherche de signaux périodiques)
     *  - Paramètre i : index
     *  - Paramètre n : nombre de points
     *  - Paramètre m : dimension du vecteur à retourner (nombre de récepteurs en DOA) */
    std::function<std::tuple<ArrayXf, ArrayXcf> (int i, int n, int m)> balayage = [](int i, int n, int m) -> std::tuple<ArrayXf, ArrayXcf>
    {
      // Calcul du "vecteur de steering"
      ArrayXf f(1);
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
    ArrayXXf var;

    /** @brief Valeurs du spectre */
    ArrayXf spectrum;

    /** @brief Nombre de sources détectées */
    int Ns = 0;
  };


  /** @brief Calcul d'une réponse générale par la méthode des sous-espaces
   *
   *  <h3>Méthode des sous-espaces</h3>
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
   *  @param config Configuration
   *  @returns La réponse générale
   *
   *  @par Bibliographie
   *  <i>Statistical signal processing and modelling</i>, M.H. Hayes, 1996
   */
  extern SubSpaceSpectrum subspace_spectrum(const MatrixXcf &R, const SubSpaceSpectrumConfig &config);



  /** @} */


}
