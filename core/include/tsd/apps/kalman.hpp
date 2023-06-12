#pragma once

#include "tsd/tsd.hpp"
#include "tsd/vue.hpp"
#include "tsd/geometrie.hpp"

namespace tsd::kalman {

  /** @addtogroup kalman
   *  @{
   */

/** @brief Modèle abtrait pour une représentation d'état.
 *
 * Modèle général, non linéaire :
 * @f[
 * \begin{align}x_n &= f(x_{n-1}) + v_n\\ y_n &= g(x_n) + w_n\end{align}
 * @f]
 *
 * Where @f$u_n@f$ is the system input,
 * @f$x_n@f$ the system state, @f$y_n@f$ the system output
 * (observation), and @f$v_n@f$ and @f$w_n@f$ are
 * white gaussian noise processes of covariance matrix @f$Q@f$ and @f$R@f$.
 *
 *  */
struct SSM
{
  /** @brief Nombre d'états */
  entier ns = 1;

  /** @brief Nombre d'observations */
  entier no = 1;

  /** @brief Nombre d'entrées */
  entier ni = 0;

  /** @brief State evolution covariance matrix (default is 0) */
  Tabf Q;

  /** @brief Observation covariance matrix (default is 0) */
  Tabf R;

  /** @brief Square roots of R & Q */
  Tabf LR, LQ;

  string nom;

  vector<string> noms_etats, noms_obs;

  string lis_nom_etat(entier i) const;
  string lis_nom_obs(entier i) const;

  /** @brief transition function (ns x 1 vector to ns x 1 vector) */
  virtual Vecf f(const Vecf &x) const = 0;

  /** @brief observation function (ns x 1 vector to no x 1 vector) */
  virtual Vecf g(const Vecf &x) const = 0;

  /** @brief Function to compute the jacobian matrix of the transition function f (ns x 1 vector to ns x ns matrix) */
  virtual Tabf Jf(const Vecf &x) const;

  /** @brief Function to compute the jacobian matric of the observation function g (ns x 1 vector to ns x no matrix) */
  virtual Tabf Jg(const Vecf &x) const;


  /** @brief  Simulation of a state / space model
   *
   * <h3>Simulation of a state / space model</h3>
   *
   * Simulate a state / space model state evolution sequence, injecting noise
   * in the state evolution rule  and in the
   * observations vectors (using the Q and R covariance matrix,
   * as specified in the input model).
   *
   * @param n  Number of iterations to compute
   * @param x0 Initial state vector (column vector), default is 0.
   * @param vin Vecteur d'entrée (contrôle) optionnel
   * @return {x,y} : state & observation vectors
   */
  tuple<Tabf, Tabf> steps(entier n, const Vecf &x0 = Vecf(), const Tabf &vin = Tabf());

  // Checks the validity of the user supplied Jacobians for a EKF model
  //
  // Calling Sequence
  //  ekf_check(ekf)
  //
  // Description
  //  This function generates an error if the supplied state evolution and observation
  //  Jacobians functions results do not match with numerically computed Jacobians.
  //
  /** @brief Vérifie la validité des Jacobienne fournies.
   *
   *  <h3>Vérification de la validité des Jacobienne fournies</h3>
   *
   *  Cette fonction vérifie les Jacobiennes fournies pour la fonction d'évaluation d'état (Jf) et la
   *  fonction d'observation (Jg) sont cohérentes avec des Jacobiennes calculées numériquement à partir
   *  de f et g.
   */
  void verifications() const;

  entier k = 0;

  /** @brief Etat en cours */
  Vecf etat;
};




/** @brief Modèle linéaire
 *
 *  @f[
 *  \begin{align}x_n &= Ax_{n-1} + Bu_n + v_n\\ y_n &= Cx_n + Du_n + w_n\end{align}
 *  @f]
 */
struct SSMLineaire: SSM
{
  // A,B,C,D: System matrix for a linear system definition
  /** @brief Matrice de process (dépendance en fonction de l'état précédent). */
  Tabf A;

  /** @brief Matrice de process (dépendance en fonction du vecteur d'entrée). */
  Tabf B;

  /** @brief Matrice d'observation (dépendance en fonction de l'état) */
  Tabf C;

  /** @brief Matrice d'observation (dépendance en fonction de l'entrée) */
  Tabf D;

  /** @brief transition function (ns x 1 vector to ns x 1 vector) */
  Vecf f(const Vecf &x) const;

  /** @brief observation function (ns x 1 vector to no x 1 vector) */
  Vecf g(const Vecf &x) const;

  /** @brief Function to compute the jacobian matrix of the transition function f (ns x 1 vector to ns x ns matrix) */
  Tabf Jf(const Vecf &x) const;

  /** @brief Function to compute the jacobian matric of the observation function g (ns x 1 vector to ns x no matrix) */
  Tabf Jg(const Vecf &x) const;
};


extern sptr<SSMLineaire> ssm_lineaire(
    const Tabf &A,
    const Tabf &C,
    const Tabf &B = Tabf(),
    const Tabf &D = Tabf());



/** @brief Dérivation numérique.
 *
 * <h3>Dérivation numérique</h3>
 *
 * Calcul numérique des dérivées partielles (matrice Jacobienne) en un point donné :
 * @f[
 * J_{ij} = \frac{df_j}{dx_i} \sim \frac{f_j(x_i + \epsilon) - f(x_i)}{\epsilon}
 * @f]
 *
 * @warning Dans la mesure du possible, préférer une dérivée analytique, qui sera en général beaucoup
 * plus précise.
 *
 */
extern Tabf jacobienne_num(fonction<Vecf(const Vecf &)> f, const Vecf &x);




/** @brief %Filtre abstrait pour retrouver la séquence d'états cachés d'un modèle SSM.
 *
 *
 * @sa filtre_kalman(), filtre_ekf()
 */
struct FiltreSSM
{
  virtual ~FiltreSSM(){}

  /** @brief Prédiction de l'état caché à partir des vecteurs d'obvervation et d'entrée */
  virtual Vecf step(const Vecf &y, const Vecf &u = Vecf()) = 0;

  // Proceed to one or more steps of the standard Kalman filter
  /** @brief Appel direct sur une matrice d'observations */
  Tabf steps(const Tabf &y, const Tabf &u = Tabf());

  sptr<SSM> ssm;
  Vecf etat;
};


/** @brief Construction d'un filtre de Kalman.
 *
 * <h3>%Filtre de Kalman (modèle linéaire)</h3>
 *
 * Cette fonction renvoie un filtre de Kalman standard, fonctionnant sur un modèle linéaire.
 *
 * @param ssm Modèle à inverser
 * @param x0  Etat initial
 * @param p0  Matrice de covariance (incercitude sur l'état initial)
 *
 * @sa filtre_ekf()
 */
extern sptr<FiltreSSM> filtre_kalman(sptr<SSMLineaire> ssm,
                                     const Vecf &x0 = Vecf(),
                                     const Tabf &p0 = Tabf());

/** @brief Construction d'un filtre de Kalman étendu (EKF).
 *
 *  <h3>%Filtre de Kalman étendu (EKF / Extended Kalman Filter)</h3>
 *
 *
 * @param ssm Modèle à inverser
 *
 * @sa filtre_kalman()
 *
 */
extern sptr<FiltreSSM> filtre_ekf(sptr<SSM> ssm);


// Random walk model
//
// Calling Sequence
// ssm = model_rw(ndim, Q, R);
// ssm = model_rw(); // Default is ndims = 1, and Q = R = 1
//
// Parameters
// ndim: Number of dimensions (1 means scalar model)
// Q: Covariance of process noise
// R: Covariance of observation noise
// ssm: Output state / space model
//
// Description
// This function create a simple scalar or vectorial random walk model, e.g. at each time interval, the state is incremented by a random value, and observed through a noisy measurement. That is:
//  <latex>$$x_n = x_{n-1} + u_n,\quad y_n = x_n+v_n$$</latex>,
// with <latex>$$u_n: N(0,R), v_n: N(0,Q)$$</latex>.
//
// Random walk model example (2D)

/** @brief Marche aléatoire.
 *
 *  <h3>Marche aléatoire</h3>
 *
 *  @param ndim Nombre de dimensions (1 = modèle scalaire)
 *  @param Q Covariance du bruit de process
 *  @param R Covariance du bruit d'observation
 *  @return Modèle en représentation d'état
 *
 *  This function create a simple scalar or vectorial random walk model, e.g. at each time interval, the state is incremented by a random value, and observed through a noisy measurement. That is:
 *  @f[
 *  x_n = x_{n-1} + u_n,\quad y_n = x_n+v_n@f]
 *  avec @f$u_n: \mathcal{N}(0,R), v_n: \mathcal{N}(0,Q)@f$.
 *
 *  @par Exemple en dimension 2 :
 *  @snippet exemples/src/ex-kalman.cc ex_modele_ma
 *  @image html ex-modele-ma.png width=800px
 *
 *  @sa ssm_simu(), model_kitagawa(), model_imu(), model_constant()
 *
 */
extern sptr<SSMLineaire> modele_marche_aleatoire(entier ndim = 1, const Tabf &Q = Tabf(), const Tabf &R = Tabf());


// Most simple state / space system: constant scalar value
//
//  R: Observation noise variance
// ssm: Output state / space model

// <refsection><title>Example</title></refsection>
// <programlisting>
//ssm = model_constant();
//[x,y] = ssm_simu(ssm, 100, 5);
//clf();
//plot(x,'b-'); plot(y,'m.');
// </programlisting>
// <mediaobject><imageobject><imagedata fileref="ex-model-constant.png" format="PNG"/></imageobject><caption><para>Constant value model example</para></caption></mediaobject>
//
/** @brief Le système en représentation d'état le plus simple : une valeur constante.
 *
 *  <h3>Valeur constante</h3>
 *
 *  Cette fonction renvoie le système en représentation d'état le plus simple : une valeur constante,
 *  observée à travers un canal bruité, c'est-à-dire :
 *  @f[
 *  x_n = cst,\quad y_n = x_n+v_n
 *  @f]
 *  avec @f$v_n: \mathcal{N}(0,R)@f$
 *
 *  @param R: Variance du bruit d'observation
 *  @return Modèle linéaire
 *
 *  @par Exemple :
 *  @snippet exemples/src/ex-kalman.cc ex_modele_constante
 *  @image html ex-modele-constante.png width=800px
 *
 *  @sa modele_marche_aleatoire(), modele_kitagawa(), modele_imu()
 */
extern sptr<SSMLineaire> modele_constante(float R = 0.1);




// Kitagawa scalar non linear model
//
// Calling Sequence
// sys = model_kitagawa();
//
// Parameters
// sys: State / space model
//
// Description
// This scalar, non-linear model is frequently used to benchmark non-linear bayesian estimation methods.
//
// <refsection><title>Example</title></refsection>
// <programlisting>
//ssm = model_kitagawa();
//[x,y] = ssm_simu(ssm, 100);
//scf(0); plot_obs(ssm, y);
//scf(1); plot_states(ssm, x);
// </programlisting>
// <mediaobject><imageobject><imagedata fileref="ex-model-kitagawa1.png" format="PNG"/></imageobject><caption><para>Kitagawa model (observation sequence)</para></caption></mediaobject>
// <mediaobject><imageobject><imagedata fileref="ex-model-kitagawa2.png" format="PNG"/></imageobject><caption><para>Kitagawa model (state sequence)</para></caption></mediaobject>
//
// See also
//   ssm_simu
//   model_imu
//   model_rw
//   model_constant
//
// Bibliography
//  G. Kitagawa, Monte Carlo filter and smoother for non-Gaussian non-linear state space models, J. Comput. Graph. Statist., vol. 5, no. 1, pp. 1–25, 1996.
//

/** @brief Modèle non linéaire de Kitagawa.
 *
 * <h3>Modèle de Kitagawa</h3>
 * Ce modèle scalaire, non linéaire, est souvent utilisé pour comparer différentes méthodes d'estimation.
 *
 * @f[
 *     f_k(x) = x/2 + \frac{25x}{1+x^2}+8\cos\left(6k/5\right)
 * @f]
 *
 * @f[
 *     g(x) = x^2/20
 * @f]
 *
 *  @par Exemple :
 *  @snippet exemples/src/ex-kalman.cc ex_modele_kitagawa
 *  @image html modele-kitagawa-obs.png width=800px
 *  @image html modele-kitagawa-etat.png width=800px
 *
 *
 * @sa modele_marche_aleatoire(), modele_constante(), modele_imu()
 *
 * @par Références
 * G. Kitagawa, Monte Carlo filter and smoother for non-Gaussian non-linear state space models, J. Comput. Graph. Statist., vol. 5, no. 1, pp. 1–25, 1996.
 *
 */
extern sptr<SSM> modele_kitagawa();


/** @brief Résolution de l'équation de Riccati discrète (DARE).
 *
 * Résolution de l'équation de Riccati discrète (DARE) :
 * @f[
 *       P=A^{T}PA-(A^{T}PC)(R+C^{T}PC)^{-1}(C^{T}PA)+Q
 * @f]
 *
 *  @param A ...
 *  @param C ...
 *  @param Q ...
 *  @param R ...
 *  @param tolerance Erreur maximale avant de sortir
 *  @param iter_max Nombre maximal d'itérations
 *  @return Tuple de deux éléments : solution P, et erreur.
 *
 * D'après https://en.wikipedia.org/wiki/Algebraic_Riccati_equation#Solution
 *
 * @sa kalman_ssg()
 */
tuple<Tabf, float> dare(const Tabf &A,
          const Tabf &C, const Tabf &Q,
          const Tabf &R,
          float tolerance = 1e-7,
          entier iter_max = 1000);

// Compute the Kalman steady-state gain
//
// ssm: Linear state-space model (see <link linkend="ssm_init">ssm_init</link>)
// A: System matrix
// C: Observation matrix
// Q: System noise covariance matrix
// R: Observation noise covariance matrix
// K: Steady-state Kalman gain
//
// Description
// Compute the Kalman steady-state gain (for a linear system, with fixed matrix) using the Riccati equation solver of Scilab.
// The steady-state gain can be used to simplify the implementation of a Kalman filter (avoiding to compute any matrix inverse).
//
// <refsection><title>Example</title></refsection>
// In this example, we compute the steady-state gain for one of the most simple state / space model: a scalar random walk.
//
// Different values of process noise are tested, while the observation noise is kept constant.
// Kalman gain as a function of process noise (with observation noise kept constant).

//
/** @brief Calcul le gain de Kalman en régime établi
 *
 * <h3>Gain de Kalman en régime établi</h3>
 *
 * Calcule du gain en régime établi pour un système linéaire (matrices fixes), à l'aide de l'équation de Riccati discrète (DARE).
 * Ce gain peut être utilisé pour simplifié l'implémentation d'un filtre de Kalman
 * (il n'est plus nécessaire d'inverser des matrices).
 *
 * @param modele Modéle linéaire
 *
 * @par Exemple
 * Dans l'exemple ci-dessous, on calcule le gain en régime établi pour un système très simple,
 * à savoir une marche aléatoire scalaire. Différentes valeurs de bruit de process sont testées,
 * en gardant le bruit d'observation constant.
 * @snippet exemples/src/ex-kalman.cc ex_kalman_ssg
 * @image html ex-kalman-ssg.png "Steady-state Kalman gain vs process noise deviation" width=800px
 *
 * @sa dare()
 */
extern Tabf kalman_ssg(sptr<const SSMLineaire> modele);

// Plot a state sequence
//  ssm: State / space model
//  x: State sequence matrix
//  xtrue: State sequence matrix (true state)
//  xestim: State sequence matrix (estimated state)
// Plot the state evolution. The state vector must be a matrix
// with each column being the full state at a given time (e.g., the number of columns is the number of samples points).
// This function can also plot the true state agains the estimated state on the same graph.
/** @brief Affichage d'une séquence d'états
 *
 * <h3>Affichage d'une séquence d'états</h3>
 *
 * @param ssm Système
 * @param x Séquence d'états (en colonne : les différents états, en ligne : les différents instants)
 * @param xe Séquence d'états estimés
 * @return La figure générée
 *
 * @sa plot_obs()
 */
extern tsd::vue::Figures  plot_etats(sptr<const SSM> ssm, const Tabf &x, const Tabf &xe = Tabf());

// Plots an observations sequence
//  ssm: State / space model
//  y: Observation sequence (a matrix for which each column is an observation vector)
//
/** @brief Affichage d'une séquence d'observations
 *
 * <h3>Affichage d'une séquence d'observations</h3>
 *
 * @param sys Système
 * @param y   Séquence d'observations
 * @return La figure générée
 *
 * @sa plot_etats()
 */
tsd::vue::Figures  plot_obs(sptr<const SSM> sys, const Tabf &y);


//   <listitem><para><emphasis role="bold">Specification of states to compute:</emphasis></para></listitem>
// <itemizedlist>
//   <listitem><para><varname>cfg.estim_orient</varname>: If true, the model includes orientation (Euler angles) estimation. </para></listitem>
//   <listitem><para><varname>cfg.estim_pos</varname>:  </para></listitem>
// </itemizedlist>
//   <listitem><para><emphasis role="bold">Specification of available sensors:</emphasis></para></listitem>
// <itemizedlist>
//   <listitem><para><varname>cfg.acc_available</varname>: Boolean specifying if accelerometer measurements are available or not</para></listitem>
//   <listitem><para><varname>cfg.gyro_available</varname>: Boolean specifying if gyroscope measurements are available or not</para></listitem>
//   <listitem><para><varname>cfg.mag_available</varname>: Boolean specifying if magnetometer (compas) measurements are available or not</para></listitem>
//   <listitem><para><varname>cfg.baro_available</varname>: Boolean specifying if barometer measurements are available or not</para></listitem>
//   <listitem><para><varname>cfg.gps_available</varname>: Boolean specifying if GPS measurements are available or not</para></listitem>
// </itemizedlist>
//   <listitem><para><emphasis role="bold">Specification of sensors noise:</emphasis></para></listitem>
// <itemizedlist>
//   <listitem><para><varname>cfg.anoise</varname>: Noise standard deviation on the accelerometers outputs (scalar)</para></listitem>
//   <listitem><para><varname>cfg.gnoise</varname>: Noise standard deviation on the gyroscopes outputs (scalar)</para></listitem>
//   <listitem><para><varname>cfg.mnoise</varname>: Noise standard deviation on the compas outputs (scalar)</para></listitem>
//   <listitem><para><varname>cfg.bnoise</varname>: Noise standard deviation on the barometer output (scalar)</para></listitem>
//   <listitem><para><varname>cfg.gpsnoise</varname>: Noise standard deviation on the GPS outputs (scalar)</para></listitem>


/** @brief Configuration pour le modèle IMU */
struct IMUConfig
{
  /** @brief Accelerometer measurements are available or not. */
  bouléen acc_available;

  /** @brief Gyroscope measurements are available or not? */
  bouléen gyro_available;

  /** @brief Magnetometer (compas) measurements are available or not? */
  bouléen mag_available;

  /** @brief Barometer measurements are available or not? */
  bouléen baro_available;

  /** @brief GPS measurements are available or not? */
  bouléen gps_available;

  /** @brief If true, the model includes position estimation. */
  bouléen estim_pos;

  /** @brief If true, the model includes orientation (Euler angles) estimation. */
  bouléen estim_orient;

  /** @brief Sampling frequency (Hz) */
  float fs = 100;


  /** @brief Noise standard deviation on the accelerometers outputs (scalar) */
  float anoise   = 1.0;

  /** @brief Noise standard deviation on the gyroscopes outputs (scalar) */
  float gnoise   = 2.2;

  /** @brief Noise standard deviation on the compas outputs (scalar) */
  float mnoise   = 0.5;

  /** @brief Noise standard deviation on the barometer outputs (scalar) */
  float bnoise   = 1.0;

  /** @brief Noise standard deviation on the GPS outputs (scalar) */
  float gpsnoise = 5.0;
};

/** @brief Modèle générique IMU (Inertial Measuring Unit).
 *
 *  <h3>Modèle générique IMU (Inertial Measuring Unit)</h3>
 *
 * Generic IMU model with Euler angles parametrization, and with the following
 * (optionnal) sensors:
 *   - Accelerometer (3 axis),
 *   - Gyroscope (3 axis),
 *   - Magnetometer / compas (3 axis),
 *   - Barometer,
 *   - GPS
 *
 * This model can be used either for simulation (see @ref SSM::steps())
 * or for orientation / position estimation from experimental data
 * (for instance using an EKF filter, see @ref filtre_ekf()).
 *
 *  The model can use as a state vector either only the object orientation, or both the
 *  object orientation and localization.
 * So it is two parts, where each part is optionnal (according to the configuration):
 *  - <b>Position:</b> @f$(x,y,z,\dot{x},\dot{y},\dot{z},\ddot{x},\ddot{y},\ddot{z})@f$,
 *  - <b>Orientation:</b> @f$(\phi,\theta,\psi,\dot{\phi},\dot{\theta},\dot{\psi})@f$.
 *
 * The Observation vector is in this order: @f$(a_x,a_y,a_z,\omega_x,\omega_y,\omega_y,m_x,m_y,m_z,p_{baro},x_{gps},y_{gps},z_{gps})@f$,
 * where each element is optionnal, according to the configuration.
 *
 *  @param cfg IMU model configuration,
 *  @return  Modèle SSM
 */
extern sptr<SSM> modele_imu(const IMUConfig &cfg);

/** @} */



}
