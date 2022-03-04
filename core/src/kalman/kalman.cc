#include "tsd/tsd.hpp"
#include "tsd/apps/kalman.hpp"
#include "tsd/vue.hpp"
#include <Eigen/Dense>



namespace tsd::kalman {


  /** D'après https://en.wikipedia.org/wiki/Algebraic_Riccati_equation#Solution */
  std::tuple<MatrixXf, float> dare(const MatrixXf &A,
      const MatrixXf &C, const MatrixXf &Q,
      const MatrixXf &R,
      float tolerance,
      int iter_max)
  {
    MatrixXf P = Q,
      At = A.transpose(),
      Ct = C.transpose(),
      Rinv = R.inverse();
    float diff = 0;
    for(auto i = 0; i < iter_max; ++i)
    {
      MatrixXf P2 = At * P * A - At * P * C * (R + Ct * P * C).inverse() * Ct * P * A + Q;
      diff = ((P2 - P).array().abs().maxCoeff());
      P = P2;
      if(diff < tolerance)
        break;
    }
    return {P, diff};
  }


sptr<SSMLineaire> ssm_lineaire(
    const Eigen::MatrixXf &A,
    const Eigen::MatrixXf &C,
    const Eigen::MatrixXf &B,
    const Eigen::MatrixXf &D)
{
  auto res = std::make_shared<SSMLineaire>();
  res->A = A;
  res->C = C;
  res->B = B;
  res->D = D;
  res->ns = A.rows();
  res->ni = B.cols();
  res->no = C.rows();

  tsd_assert_msg(A.cols() == A.rows(), "SSM linéaire : la matrice A devrait être carrée.");
  tsd_assert_msg(A.cols() == C.cols(), "SSM linéaire / matrice C : nombre de colonnes invalide.");

  return res;
}

MatrixXf jacobienne_num(std::function<VectorXf(const VectorXf &)> f, const VectorXf &x)
{
  VectorXf y = f(x);
  int ne = x.rows(), ns = y.rows();
  MatrixXf J(ns, ne);

  float ε = 1e-5;

  for(auto i = 0; i < ne; i++)
  {
    VectorXf x2 = x;
    x2(i) += ε;
    VectorXf y2 = f(x2);

    J.col(i) = (y2 - y) / ε;
  }

  return J;
}

std::string SSM::lis_nom_etat(int i) const
{
  if(i < (int) noms_etats.size())
    return noms_etats[i];
  else
    return fmt::format("Etat {}", i);
}

std::string SSM::lis_nom_obs(int i) const
{
  if(i < (int) noms_obs.size())
    return noms_obs[i];
  else
    return fmt::format("Obs {}", i);
}

MatrixXf SSM::Jf(const VectorXf &x) const
{
  return jacobienne_num(std::bind(&SSM::f, this, std::placeholders::_1), x);
}

MatrixXf SSM::Jg(const VectorXf &x) const
{
  return jacobienne_num(std::bind(&SSM::g, this, std::placeholders::_1), x);
}

VectorXf SSMLineaire::f(const VectorXf &x) const
{
  return A * x; // u ?
}

VectorXf SSMLineaire::g(const VectorXf &x) const
{
  return C * x;
}

MatrixXf SSMLineaire::Jf(const VectorXf &x) const
{
  return A;
}

MatrixXf SSMLineaire::Jg(const VectorXf &x) const
{
  return C;
}

std::tuple<ArrayXXf, ArrayXXf> SSM::steps(int n, const VectorXf &x0, const ArrayXXf &vin_)
{
  if(LQ.rows() == 0)
  {
    LQ = Q;
    LR = R;
    if(!Q.isZero())
      LQ = Q.llt().matrixL();
    if(!R.isZero())
      LR = R.llt().matrixL();

    msg("R = {}", R);
    msg("LR = {}", LR);
  }


  MatrixXf x = MatrixXf::Zero(ns, n);
  MatrixXf y = MatrixXf::Zero(no, n);

  if(x0.rows() > 0)
  {
    if(x0.rows() != ns)
      echec("ssm_simu: nb states = {}, but got x0 of size {}.", ns, x0.rows());
    etat = x0;
  }
  if(etat.rows() == 0)
    etat = VectorXf::Zero(ns);
  x.col(0) = etat;

  MatrixXf vin = MatrixXf::Zero(ni, n);

  if(vin_.size() > 0)
  {
    if((vin_.rows() != ni) || (vin_.cols() != n))
      echec("ssm_simu: nb input = {}, nb spl = {}, but got vin of size {}*{}",
          ni, n, vin.rows(), vin.cols());
    vin = vin_.matrix();
  }

  auto lmodele = dynamic_cast<SSMLineaire *>(this);


  if(lmodele)
  {
    if(ni > 0)
    {
      y.col(0) = lmodele->C * x.col(0) + lmodele->D * vin.col(0);

      for(auto i = 1; i < n; i++)
      {
        auto V = LQ * randn(ns).matrix();
        x.col(i) = lmodele->A * x.col(i-1) + lmodele->B * vin.col(i) + V;
        auto W = LR * randn(no).matrix();
        y.col(i) = lmodele->C * x.col(i) + lmodele->D * vin.col(i)  + W;
      }
    }
    else
    {
      y.col(0) = lmodele->C * x.col(0);

      for(auto i = 1; i < n; i++)
      {
        VectorXf V = LQ * randn(ns).matrix();
        x.col(i) = lmodele->A * x.col(i-1) + V;
        VectorXf W = LR * randn(no).matrix();
        y.col(i) = lmodele->C * x.col(i)  + W;
      }
    }
  }
  else
  {
    //k = 0;
    //y.col(0) = g(x.col(0));

    for(auto i = 0; i < n; i++)
    {
      VectorXf W = LR * randn(no).matrix();
      y.col(i) = g(etat) + W;
      VectorXf V = LQ * randn(ns).matrix();
      etat = f(etat) + V;

      if(i + 1 < n)
        x.col(i+1) = etat;

      k++;
    }
  }
  return {x.array(), y.array()};
}

void SSM::verifications() const
{
  msg("Vérification système SSM [{}]...", nom);
  msg("ns={}, ni={}, no={}", ns, ni, no);

  VectorXf x = randn(ns).matrix();

  {
    // (1) Vérif jacobienne de g (fonction d'observation)
    MatrixXf MJg = Jg(x);

    if(MJg.rows() != no)
      echec("Dimension Jg incohérente ({} != {}).", MJg.rows(), no);
    if(MJg.cols() != ns)
      echec("Dimension Jg incohérente ({} != {}).", MJg.cols(), ns);


    MatrixXf Jgnum = jacobienne_num(std::bind(&SSM::g, this, std::placeholders::_1), x);//list(ekf.ssm.g, ekf.ssm), x);

    float eg   = (MJg - Jgnum).norm();
    msg("Erreur jacobienne g : {}.", eg);

    if(eg > 1e-4)
    {
      int m, k;
      (MJg - Jgnum).array().abs().maxCoeff(&m, &k);
      msg("Erreur max sur Jg: {}, {}", m, k);
    }
  }

  {
    MatrixXf MJf = Jf(x);

    if(MJf.rows() != ns)
      echec("Dimension Jf incohérente.");
    if(MJf.cols() != ns)
      echec("Dimension Jf incohérente.");

    MatrixXf Jfnum = jacobienne_num(std::bind(&SSM::f, this, std::placeholders::_1), x);

    float ef = (MJf - Jfnum).norm();
    msg("Erreur jacobienne f : {}.", ef);

    if(ef > 1e-4)
    {
      int m, k;
      (MJf - Jfnum).array().abs().maxCoeff(&m, &k);
      msg("Erreur max sur Jf: {}, {}", m, k);
    }
  }

  msg("ok.");
}



struct FiltreEKF : FiltreSSM
{
  VectorXf vobs;
  MatrixXf P, K;
  int k;
  bool mode_debug = false;

  FiltreEKF(sptr<SSM> ssm)
  {
    this->ssm = ssm;

    ssm->verifications();

    etat = VectorXf::Zero(ssm->ns);
    // Matrice de covariance de l'état initial
    P = MatrixXf::Identity(ssm->ns,ssm->ns);

    vobs = VectorXf::Zero(ssm->no);
    k = 0;
    //ekf.ssm.k = 0;
    mode_debug = false;

    //if(argn(2) > 1)
    //    ekf.state = varargin(1);
    //end;

    //if(argn(2) > 2)
    //    ekf.P = varargin(2);
    //end;

    /*if(~isequal(size(ekf.state), [ssm.ns, 1]))
        error(sprintf("ekf_init: invalid initial state vector size (should be a %d * 1 vector, and is %d * %d.", ssm.ns, size(ekf.state,1), size(ekf.state,2)));
    end;

    if(~isequal(size(ekf.P), [ssm.ns, ssm.ns]))
        error("ekf_init: invalid initial state covariance matrix size.");
    end;*/
  }

  /*
   *     nspl = size(vobs,2);

      // Appel direct sur une matrice d'observations
      if(nspl > 1)

          s = zeros(ekf.ssm.ns,nspl);
          err = zeros(ekf.ssm.ns,nspl);
          //t00 = tic();
          for i = 2:nspl
              vobs_ = vobs(:,i);
              if(argn(2) > 2)
                 dts = varargin(1);
                [ekf,s_] = ekf_process(ekf,vobs_,dts(i));
              else
                [ekf,s_] = ekf_process(ekf,vobs_);
              end;
              s(:,i) = s_;
              //err(:,i) = xv(:,i) - s_;
              if(modulo(i,64) == 0)
                  //if(toc() - t00 > 1)
                    //  t00 = tic();
                      printf("ekf_process: itération %d/%d.\n", i, nspl);
                  //end;
              end;
          end
          return;
      end;
   */


  //function [ekf,s] = ekf_process(ekf,vobs,varargin)
  VectorXf step(const VectorXf &y, const VectorXf &u)//const VectorXf vobs, float dt)
  {
  // Proceed to one or more steps of the extended Kalman filter
  //
  // Calling Sequence
  // [ekf,x] = ekf_process(ekf,y)
  // [ekf,x] = ekf_process(ekf,y,dt)
  //
  // Parameters
  // ekf: EKF filter object (see <link linkend="ekf_init">ekf_init</link>)
  // y: Observations vector or matrix (one column for each time step)
  // dt: optional time interval(s) (scalar  or vector)
  // x: Output state vector or matrix (one column for each time step)
  //
  // Description
  // With the Extended Kalman Filter (EKF), the Kalman equations can be applied to a non-linear model (suboptimal solution).
  // To do this, the system matrix <latex>$A$</latex> (matrix governing system evolution) and <latex>$C$</latex> (matrix governing system observations) are taken as the local derivative of the system and observation functions.
  //
  // If the user has not provided these Jacobians in <link linkend="ekf_init">ekf_init</link>, then numerical derivatives are computed (using the Scilab build-in function <link linkend="numderivative">numderivative</link>).
  //
  // See also
  // ekf_init
  // ekf_check
  //
  // Authors
  //  J.A., full documentation available on <ulink url="http://www.tsdconseil.fr/log/bayes">http://www.tsdconseil.fr/log/bayes</ulink>

      //if(ekf.ssm.subtype == 'linear')
      //    [ekf,s] = kalman_process(ekf,vobs);
      //    return;
      //end;



      // Mise à jour intervalle de temps ?
      // (pour échantillonnage non régulier)
      //if(argn(2) > 2)
      //    ekf.ssm.dt = varargin(1);
      //end;

      // Jacobienne de la fonction de transition (à partir de l'état précédent)
      //MatrixXf Ak;
      //if(!ssm.Jf)
      //  Ak = jacobienne_num(ssm.f, state);
      //else
      MatrixXf Ak = ssm->Jf(etat);


      // Estimation a priori de l'état suivant
      VectorXf xm = ssm->f(etat);

      if(mode_debug)
        msg("etat predit : {}", xm);

      // Jacobienne de la fonction d'observation
      MatrixXf Ck = ssm->Jg(xm);
      //if(ssm.Jg)
      //   Ck = ssm.Jg(xm);
      //else
      //   Ck = jacobienne_num(ssm.g, xm);

      if(mode_debug)
        msg("Ck : {}", Ck);

      // Observation a priori (à partir du nouvel état !)
      VectorXf vobs_th = ssm->g(xm);

      if(mode_debug)
      {
        msg("obs predit : {}", vobs_th);
        msg("obs reelle : {}", vobs);
      }

      // A priori estimation of the error covariance matrix
      MatrixXf Pm = Ak * P * Ak.transpose() + ssm->Q;

      if(mode_debug)
        msg("P predit : {}", Pm);

      // Optimal gain of the filter
      if((Pm.isZero()) && (ssm->R.isZero()))
          K = Ck.transpose().lu().inverse();
      else
        K = (Pm*Ck.transpose()) * (Ck * Pm * Ck.transpose() + ssm->R).lu().inverse();

      if(mode_debug)
        msg("K : {}", K);

      // A posteriori estimation of the state variables
      VectorXf x = xm + K * (y - vobs_th);

      if(mode_debug)
        msg("etat post : {}", x);

      etat = x;

      // A posteriori estimation of the error covariance matrix:
      P = (MatrixXf::Identity(ssm->ns,ssm->ns)-K*Ck)*Pm;

      if(mode_debug)
      {
        msg("K*Ck :   {}", K*Ck);
        msg("Gauche : {}", MatrixXf::Identity(ssm->ns,ssm->ns)-K*Ck);
        msg("P post : {}", P);
      }
      k++;

      return x;
  }

};


struct FiltreKalman : FiltreSSM
{
  MatrixXf P;


  // Initialization of a standard Kalman filter object
  FiltreKalman(sptr<SSMLineaire> &ssm, const Eigen::VectorXf &x0, const Eigen::MatrixXf &p0)
  {
    this->ssm = ssm;
    if(x0.rows() == 0)
      etat = Eigen::VectorXf::Zero(ssm->ns);
    else
      etat = x0;
    // Matrice de covariance de l'état initial
    if(p0.rows() == 0)
      P = Eigen::MatrixXf::Identity(ssm->ns,ssm->ns);
    else
      P = p0;
    tsd_assert_msg(etat.rows() == ssm->ns,
        "Filtre de Kalman : nombre de lignes incorrect pour le vecteur état initial ({}, alors que ns={}).", etat.rows(), ssm->ns);
    tsd_assert_msg((P.rows() == ssm->ns) && (P.cols() == ssm->ns),
        "Filtre de Kalman : la matrice de covariance sur l'état initial devrait être de dimension {}x{}.", ssm->ns, ssm->ns);
  }


  // y: Current observations (column vector or matrix with one column for each time step)
  // u: Current inputs (column vector or matrix with one column for each time step)
  // x0: Initial state estimate (default is 0 or previous state if kalman_process already called)
  // p0: Covariance matrix on the initial state (default is eye matrix, or previously estimated covariance matrix is kalman_process already called)
  // retourne x : Current estimate of the state (column vector or matrix with one column for each time step)
  VectorXf step(const VectorXf &y, const VectorXf &u)
  {
    auto sl = std::dynamic_pointer_cast<SSMLineaire>(ssm);

    // Estimation a priori de l'état suivant
    VectorXf xm = sl->A * etat;

    if(sl->B.rows() > 0)
      xm += sl->B * u;

    // Observation a priori
    VectorXf vobs_th = sl->C * xm;

    if(sl->D.rows())
      vobs_th += sl->D * u;

    // A priori estimation of the error covariance matrix
    MatrixXf Pm = sl->A *  P * sl->A.transpose() + sl->Q;

    MatrixXf K;

    // Optimal gain of the filter
    if((Pm.isZero()) && (sl->R.isZero()))
        K = sl->C.inverse();
    else
      K = (Pm * sl->C.transpose()) * (sl->C * Pm * sl->C.transpose() + sl->R).lu().inverse();

    // A posteriori estimation of the state variables
    VectorXf x = xm + K * (y - vobs_th);

    etat = x;

    // A posteriori estimation of the error covariance matrix:
    P = (Eigen::MatrixXf::Identity(sl->ns,sl->ns) - K * sl->C) * Pm;

    return etat.array();
  }
};


ArrayXXf FiltreSSM::steps(const ArrayXXf &y, const ArrayXXf &u)
{
  auto nspl = y.cols();


  if(u.rows() != ssm->ni)
  {
    msg_erreur("kalman_process: nb inputs should be {}, got {}.", ssm->ni, u.rows());
    return ArrayXXf();
  }
  if((ssm->ni > 0) && (u.cols() != nspl))
  {
    msg_erreur("kalman_process: not as many input columns as observation columns ({} vs {}).",
        u.cols(), nspl);
    return ArrayXXf();
  }

  MatrixXf s = MatrixXf::Zero(ssm->ns, nspl);
  s.col(0) = etat;
  msg("Kalman : traitement {} échantillons...", nspl);
  if(ssm->ni > 0)
  {
    for(auto i = 1; i < nspl; i++)
      s.col(i) = step(y.col(i).matrix(), u.col(i).matrix());
  }
  else
  {
    for(auto i = 1; i < nspl; i++)
      s.col(i) = step(y.col(i).matrix());
  }
  return s;
}

sptr<FiltreSSM> filtre_kalman(sptr<SSMLineaire> ssm, const Eigen::VectorXf &x0, const Eigen::MatrixXf &p0)
{
  return std::make_shared<FiltreKalman>(ssm, x0, p0);
}

sptr<FiltreSSM> filtre_ekf(sptr<SSM> ssm)
{
  return std::make_shared<FiltreEKF>(ssm);
}



sptr<SSMLineaire> modele_marche_aleatoire(int ndim, const MatrixXf &Q, const MatrixXf &R)
{
  auto ssm = std::make_shared<SSMLineaire>();

  ssm->ni = 0;
  ssm->no = ndim;
  ssm->ns = ndim;

  ssm->Q = Q;
  ssm->R = R;

  if(Q.rows() == 0)
    ssm->Q = MatrixXf::Identity(ndim, ndim);
  if(R.rows() == 0)
    ssm->R = MatrixXf::Identity(ndim, ndim);


  // Vecteur d'état : x
  ssm->A = MatrixXf::Identity(ndim,ndim);
  ssm->C = MatrixXf::Identity(ndim,ndim);

  //ssm = ssm_init(A,B,C,D,Q,R);
  ssm->nom = "marche aléatoire";

  return ssm;
}


sptr<SSMLineaire> modele_constante(float R)
{
  auto res = modele_marche_aleatoire(1, MatrixXf::Zero(1,1), R * MatrixXf::Ones(1,1));
  res->nom = "valeur constante";
  return res;
}


struct ModeleKitagawa: SSM
{
  VectorXf f(const VectorXf &x_) const
  {
    auto x = x_.array();
    return (x/2 + (25*x)/(1+x*x)+8*cos(1.2*k)).matrix();
  }
  VectorXf g(const VectorXf &x) const
  {
    VectorXf y(1);
    y(0) = x(0) * x(0) / 20;
    return y;
  }

  ModeleKitagawa()
  {
    ni = 0;
    ns = 1;
    no = 1;
    Q  = 10 * MatrixXf::Ones(1,1);
    R  = MatrixXf::Ones(1,1);
    nom = "kitagawa";
  }
};

sptr<SSM> modele_kitagawa()
{
  return std::make_shared<ModeleKitagawa>();
}




MatrixXf kalman_ssg(sptr<const SSMLineaire> modele)
{


    // Bruit de process null
    // => Gain asymptotique nul
    if(modele->Q.isZero())
    {
      msg_avert("kalman_ssg: process noise is null => zero steady state kalman gain.");
      return MatrixXf::Zero(modele->ns, modele->no);
    }
    else
    {
      auto [P, eps] = dare(modele->A.transpose(), modele->C.transpose(), modele->Q, modele->R, 1e-6, 1000);
      //return (modele->R + modele->C.transpose()*P*modele->C).inverse() * modele->C.transpose() * P * modele->A;
      return P * modele->C.transpose() * (modele->R + modele->C*P*modele->C.transpose()).inverse();
    }
}




tsd::vue::Figures plot_etats(sptr<const SSM> ssm, const ArrayXXf &x, const ArrayXXf &xe)
{
  tsd::vue::Figures res;
  bool has_xe = xe.rows() > 0;

  for(auto i = 0; i < ssm->ns; i++)
  {
    auto f = res.subplot();
    if(has_xe)
    {
      f.plot(x.row(i),"bo-", "Vrai");
      f.plot(xe.row(i),"gs-", "Estimé");
    }
    else
      f.plot(x.row(i),"bo-");
    f.titre(ssm->lis_nom_etat(i));
    if(has_xe)
    {
      auto f = res.subplot();
      f.plot(xe.row(i)-x.row(i),"r-","Erreur");
    }
  }
  return res;
}


tsd::vue::Figures plot_obs(sptr<const SSM> sys, const ArrayXXf &y)
{
  tsd::vue::Figures res;
  //f.titre("Observations");
  for(auto i = 0; i < sys->no; i++)
  {
    auto f = res.subplot();
    f.plot(y.row(i),"b-o");
    f.titre(sys->lis_nom_obs(i));
  }
  return res;
}





}
