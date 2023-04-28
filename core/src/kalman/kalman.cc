#include "tsd/tsd.hpp"
#include "tsd/apps/kalman.hpp"
#include "tsd/vue.hpp"
#include "tsd/eig-util.hpp"
#include <Eigen/Dense>



namespace tsd::kalman {


  /** D'après https://en.wikipedia.org/wiki/Algebraic_Riccati_equation#Solution */
  tuple<Tabf, float> dare(const Tabf &A,
      const Tabf &C, const Tabf &Q,
      const Tabf &R,
      float tolerance,
      entier iter_max)
  {
    Tabf P = Q,
        At    = A.transpose(),
        Ct    = C.transpose(),
        Rinv  = R.inverse_matricielle();
    float diff = 0;
    pour(auto i = 0; i < iter_max; ++i)
    {
      soit P2 = At * P * A - At * P * C * (R + Ct * P * C).inverse_matricielle() * Ct * P * A + Q;
      diff = (abs(P2 - P)).maxCoeff();
      P = P2;
      si(diff < tolerance)
        break;
    }
    retourne {P, diff};
  }


sptr<SSMLineaire> ssm_lineaire(
    const Tabf &A,
    const Tabf &C,
    const Tabf &B,
    const Tabf &D)
{
  soit res = std::make_shared<SSMLineaire>();
  res->A = A;
  res->C = C;
  res->B = B;
  res->D = D;
  res->ns = A.rows();
  res->ni = B.cols();
  res->no = C.rows();

  tsd_assert_msg(A.cols() == A.rows(), "SSM linéaire : la matrice A devrait être carrée.");
  tsd_assert_msg(A.cols() == C.cols(), "SSM linéaire / matrice C : nombre de colonnes invalide.");

  if(!(B.est_vide()))
  {
    tsd_assert_msg(C.rows() == D.rows(), "SSM linéaire : les matrices C et D devraientt avoir le même nombre de lignes.");
  }

  retourne res;
}

Tabf jacobienne_num(std::function<Vecf(const Vecf &)> f, const Vecf &x)
{
  soit y = f(x);
  soit ne = x.rows(), ns = y.rows();
  soit ε = 1e-5f;
  Tabf J(ns, ne);

  pour(auto i = 0; i < ne; i++)
  {
    soit x2 = x.clone();
    x2(i) += ε;
    soit y2 = f(x2);
    soit d = (y2 - y) / ε;
    J.col(i) = d;
  }

  retourne J;
}

string SSM::lis_nom_etat(entier i) const
{
  si(i < (entier) noms_etats.size())
    retourne noms_etats[i];
  sinon
    retourne sformat("Etat {}", i);
}

string SSM::lis_nom_obs(entier i) const
{
  si(i < (entier) noms_obs.size())
    retourne noms_obs[i];
  sinon
    retourne sformat("Obs {}", i);
}

Tabf SSM::Jf(const Vecf &x) const
{
  retourne jacobienne_num(std::bind(&SSM::f, this, std::placeholders::_1), x);
}

Tabf SSM::Jg(const Vecf &x) const
{
  retourne jacobienne_num(std::bind(&SSM::g, this, std::placeholders::_1), x);
}



Vecf SSMLineaire::f(const Vecf &x) const
{
  retourne A.matprod(x);
}

Vecf SSMLineaire::g(const Vecf &x) const
{
  retourne C.matprod(x);
}

Tabf SSMLineaire::Jf(const Vecf &x) const
{
  retourne A;
}

Tabf SSMLineaire::Jg(const Vecf &x) const
{
  retourne C;
}

tuple<Tabf, Tabf> SSM::steps(entier n, const Vecf &x0, const Tabf &vin_)
{
  si(LQ.rows() == 0)
  {
    LQ = Q.clone();
    LR = R.clone();
    si(!Q.est_nul())
      LQ = etab2tab(Eigen::MatrixXf(tab2etab(Q).llt().matrixL()));
    si(!R.est_nul())
      LR = etab2tab(Eigen::MatrixXf(tab2etab(R).llt().matrixL()));
  }


  soit x = Tabf::zeros(ns, n),
       y = Tabf::zeros(no, n);

  si(x0.rows() > 0)
  {
    si(x0.rows() != ns)
      echec("ssm_simu: nb states = {}, but got x0 of size {}.", ns, x0.rows());
    etat = x0;
  }
  si(etat.rows() == 0)
    etat = Vecf::zeros(ns);
  x.col(0) = etat;

  soit vin = Tabf::zeros(ni, n);
  si(vin_.rows() > 0)
  {
    si((vin_.rows() != ni) || (vin_.cols() != n))
      echec("ssm_simu: nb input = {}, nb spl = {}, but got vin of size {}*{}",
          ni, n, vin.rows(), vin.cols());
    vin = vin_;
  }

  soit lmodele = dynamic_cast<SSMLineaire *>(this);


  si(lmodele)
  {
    si(ni > 0)
    {
      y.col(0) = (lmodele->C).matprod(x.col(0));// +
      si(!(lmodele->D.est_vide()))
        y.col(0) += (lmodele->D).matprod(vin.col(0));

      pour(auto i = 1; i < n; i++)
      {
        soit V = LQ.matprod(randn(ns));
        x.col(i) = lmodele->A.matprod(x.col(i-1)) + V;

        si(!(lmodele->B.est_vide()))
          x.col(i) += (lmodele->B).matprod(vin.col(i));

        soit W = LR.matprod(randn(no));
        y.col(i) = lmodele->C.matprod(x.col(i)) + W;

        si(!(lmodele->D.est_vide()))
          y.col(i) += (lmodele->D).matprod(vin.col(i));
      }
    }
    sinon
    {
      y.col(0) = lmodele->C.matprod(x.col(0));

      pour(auto i = 1; i < n; i++)
      {
        soit V = LQ.matprod(randn(ns));
        x.col(i) = lmodele->A.matprod(x.col(i-1)) + V;
        soit W = LR.matprod(randn(no));
        y.col(i) = lmodele->C.matprod(x.col(i))  + W;
      }
    }
  }
  sinon
  {
    //k = 0;
    //y.col(0) = g(x.col(0));

    pour(auto i = 0; i < n; i++)
    {
      soit W = LR.matprod(randn(no));
      y.col(i) = g(etat) + W;
      soit V = LQ.matprod(randn(ns));
      etat = f(etat) + V;

      si(i + 1 < n)
        x.col(i+1) = etat;

      k++;
    }
  }
  retourne {x, y};
}

void SSM::verifications() const
{
  msg("Vérification système SSM [{}]...", nom);
  msg("ns={}, ni={}, no={}", ns, ni, no);

  soit x = randn(ns);

  {
    // (1) Vérif jacobienne de g (fonction d'observation)
    soit MJg = Jg(x);

    si(MJg.rows() != no)
      echec("Dimension Jg incohérente ({} != {}).", MJg.rows(), no);
    si(MJg.cols() != ns)
      echec("Dimension Jg incohérente ({} != {}).", MJg.cols(), ns);


    soit Jgnum = jacobienne_num(std::bind(&SSM::g, this, std::placeholders::_1), x);

    soit eg   = (MJg - Jgnum).norme();
    msg("Erreur jacobienne g : {}.", eg);

    si(eg > 1e-4)
    {
      soit [m, k] = (abs(MJg - Jgnum)).index_max();
      msg("Erreur max sur Jg: {}, {}", m, k);
    }
  }

  {
    soit MJf = Jf(x);

    si(MJf.rows() != ns)
      echec("Dimension Jf incohérente.");
    si(MJf.cols() != ns)
      echec("Dimension Jf incohérente.");

    soit Jfnum = jacobienne_num(std::bind(&SSM::f, this, std::placeholders::_1), x);

    float ef = (MJf - Jfnum).norme();
    msg("Erreur jacobienne f : {}.", ef);

    si(ef > 1e-4)
    {
      soit [m, k] = (abs(MJf - Jfnum)).index_max();
      msg("Erreur max sur Jf: {}, {}", m, k);
    }
  }

  msg("ok.");
}



struct FiltreEKF : FiltreSSM
{
  Vecf vobs;
  Tabf P, K;
  entier k;
  bouléen mode_debug = non;

  FiltreEKF(sptr<SSM> ssm)
  {
    this->ssm = ssm;

    ssm->verifications();

    etat = Vecf::zeros(ssm->ns);
    // Matrice de covariance de l'état initial
    P = Tabf::eye(ssm->ns);

    vobs = Vecf::zeros(ssm->no);
    k = 0;
    //ekf.ssm.k = 0;
    mode_debug = non;

    //si(argn(2) > 1)
    //    ekf.state = varargin(1);
    //end;

    //si(argn(2) > 2)
    //    ekf.P = varargin(2);
    //end;

    /*si(~isequal(size(ekf.state), [ssm.ns, 1]))
        error(sprintf("ekf_init: invalid initial state vector size (should be a %d * 1 vector, and is %d * %d.", ssm.ns, size(ekf.state,1), size(ekf.state,2)));
    end;

    si(~isequal(size(ekf.P), [ssm.ns, ssm.ns]))
        error("ekf_init: invalid initial state covariance matrix size.");
    end;*/
  }

  /*
   *     nspl = size(vobs,2);

      // Appel direct sur une matrice d'observations
      si(nspl > 1)

          s = zeros(ekf.ssm.ns,nspl);
          err = zeros(ekf.ssm.ns,nspl);
          //t00 = tic();
          pour i = 2:nspl
              vobs_ = vobs(:,i);
              si(argn(2) > 2)
                 dts = varargin(1);
                [ekf,s_] = ekf_process(ekf,vobs_,dts(i));
              sinon
                [ekf,s_] = ekf_process(ekf,vobs_);
              end;
              s(:,i) = s_;
              //err(:,i) = xv(:,i) - s_;
              si(modulo(i,64) == 0)
                  //si(toc() - t00 > 1)
                    //  t00 = tic();
                      printf("ekf_process: itération %d/%d.\n", i, nspl);
                  //end;
              end;
          end
          retourne;
      end;
   */


  //function [ekf,s] = ekf_process(ekf,vobs,varargin)
  Vecf step(const Vecf &y, const Vecf &u)//const VectorXf vobs, float dt)
  {
  // Proceed to one or more steps of the extended Kalman filter
  //
  // Calling Sequence
  // [ekf,x] = ekf_process(ekf,y)
  // [ekf,x] = ekf_process(ekf,y,dt)
  //
  // Parameters
  // ekf: EKF filter object (see <link linkend="ekf_init">ekf_init</link>)
  // y: Observations vector or matrix (one column pour each time step)
  // dt: optional time interval(s) (scalar  or vector)
  // x: Output state vector or matrix (one column pour each time step)
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

      //si(ekf.ssm.subtype == 'linear')
      //    [ekf,s] = kalman_process(ekf,vobs);
      //    retourne;
      //end;



      // Mise à jour intervalle de temps ?
      // (pour échantillonnage non régulier)
      //si(argn(2) > 2)
      //    ekf.ssm.dt = varargin(1);
      //end;

      // Jacobienne de la fonction de transition (à partir de l'état précédent)
      //MatrixXf Ak;
      //si(!ssm.Jf)
      //  Ak = jacobienne_num(ssm.f, state);
      //sinon
      soit Ak = ssm->Jf(etat);


      // Estimation a priori de l'état suivant
      soit xm = ssm->f(etat);

      si(mode_debug)
        msg("etat predit : {}", xm);

      // Jacobienne de la fonction d'observation
      soit Ck = ssm->Jg(xm);
      //si(ssm.Jg)
      //   Ck = ssm.Jg(xm);
      //sinon
      //   Ck = jacobienne_num(ssm.g, xm);

      si(mode_debug)
        msg("Ck : {}", Ck);

      // Observation a priori (à partir du nouvel état !)
      soit vobs_th = ssm->g(xm);

      si(mode_debug)
      {
        msg("obs predit : {}", vobs_th);
        msg("obs reelle : {}", vobs);
      }

      // A priori estimation of the error covariance matrix
      soit Pm = Ak.matprod(P.matprod(Ak.transpose())) + ssm->Q;

      si(mode_debug)
        msg("P predit : {}", Pm);

      // Optimal gain of the filter
      si((Pm.est_nul()) && (ssm->R.est_nul()))
          K = Ck.transpose().inverse_matricielle();
      sinon
        K = (Pm.matprod(Ck.transpose())).matprod((Ck.matprod((Pm.matprod(Ck.transpose())) + ssm->R).inverse_matricielle()));

      si(mode_debug)
        msg("K : {}", K);

      // A posteriori estimation of the state variables
      soit x = xm + K.matprod((y - vobs_th));

      si(mode_debug)
        msg("etat post : {}", x);

      etat = x;

      // A posteriori estimation of the error covariance matrix:
      P = (Tabf::eye(ssm->ns) - K.matprod(Ck)).matprod(Pm);

      si(mode_debug)
      {
        //msg("K*Ck :   {}", K*Ck);
        //msg("Gauche : {}", MatrixXf::Identity(ssm->ns,ssm->ns)-K*Ck);
        //msg("P post : {}", P);
      }
      k++;

      retourne x;
  }

};


struct FiltreKalman : FiltreSSM
{
  Tabf P;


  // Initialization of a standard Kalman filter object
  FiltreKalman(sptr<SSMLineaire> &ssm, const Vecf &x0, const Tabf &p0)
  {
    this->ssm = ssm;
    si(x0.rows() == 0)
      etat = Vecf::zeros(ssm->ns);
    sinon
      etat = x0;
    // Matrice de covariance de l'état initial
    si(p0.rows() == 0)
      P = Tabf::eye(ssm->ns);
    sinon
      P = p0;
    tsd_assert_msg(etat.rows() == ssm->ns,
        "Filtre de Kalman : nombre de lignes incorrect pour le vecteur état initial ({}, alors que ns={}).", etat.rows(), ssm->ns);
    tsd_assert_msg((P.rows() == ssm->ns) && (P.cols() == ssm->ns),
        "Filtre de Kalman : la matrice de covariance sur l'état initial devrait être de dimension {}x{}.", ssm->ns, ssm->ns);
  }


  // y: Current observations (column vector or matrix with one column pour each time step)
  // u: Current inputs (column vector or matrix with one column pour each time step)
  // x0: Initial state estimate (default is 0 or previous state si kalman_process already called)
  // p0: Covariance matrix on the initial state (default is eye matrix, or previously estimated covariance matrix is kalman_process already called)
  // retourne x : Current estimate of the state (column vector or matrix with one column pour each time step)
  Vecf step(const Vecf &y, const Vecf &u)
  {
    soit sl = std::dynamic_pointer_cast<SSMLineaire>(ssm);

    // Estimation a priori de l'état suivant
    soit xm = sl->A.matprod(etat);

    si(sl->B.rows() > 0)
      xm += sl->B.matprod(u);

    // Observation a priori
    soit vobs_th = sl->C.matprod(xm);

    si(sl->D.rows())
      vobs_th += sl->D.matprod(u);

    // A priori estimation of the error covariance matrix
    soit Pm = sl->A.matprod(P.matprod((sl->A.transpose()))) + sl->Q;

    Tabf K;

    // Optimal gain of the filter
    si((Pm.est_nul()) && (sl->R.est_nul()))
        K = sl->C.inverse_matricielle();
    sinon
      K = (Pm.matprod(sl->C.transpose())).matprod((sl->C.matprod(Pm.matprod(sl->C.transpose())) + sl->R).inverse_matricielle());

    // A posteriori estimation of the state variables
    soit x = xm + K.matprod(y - vobs_th);

    etat = x;

    // A posteriori estimation of the error covariance matrix:
    P = (Tabf::eye(sl->ns) - K.matprod(sl->C)).matprod(Pm);

    retourne etat;
  }
};


Tabf FiltreSSM::steps(const Tabf &y, const Tabf &u)
{
  soit nspl = y.cols();


  si(u.rows() != ssm->ni)
  {
    msg_erreur("kalman_process: nb inputs should be {}, got {}.", ssm->ni, u.rows());
    retourne {};
  }
  si((ssm->ni > 0) && (u.cols() != nspl))
  {
    msg_erreur("kalman_process: not as many input columns as observation columns ({} vs {}).",
        u.cols(), nspl);
    retourne {};
  }

  soit s = Tabf::zeros(ssm->ns, nspl);
  s.col(0) = etat;
  msg("Kalman : traitement {} échantillons...", nspl);
  si(ssm->ni > 0)
  {
    pour(auto i = 1; i < nspl; i++)
      s.col(i) = step(y.col(i), u.col(i));
  }
  sinon
  {
    pour(auto i = 1; i < nspl; i++)
      s.col(i) = step(y.col(i));
  }
  retourne s;
}

sptr<FiltreSSM> filtre_kalman(sptr<SSMLineaire> ssm, const Vecf &x0, const Tabf &p0)
{
  retourne std::make_shared<FiltreKalman>(ssm, x0, p0);
}

sptr<FiltreSSM> filtre_ekf(sptr<SSM> ssm)
{
  retourne std::make_shared<FiltreEKF>(ssm);
}



sptr<SSMLineaire> modele_marche_aleatoire(entier ndim, const Tabf &Q, const Tabf &R)
{
  soit ssm = std::make_shared<SSMLineaire>();

  ssm->ni = 0;
  ssm->no = ndim;
  ssm->ns = ndim;

  ssm->Q = Q;
  ssm->R = R;

  si(Q.rows() == 0)
    ssm->Q = Tabf::eye(ndim);
  si(R.rows() == 0)
    ssm->R = Tabf::eye(ndim);


  // Vecteur d'état : x
  ssm->A = Tabf::eye(ndim);
  ssm->C = Tabf::eye(ndim);
  ssm->nom = "marche aléatoire";

  retourne ssm;
}


sptr<SSMLineaire> modele_constante(float R)
{
  soit res = modele_marche_aleatoire(1, Tabf::zeros(1,1), R * Tabf::ones(1,1));
  res->nom = "valeur constante";
  retourne res;
}


struct ModeleKitagawa: SSM
{
  Vecf f(const Vecf &x) const
  {
    retourne (x/2 + (25*x)/(1+x*x)+8*cos(1.2*k));
  }
  Vecf g(const Vecf &x) const
  {
    Vecf y(1);
    y(0) = x(0) * x(0) / 20;
    retourne y;
  }

  ModeleKitagawa()
  {
    ni = 0;
    ns = 1;
    no = 1;
    Q  = 10 * Tabf::ones(1,1);
    R  = Tabf::ones(1,1);
    nom = "kitagawa";
  }
};

sptr<SSM> modele_kitagawa()
{
  retourne std::make_shared<ModeleKitagawa>();
}




Tabf kalman_ssg(sptr<const SSMLineaire> modele)
{


    // Bruit de process null
    // => Gain asymptotique nul
    si(modele->Q.est_nul())
    {
      msg_avert("kalman_ssg: process noise is null => zero steady state kalman gain.");
      retourne Tabf::zeros(modele->ns, modele->no);
    }
    sinon
    {
      soit [P, eps] = dare(modele->A.transpose(), modele->C.transpose(), modele->Q, modele->R, 1e-6, 1000);
      retourne P.matprod(modele->C.transpose().matprod(modele->R + modele->C.matprod(P.matprod(modele->C.transpose()))).inverse_matricielle());
    }
}




tsd::vue::Figures plot_etats(sptr<const SSM> ssm, const Tabf &x, const Tabf &xe)
{
  tsd::vue::Figures res;
  bouléen has_xe = xe.rows() > 0;

  pour(auto i = 0; i < ssm->ns; i++)
  {
    soit f = res.subplot();
    si(has_xe)
    {
      f.plot(x.row(i),"bo-", "Vrai");
      f.plot(xe.row(i),"gs-", "Estimé");
    }
    sinon
      f.plot(x.row(i),"bo-");
    f.titre(ssm->lis_nom_etat(i));
    si(has_xe)
    {
      soit f = res.subplot();
      f.plot(xe.row(i)-x.row(i),"r-","Erreur");
    }
  }
  retourne res;
}


tsd::vue::Figures plot_obs(sptr<const SSM> sys, const Tabf &y)
{
  tsd::vue::Figures res;
  //f.titre("Observations");
  pour(auto i = 0; i < sys->no; i++)
  {
    soit f = res.subplot();
    f.plot(y.row(i),"b-o");
    f.titre(sys->lis_nom_obs(i));
  }
  retourne res;
}





}
