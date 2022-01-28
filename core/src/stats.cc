#include "tsd/tsd.hpp"
#include "tsd/stats.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/figure.hpp"

using namespace tsd::fourier;
using namespace tsd::vue;


namespace tsd::stats {


  // D'après AN2197 "Implementing l-d on starcore"
  ArrayXf levinson_reel(const ArrayXf &r)
  {
    int n = r.rows();
    ArrayXf a(n);
    // R[i] - autocorrelation coefficients
    // A[i] - filter coefficients
    // K - reflection coefficients
    // Alpha - prediction gain
    a(0) = 1;
    float k = - r(1) / r(0);
    a(1) = k;
    float α = r(0) * (1-k*k);

    for(auto i = 2; i <= n - 1; i++)
    {
      float s = (r.segment(1,i-1) * a.segment(1, i-1).reverse()).sum() + r(i);
      k = -s / α;
      a.segment(1, i-1) += k * a.segment(1, i-1).reverse().eval();
      a(i) = k;
      α *= (1-k*k);
    }
    return a;
  }


  // D'après wikipédia
  // https://en.wikipedia.org/wiki/Levinson_recursion
  ArrayXf levinson(const ArrayXf &l1, const ArrayXf &c1, const ArrayXf &y)
  {
    int n = l1.rows();

    //msg("y = {}", y.transpose());

    tsd_assert_msg(n == c1.rows(), "Levinson : ligne et colonne doivent être de dimensions identiques.");
    tsd_assert_msg(n == y.rows(),  "Levinson : dimensions invalides.");
    tsd_assert_msg(l1(0) == c1(0), "Levinson : vecteurs invalides.");

    // t(n-1), t(n-2), ... , t(1)    (t(0) exclu)
    Eigen::VectorXf c1r = c1.tail(n-1).reverse().matrix();
    // t(-1), t(-2), ..., t(-(n-1))
    Eigen::VectorXf l1r = l1.tail(n-1).matrix();

    // Calcul des vecteurs avant et arrière
    // et de la solution
    Eigen::VectorXf f(n), b(n);
    f(0) = 1 / l1(0);
    b(0) = 1 / c1(0);

    Eigen::VectorXf x = Eigen::VectorXf::Zero(n);
    x(0) = y(0) / l1(0);

    // i = dimension nouveau vecteur
    for(auto i = 2; i <= n; i++)
    {

      Eigen::VectorXf fb(i), bb(i);
      fb.head(i-1)  = f.head(i-1);
      fb(i-1)       = 0;
      bb(0)         = 0;
      bb.tail(i-1)  = b.head(i-1);

      //float ef = M.row(i).head(i-1).dot(f);
      //float eb = M.row(0).segment(1,i-1).dot(b);

      float ef = c1r.tail(i-1).dot(f.head(i-1));
      float eb = l1r.head(i-1).dot(b.head(i-1));
      float fr = 1/(1 - eb * ef);

      // Mise à jour forward / backward
      f.head(i) = fr * fb - ef * fr * bb;
      b.head(i) = fr * bb - eb * fr * fb;

      float ex = c1r.tail(i-1).dot(x.head(i-1));

      x.head(i) += (y(i-1) - ex) * b.head(i);

#     if 0
      // Check des vecteurs f / b
      MatrixXf Ri(i, i);
      for(auto k = 0; k < i; k++)
      {
        for(auto l = 0; l < i; l++)
        {
          if(k >= l)
            Ri(k, l) = c1(k-l);
          else
            Ri(k, l) = l1(l-k);
        }
      }
      Eigen::VectorXf Rf = Ri * f.head(i);
      Eigen::VectorXf Rb = Ri * b.head(i);
      Eigen::VectorXf Rx = Ri * x.head(i);
      msg("Rf[{}] = {}", i, Rf.transpose());
      msg("Rb[{}] = {}", i, Rb.transpose());
      msg("Rx[{}] = {}", i, Rx.transpose());
#     endif
    }

    return x.array();
  }









  std::tuple<ArrayXf, ArrayXf> lpc(const ArrayXf &x, int p)
  {
    //ArrayXf r = xcorr(x, x, p).real();
    auto [lags, cr] = xcorr(x, x, p);
    ArrayXf a = levinson_reel(cr.tail(p).real());
    // Comment obtenir e ?
    // (1 - a(z)) * x = e
    // x_n = a_1 x_n-1 + a2 x_n-2 + ... + e_n
    // e est tout simplement l'erreur de prédiction

    // e = filtrage RIF des xi

    auto f = tsd::filtrage::filtre_rif<float,float>(a);

    ArrayXf e = f->step(x);

    return {a, e};
    //ArrayXf num = ArrayXf::Ones(1);
    //FRat<float>::rii(num, a);
  }

  SubSpaceSpectrum subspace_spectrum(const MatrixXcf &R, const SubSpaceSpectrumConfig &config)
  {
    int Ns = config.Ns;
    Eigen::SelfAdjointEigenSolver<MatrixXcf> eig(R.matrix());

    int m = R.rows(); // Nombre d'antennes en DOA

    if(Ns == -1)
    {
      // Dans l'ordre de valeur absolue croissante
      ArrayXcf l = eig.eigenvalues().array().reverse();
      ArrayXf la = l.abs();

      ArrayXf vmdl = ArrayXf::Zero(m-1);

      // Calcul du critère MDL pour déterminer le nombre de sources
      // D'après "Williams, D.B. “Detection: Determining the Number of Sources”
      // Digital Signal Processing Handbook
      for(auto Ns = 1; Ns < m; Ns++)
      {
        // N est le nombre d'observations en DOA
        int N = m;//n; // à confirmer
        // Produit des valeurs propres dans le bruit
        // (attention : tends vers zéro très vite ?)
        // -> Calcul de la somme des logs
        float somme = la.tail(m - Ns).log().sum();
        float mdl = -N * (somme - (Ns - m) * (-log(m-Ns) + somme)) + 0.5 * Ns * (2*m-Ns) * log(N);

        vmdl(Ns-1) = mdl;
      }

      //msg("mdl = {}", vmdl);

      vmdl.maxCoeff(&Ns);
      Ns++;

      if(config.debug_actif)
      {
        Figures f;
        f.subplot().plot(10 * la.log10(), "", "Valeurs propres (log)");
        auto s = f.subplot();
        s.plot(linspace(1, m-1, m-1), vmdl);
        s.titres("Critère MDL", "Nombre de sources", "MDL");
        //f.enregistrer("./res/st/mdl.png");
        f.afficher("MDL");
      }
    }

    tsd_assert_msg(Ns < m, "Le nombre de sources doit être inférieur à la dimension de la matrice de covariance (ici Ns={}, m={}).", Ns, m);

    //Nsm = 10;//nfreqs;
    //Nsm = nfreqs;

    msg("Nb sources estimée : {}", Ns);




    //ArrayXf fr = linspace(-0.5f, 0.5f - 1.0f / Nf, Nf);
    ArrayXXf P;
    ArrayXf Smu(config.Nf);//, Smu2(Nf);
    for(auto i = 0; i < config.Nf; i++)
    {
      // Calcul du "vecteur de steering" (ici simple exponentielle complexe)
      /*ArrayXcf*/
      auto [prm, v] = config.balayage(i, config.Nf, m);//sigexp(fr(i), m);

      if(i == 0)
        P.resize(config.Nf, prm.rows());

      P.row(i) = prm.transpose();

      float S = 0;
      for(auto j = 0; j < m - Ns; j++)
      {
        // MUSIC
        S  += std::norm(v.matrix().dot(eig.eigenvectors().col(j)));
        //S2 += (1.0f / eig.eigenvalues()(j)) * std::norm(v.matrix().dot(eig.eigenvectors().col(j)));
      }
      Smu(i) = 10 * log10(1.0f / S);
    }



    return {P, Smu, Ns};
  }

}


namespace tsd::fourier {
  std::tuple<ArrayXf, ArrayXf> psd_subspace(const ArrayXcf &x, int Ns, int Nf, int m)
  {
    tsd::stats::SubSpaceSpectrumConfig config;
    config.Ns = Ns;
    config.Nf = Nf;

    if(m == 0)
      m = x.rows();

    MatrixXcf R = tsd::stats::covmtx(x, m);
    auto r = tsd::stats::subspace_spectrum(R, config);
    return {r.var.col(0), r.spectrum};
  }


}
