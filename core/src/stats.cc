#include "tsd/tsd.hpp"
#include "tsd/stats.hpp"
#include "tsd/fourier.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"

#include "tsd/eig-util.hpp"
#include "Eigen/Eigenvalues"

using namespace tsd::fourier;
using namespace tsd::vue;


namespace tsd::stats {


  // D'après AN2197 "Implementing l-d on starcore"
  Vecf levinson_reel(const Vecf &r)
  {
    soit n = r.rows();
    Vecf a(n);
    // R[i] - autocorrelation coefficients
    // A[i] - filter coefficients
    // K - reflection coefficients
    // Alpha - prediction gain
    a(0) = 1;
    soit k = - r(1) / r(0);
    a(1) = k;
    soit α = r(0) * (1-k*k);

    pour(auto i = 2; i <= n - 1; i++)
    {
      soit s = (r.segment(1,i-1) * a.segment(1, i-1).reverse()).somme() + r(i);
      k = -s / α;
      a.segment(1, i-1) += k * a.segment(1, i-1).reverse();
      a(i) = k;
      α *= (1-k*k);
    }
    retourne a;
  }


  // D'après wikipédia
  // https://en.wikipedia.org/wiki/Levinson_recursion
  Vecf levinson(const Vecf &l1, const Vecf &c1, const Vecf &y)
  {
    soit n = l1.rows();

    tsd_assert_msg(n == c1.rows(), "Levinson : ligne et colonne doivent être de dimensions identiques.");
    tsd_assert_msg(n == y.rows(),  "Levinson : dimensions invalides.");
    tsd_assert_msg(l1(0) == c1(0), "Levinson : vecteurs invalides.");

    // t(n-1), t(n-2), ... , t(1)    (t(0) exclu)
    soit c1r = c1.tail(n-1).reverse();
    // t(-1), t(-2), ..., t(-(n-1))
    soit l1r = l1.tail(n-1);

    // Calcul des vecteurs avant et arrière
    // et de la solution
    Vecf f(n), b(n);
    f(0) = 1 / l1(0);
    b(0) = 1 / c1(0);

    Vecf x = Vecf::zeros(n);
    x(0) = y(0) / l1(0);

    // i = dimension nouveau vecteur
    pour(auto i = 2; i <= n; i++)
    {

      Vecf fb(i), bb(i);
      fb.head(i-1)  = f.head(i-1);
      fb(i-1)       = 0;
      bb(0)         = 0;
      bb.tail(i-1)  = b.head(i-1);

      soit ef = c1r.tail(i-1).dot(f.head(i-1));
      soit eb = l1r.head(i-1).dot(b.head(i-1));
      soit fr = 1/(1 - eb * ef);

      // Mise à jour forward / backward
      f.head(i) = fr * fb - ef * fr * bb;
      b.head(i) = fr * bb - eb * fr * fb;

      soit ex = c1r.tail(i-1).dot(x.head(i-1));

      x.head(i) += (y(i-1) - ex) * b.head(i);

#     if 0
      // Check des vecteurs f / b
      MatrixXf Ri(i, i);
      pour(auto k = 0; k < i; k++)
      {
        pour(auto l = 0; l < i; l++)
        {
          si(k >= l)
            Ri(k, l) = c1(k-l);
          sinon
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

    retourne x;
  }









  std::tuple<Vecf, Vecf> lpc(const Vecf &x, entier p)
  {
    //ArrayXf r = xcorr(x, x, p).real();
    soit [lags, cr] = xcorr(x, x, p);
    soit a = levinson_reel(real(cr.tail(p)));
    // Comment obtenir e ?
    // (1 - a(z)) * x = e
    // x_n = a_1 x_n-1 + a2 x_n-2 + ... + e_n
    // e est tout simplement l'erreur de prédiction

    // e = filtrage RIF des xi

    soit f = tsd::filtrage::filtre_rif<float,float>(a);

    soit e = f->step(x);

    retourne {a, e};
    //ArrayXf num = ArrayXf::Ones(1);
    //FRat<float>::rii(num, a);
  }

  SubSpaceSpectrum subspace_spectrum(const Tabcf &R, const SubSpaceSpectrumConfig &config)
  {
    soit Ns = config.Ns, m = R.rows(); // m = nombre d'antennes en DOA;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcf> eig(tab2etab(R).matrix());

    si(Ns == -1)
    {
      // Dans l'ordre de valeur absolue croissante
      Eigen::MatrixXcf eg = eig.eigenvalues().reverse();
      Veccf l    = etab2tab(eg);
      soit la   = abs(l);
      soit vmdl = Vecf::zeros(m-1);

      // Calcul du critère MDL pour déterminer le nombre de sources
      // D'après "Williams, D.B. “Detection: Determining the Number of Sources”
      // Digital Signal Processing Handbook
      pour(auto Ns = 1; Ns < m; Ns++)
      {
        // N est le nombre d'observations en DOA
        entier N = m;//n; // à confirmer
        // Produit des valeurs propres dans le bruit
        // (attention : tends vers zéro très vite ?)
        // -> Calcul de la somme des logs
        float somme = log(la.tail(m - Ns)).somme();
        float mdl = -N * (somme - (Ns - m) * (-log(m-Ns) + somme)) + 0.5 * Ns * (2*m-Ns) * log(N);

        vmdl(Ns-1) = mdl;
      }

      //msg("mdl = {}", vmdl);

      Ns = vmdl.index_max() + 1;

      si(config.debug_actif)
      {
        Figures f;
        f.subplot().plot(pow2db(la), "", "Valeurs propres (log)");
        soit s = f.subplot();
        s.plot(linspace(1, m-1, m-1), vmdl);
        s.titres("Critère MDL", "Nombre de sources", "MDL");
        //f.enregistrer("./res/st/mdl.png");
        f.afficher("MDL");
      }
    }

    tsd_assert_msg(Ns < m,
        "Le nombre de sources doit être inférieur à la dimension de la matrice de covariance (ici Ns={}, m={}).",
        Ns, m);

    msg("Nb sources estimée : {}", Ns);

    //ArrayXf fr = linspace(-0.5f, 0.5f - 1.0f / Nf, Nf);
    Tabf P;
    Vecf Smu(config.Nf);//, Smu2(Nf);
    pour(auto i = 0; i < config.Nf; i++)
    {
      // Calcul du "vecteur de steering" (ici simple exponentielle complexe)
      soit [prm, v] = config.balayage(i, config.Nf, m);//sigexp(fr(i), m);

      si(i == 0)
        P.resize(config.Nf, prm.rows());

      P.set_row(i, prm);

      float S = 0;
      pour(auto j = 0; j < m - Ns; j++)
      {
        // MUSIC
        S  += std::norm(vec2evec(v).matrix().dot(eig.eigenvectors().col(j)));
        //S2 += (1.0f / eig.eigenvalues()(j)) * std::norm(v.matrix().dot(eig.eigenvectors().col(j)));
      }
      Smu(i) = pow2db(1 / S);
    }



    retourne {P, Smu, Ns};
  }

}


namespace tsd::fourier {
  std::tuple<Vecf, Vecf> psd_sousesp(const Veccf &x, entier Ns, entier Nf, entier m)
  {
    tsd::stats::SubSpaceSpectrumConfig config;
    config.Ns = Ns;
    config.Nf = Nf;

    si(m == 0)
      m = x.rows();

    soit R = tsd::stats::covmtx(x, m);
    soit r = tsd::stats::subspace_spectrum(R, config);
    retourne {r.var.col(0), r.spectrum};
  }


}
