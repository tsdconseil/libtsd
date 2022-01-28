#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/fourier.hpp"
#include "tsd/figure.hpp"
#include <Eigen/QR>
#include <Eigen/SVD>
#include <set>

using namespace Eigen;
using namespace std;
using namespace tsd::vue;

namespace tsd::filtrage {




template<typename T>
Vecteur<T> rifamp_naif(const Vecteur<T> &h, int L)
{
  int N = h.rows();
  int M = (N - 1) / 2;
  tsd_assert_msg((N & 1) == 1, "rifamp_naif : nombre de coefficients impair attendu (N = {})", N);

  Vecteur<T> A(L);
  for(auto i = 0; i < L; i++)
  {
    double omega = (i * π) / L;
    A(i) = h(M);
    for(auto k = 0; k < M; k++)
      A(i) += 2 * h(k) * cos((M-k) * omega);
  }
  return A;
}


// TODO : template
std::tuple<ArrayXf, ArrayXf> rifamp(const ArrayXf &h, int L, bool symetrique)
{
  // Si N (nb coefs) > 2 * L (nb points de sortie désirés) :
  //  comment faire ?
  int N = h.rows();
  int M = (N-1)/2;
  int dim_fft = max(2*L, N);

  ArrayXcf h2 = ArrayXcf::Zero(dim_fft);
  h2.head(N) = h;

  ArrayXcf H = tsd::fourier::fft(h2);

  // select [0,π]
  H = H.head(dim_fft/2).eval();
  // frequency grid
  auto w = linspace(0,dim_fft/2-1,dim_fft/2) * π / (dim_fft / 2);

  // Corrige le facteur de normalisation de la fonction FFT
  H *= sqrt((float)dim_fft);
  H *= polar(M * w); // type I ou II

  // type III ou IV
  if(!symetrique)
    H *= -complex<float>(0,1);

  tsd_assert(H.rows() == (int) L);

  auto fr = linspace(0, 0.5 - (0.5 / L), L);

  return {fr, H.real()}; // discard zero imaginary part
}

//////////////////
/*
function maxind = locmax(x)
% function maxind = locmax(x)
% finds indices of local maxima of data x

x = x(:);
n = length(x);
maxind = find(x > [x(1)-1;x(1:n-1)] & x > [x(2:n);x(n)-1]);
*/


vector<int> locmax(IArrayXf x)
{
  auto n = x.rows();
  vector<int> res;
  for(auto i = 1; i + 1 < n; i++)
  {
    if((x(i) > x(i-1)) && (x(i) > x(i+1)))
      res.push_back(i);
  }
  return res;
}


vector<int> etap(const ArrayXd &E, const vector<int> &idx, int R, const set<int> &to_keep)
{
  // Ensuring The Alternation Property
  // E : error vector
  // v : index of original E

  // Retourne une liste d'index avec des extrema positifs et négatifs alternés

  vector<int> rmv;

  vector<int> v;
  auto xe = E(idx[0]);
  auto xv = idx[0];

  //ArrayXd Ep = E;
  //auto v1 = sousvec(Ep, idx);

  for(auto k = 1u; k < idx.size(); k++)
  {
    // On ne peut plus enlever de point intermédiaire
    if((int) idx.size() - (int) rmv.size() < R + 2)
    {
      v.push_back(xv);
      xe = E(idx[k]);
      xv = idx[k];
    }


    // Signe identique ?
    else if(copysign(1.0f, E(idx[k])) == copysign(1.0f, xe))
    {
      // Attention, il ne faut pas supprimer les extréma

      if(/*(to_keep.find(idx[k]) != to_keep.end()) || */abs(E(idx[k])) > abs(xe))
      {
        rmv.push_back(xv);
        xe = E(idx[k]);
        xv = idx[k];
      }
      else
      {
        rmv.push_back(idx[k]);
      }
    }
    else
    {
      v.push_back(xv);
      xe = E(idx[k]);
      xv = idx[k];
    }
  }
  v.push_back(xv);


  //auto v2 = sousvec(Ep, v);

  //msg("Erreur avant alt:\n{}", v1);
  //msg("Erreur après alt:\n{}", v2);


#  if 0
  //auto i = 0;
  //tsd_assert(rmv.size() >= R - v.size());
  while(v.size() < (unsigned int) R)
  {
    /*if(i >= (int) rmv.size())
    {
      erreur("Echec : idx orig = %d, R = %d, rmv.size() = %d", idx.size(), R, rmv.size());
      break;
    }
    v.push_back(rmv[i++]);*/
    if(rmv.size() == 0)
    {
      msg_erreur("Echec.");
      //res = -1;
      break;
    }
    v.push_back(rmv[0]);//i++]);
  }
# endif



  return v;
}



struct DesignRifEqConfig
{
  // Nb coefs
  int n;
  // Critère d'arrêt
  double ε = 1e-5;//1e-8
  ArrayXf D;
  ArrayXf W;
  bool debug;
};


// Implémentation d'après EL713 / FIR Chebychef
ArrayXf design_rif_eq(int n, IArrayXf D, IArrayXf W, bool debug)
{
  int res = 0;
  ArrayXd Wd = W.cast<double>();
  ArrayXd Dd = D.cast<double>();

  if((n & 1) == 0)
  {
    msg_erreur("design_rif_eq : ncoefs doit être impair (n = {}).", n);
    return ArrayXf();
  }

  // Critère d'arrêt
  const auto ε = 1e-5;//1e-8

  // Note : avec design_rif_freq : n = 2 * m - 1

  // n = 2m+1, m = (n-1)/2
  const int m = n / 2;

  // Nombre de points spécifiés
  const int l = W.rows();

  // r = dimension de l'espace de référence
  const int r = m + 2;

  if(D.rows() != l)
  {
    msg_erreur("design_rif_eq : D et W doivent avoir les mêmes dimensions (ici {} vs {}).", D.rows(), l);
    return ArrayXf();
  }

  if(debug)
    msg("rifeq : ncoefs: n={} (n=2m+1, avec m={}), dim esp réf : r={}, dim gabarit : l={}...", n, m, r, l);

  // Initialise l'ensemble de référence (à peu près équi-réparti, sur les points de poids non nul)
  vector<int> f = trouve(W > ε);
  set<int> bords;
  // Les intervalles de poids non nuls
  vector<pair<int, int>> bandes;

  {
    int debut = 0;
    for(auto i = 1; i + 1 < l; i++)
    {
      // Fin de la bande passante
      if((W(i) > ε) && (W(i+1) < ε))
      {
        tsd_assert(debut != -1);
        bandes.push_back({debut, i});
        bords.insert(i);
        debut = -1;
      }
      // Fin de la bande passante
      else if((W(i-1) < ε) && (W(i) > ε))
      {
        debut = i;
        bords.insert(i);
      }
    }
    tsd_assert(debut != -1);
    bandes.push_back({debut, l-1}); // Dernière bande
  }

  if(debug)
  {
    for(auto i = 0u; i < bandes.size(); i++)
    {
      auto &b = bandes[i];
      msg("  Bande {} : freq {} - {}, gain moyen = {}", i, (1.0 * b.first) / l, (1.0 * b.second) / l,
          D.segment(b.first, b.second - b.first + 1).mean());
    }
  }


  // Il faut se débrouiller pour avoir le début et la fin de la bande coupée
  vector<int> k(r);
  int offset = 0;
  for(auto [a,b] : bandes)
  {
    auto n = b - a + 1;
    auto ratio = (1.0f * n) / l;

    int cnt = (int) round(ratio * r);

    if(b == l-1)
      cnt = r - offset;

    if(cnt + offset > r)
      cnt = r - offset;
    if(cnt > 0)
    {
      auto ls = linspace(a, b, cnt);
      for(auto i = 0; i < cnt; i++)
        k[offset+i] = round(ls(i));
      offset += cnt;
    }
  }

  ArrayXd h = ArrayXd::Zero(n);

  tsd_assert((int) n == 2 * m + 1);

  double δ = 0;

  auto itr = 0;

  if(debug)
  {
    Figures f;
    f.subplot().plot(W, "-b", "Poids");
    f.subplot().plot(D, "-b", "Réponse voulue");
    f.afficher("Remez spec");
  }

  //Figures fig, fig2, fig3;

  Eigen::VectorXd w  = linspace(0,l-1,l).cast<double>() * π / l;
  Eigen::MatrixXd mm = linspace(0,m,m+1).transpose().cast<double>();

  ArrayXd s(r);
  for(auto i = 0; i < r; i++)
  {
    if((i & 1) == 0)
      s(i) = 1;
    else
      s(i) = -1;
  }

  auto nitr_max = 100;

  ArrayXd err;
  ArrayXf verr = ArrayXf::Zero(nitr_max);

  for(itr = 0; itr < nitr_max; itr++)
  {
    Figures fig;

    tsd_assert((int) k.size() == r);

    // Construit la matrice d'interpolation
    // W(k) : R+1x1, m : M+1x1
    Eigen::MatrixXd Mat(r, m+2);
    Eigen::VectorXd wk = sousvec(w,k);
    Mat.leftCols(m+1) = (wk * mm).array().cos();
    ArrayXd Wdk = sousvec(Wd,k);
    Mat.rightCols(1)  = s / Wdk; // TO CHECK

    Eigen::MatrixXd Dk = sousvec(Dd,k);
    ArrayXd x = Mat.colPivHouseholderQr().solve(Dk).array();

    tsd_assert(x.rows() == (int) (m+2));

    // Sorte de fftshift
    h.head(m) = x.segment(1, m).reverse();
    h(m)      = 2 * x(0);
    h.tail(m) = x.segment(1, m);
    h /= 2;

    if(debug)// && (itr < nr*nc))
    {
      fig.subplot().plot(h.cast<float>(), "", "Filtre interpolé (coefs)");
    }

    // Dernier terme = E = Erreur au niveau de chaque point utilisé
    // Note : probablement très inférieur aux erreurs intermédiaires, sauf
    // si on arrive à la convergence
    δ = abs(x(m+1));

    //ArrayXf A = rifamp(h.cast<float>(), l);
    ArrayXd A = rifamp_naif(h, l);

    // Erreur pondérée
    err = (A - Dd) * Wd;
    // Maximum de l'erreur
    auto emax_general = err.abs().maxCoeff();

    {
      // Si convergence, emax_general ~ δ
      auto err_max_rel = (emax_general-δ) / δ;

      verr(itr) = emax_general;

      if(debug)
        msg(" errmax={:.3g}, δ={:.3g}, rel={:.3g}", emax_general, δ, err_max_rel);
      if(err_max_rel < 1e-3)
      {
        if(debug)
          msg(" convergence ok : sortie.");
        break;
      }
    }

    // Mise à jour de l'ensemble des points d'interpolation
    vector<int> k2;
    ArrayXd ae = err.abs();

    if(W(0) > ε)
      k2.push_back(0);

    for(auto i = 1; i + 1 < l; i++)
    {
      // Fin de la bande passante
      if((W(i) > ε) && (W(i+1) < ε))
        k2.push_back(i);
      // Fin de la bande passante
      else if((W(i-1) < ε) && (W(i) > ε))
        k2.push_back(i);
      else if(W(i) < ε)
        continue;
      // Maximum local
      else if((ae(i) > ae(i-1)) && (ae(i) > ae(i+1)) && (ae(i) >= δ*0.9))
        k2.push_back(i);
    }

    if(W(l-1) > ε)
      k2.push_back(l - 1);

    /*while((int) k2.size() < m + 2)
    {

    }*/

    if(debug)// && (itr < nr*nc))
    {
      //auto s = fig2.subplot();
      //s.plot(err.cast<float>());
      //for(auto i : k2)
        //s.canva().marqueur({(float) i, (float) err(i)}, Marqueur::CERCLE, 9);
      fig.subplot().plot(A.cast<float>(), "", "A");
    }


    if(debug)
    {
      msg("Itération n°={}, δ={} (err nv pts interp) , erreur générale max={}, nb max locaux (avant vérif altern) : {}, R={}",
          itr, δ, emax_general, k2.size(), r);
      //Figures fs;
      auto f = fig.subplot();
      f.plot(err.cast<float>(), "", "Erreur");
      for(auto i : k2)
        f.canva().marqueur({(float) i, (float) err(i)}, Marqueur::CERCLE, 9);

      f.canva().set_couleur(Couleur::Rouge);
      for(auto i : k)
        f.canva().marqueur({(float) i, (float) err(i)}, Marqueur::CERCLE, 5);


      f = fig.subplot();
      f.plot((err.abs() /*+ 1e-20f*/).log().cast<float>(), "", "Erreur (log)");


      for(auto i : k2)
        f.canva().marqueur({(float) i, (float) log(abs(err(i)))}, Marqueur::CERCLE, 9);

      f.canva().set_couleur(Couleur::Rouge);
      for(auto i : k)
        f.canva().marqueur({(float) i, (float) log(abs(err(i)))}, Marqueur::CERCLE, 5);

      f = fig.subplot();
      f.plot(W, "", "Pondération");

      fig.afficher(fmt::format("Itération {}", itr));
    }

    if(δ == 0)
    {
      msg_erreur("δ = 0");
      res = -1;
      break;
    }

    if((int) k2.size() < m + 2)
    {
      res = -1;
      msg_erreur("Pas assez de maxima trouvés. Maxima : {} éléments, m = {}.", k2.size(), m);
      //k2.push_back(l-1);
      break;
    }

    auto oldk = k2;

    // ensure the alternation property
    k2 = etap(err, k2, r, bords);

    if(debug)
    {
      // if newk is too large, remove points until size is correct
      msg("apres vérif alternance : nb max = {}, R = {}", k2.size(), r);
    }

    if((int) k2.size() < m + 2)
    {
      res = -1;
      msg_erreur("Pas assez de maxima trouvés après vérif alternance. Maxima : {} éléments, m = {}.", k2.size(), m);
      break;
    }

    // Tant qu'on peut enlever 2 points,
    // on enlève un couple contigue avec erreur totale minimale
    while((int) k2.size() >= r + 2)
    {
      int N = k2.size();
      float best = 1e10;
      int bestj = 0;


      ArrayXd ae = err.abs();
      for(auto j = 1; j + 1 < N; j++)
      {
        //if((ae(k2[j]) < best) && ((ae(k2[j+1]) < best)))
        if((ae(k2[j])  + ae(k2[j+1]) < best))
        {
          best = ae(k2[j]) + ae(k2[j+1]);
          //best = min(ae(k2[j]), ae(k2[j+1]));
          bestj = j;
        }
      }

      if(ae(k2[N-1]) < best)
      {
        k2.erase(k2.begin() + (N - 1));
      }
      else if(bestj == 0)
      {
        k2.erase(k2.begin());
      }
      else
      {
        k2.erase(k2.begin() + bestj);
        k2.erase(k2.begin() + bestj);
      }
    }

    //
    while((int) k2.size() > r)
    {
      // il faut enlever 1 point
      if(abs(err(k2[0])) < abs(err(k2[k2.size()-1])))
        k2.erase(k2.begin());
      else
        k2.pop_back();
    }

    if(debug)// && (itr < 8))
    {
      for(auto i : k2)
        fig.gcf().canva().marqueur({(float) i, (float) err(i)}, Marqueur::ETOILE, 9);
    }

    if(debug)
      msg("newk : {} éléments.", k2.size());

    // --------------- Check Convergence -------------------------------
    ArrayXd tmp1 = sousvec(err,k2);
    auto emax = tmp1.abs().maxCoeff();
    auto err_max = (emax-δ) / δ;

    verr(itr) = emax;

    if(debug)
      msg(" errmax = {:.3g}, δ = {:.3g}, rel = {:.3g}", emax, δ, err_max);
    if(err_max < 1e-3)//epsilon)
    {
      if(debug)
        msg(" convergence ok : sortie.");
      break;
    }

    if((int) k2.size() < m + 2)
    {
      res = -1;
      msg_erreur(" maxima : {} éléments, m = {}.", k2.size(), m);
      //cout << "err(oldk) = " << sousvec(err,oldk) << "\n";
      //cout << "err(newk) = " << sousvec(err,newk) << "\n";
      break;
    }

    k = k2;
  }

  if(debug)
  {
    /*fig.afficher("Remez debug - h (coefficients)");
    fig2.afficher("Remez debug - Erreur");
    fig3.afficher("Remez debug - A (réponse fréquentielle)");*/

    {
      Figure f;
      f.plot(verr, "", "Erreur = f(itr)");
      f.afficher("Evo erreur");
    }

    {
      Figure f;
      f.plot(err.cast<float>());
      f.afficher("Erreur finale");
    }

  }

  if(debug)
    msg("fin rifeq, itr={}, δ={:e}", itr, δ);

  if(res)
    return ArrayXf();

  return h.cast<float>();
}

/*template<typename D1, typename D2>
auto lsq(const Eigen::MatrixBase<D1> &A, const Eigen::MatrixBase<D2> &b)
{
  //Eigen::VectorXf x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  //Eigen::VectorXf x = A.colPivHouseholderQr().solve(b);
  //return x;
  return A.colPivHouseholderQr().solve(b);
}*/


template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> lsq(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &b)
{
  //Eigen::VectorXf x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  //Eigen::VectorXf x = A.colPivHouseholderQr().solve(b);
  //return x;
  return A.colPivHouseholderQr().solve(b);
}


template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> lsi(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &b,
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &w,
    bool debug = false)
{
  int n = A.rows();//, m = A.cols();

  tsd_assert(n == b.rows());

  using M = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
  using V = Eigen::Matrix<T, Eigen::Dynamic, 1>;

  M Ai = A;
  V bi = b;
  V x = lsq(A, b);

  V p = V::Ones(n);

  if(debug)
  {
    Figures f;
    f.subplot().plot(b.array().template cast<float>(), "g-", "Vecteur cible");
    f.subplot().plot(w.array().template cast<float>(), "a-", "Pondération");
    f.afficher("LSI - Vecteur cible");
  }

  for(auto itr = 0; itr < 100; itr++)
  {
    V e = A * x - b;
    float e2 = e.norm();

    V Ea = (e.cwiseAbs().array() * w.array()).matrix();

    float emax = Ea.maxCoeff();
    float emin = Ea.minCoeff();

    //msg("Vecteur d'errreur : {}", e.transpose());

    if(debug)// && ((itr % 4) == 0))
    {
      msg("Itération {} : norme 2 = {}, norme inf = {}, emin = {}", itr, e2, emax, emin);
      //ArrayXf tmp = Ea;//e.array().abs();
      Figures f;
      f.subplot().plot((A * x).array().template cast<float>(), "-g", "Ax");
      f.subplot().plot(Ea.array().template cast<float>(), "-r", "Erreur pondérée");
      f.subplot().plot(p.array().template cast<float>(), "-b", "Pondération L2");
      f.afficher(fmt::format("LSI - ITR {}",  itr));
    }

    if((emax - emin) / emax < 1e-5)
    {
      msg("emax - emin = {} -> exit.", emax - emin);
      break;
    }

    double α = 0.5;
    p = (p.array() * (1 + (α / (emax - emin)) * (Ea.array() - emin))).matrix().eval();

    /*for(auto i = 0; i < n; i++)
    {
      float l;
      float eai = Ea(i);//std::abs(e(i));

      float alpha = 10;//0.4;
      //l = 1 + alpha * (eai - emin) / (emax - emin);

      //p(i) += alpha * (eai - emin) / (emax - emin);

      p(i) *= 1 + alpha * (eai - emin) / (emax - emin);
    }*/

    p /= p.sum();

    Ai = p.asDiagonal() * A;
    bi = (b.array() * p.array()).matrix();

    /*for(auto i = 0; i < n; i++)
    {
      Ai.row(i) = A.row(i) * p(i);
      bi(i)     = b(i) * p(i);
    }*/

    x = lsq(Ai, bi);
  }
  msg("Fin lsi.");
  return x;
}

#if 0
Eigen::VectorXf lsi(const Eigen::MatrixXf &A,
                    const Eigen::VectorXf &b,
                    const ArrayXf &w,
                    bool debug = false)
{
  int n = A.rows();//, m = A.cols();

  tsd_assert(n == b.rows());

  Eigen::MatrixXf Ai = A;
  Eigen::VectorXf bi = b;

  Eigen::VectorXf x = lsq(A, b);

  ArrayXf p = ArrayXf::Ones(n);

  if(debug)
  {
    Figures f;
    f.subplot().plot(b.array(), "g-", "Vecteur cible");
    f.subplot().plot(w, "a-", "Pondération");
    f.afficher("LSI - Vecteur cible");
  }

  for(auto itr = 0; itr < 100; itr++)
  {
    Eigen::VectorXf e = A * x - b;
    float e2 = e.norm();

    ArrayXf Ea = e.cwiseAbs().array() * w;

    float emax = Ea.maxCoeff();
    float emin = Ea.minCoeff();

    //msg("Vecteur d'errreur : {}", e.transpose());

    if(debug)// && ((itr % 4) == 0))
    {
      msg("Itération {} : norme 2 = {}, norme inf = {}, emin = {}", itr, e2, emax, emin);
      //ArrayXf tmp = Ea;//e.array().abs();
      Figures f;
      f.subplot().plot((A * x).array(), "-g", "Ax");
      f.subplot().plot(Ea, "-r", "Erreur pondérée");
      f.subplot().plot(p, "-b", "Pondération L2");
      f.afficher(fmt::format("LSI - ITR {}",  itr));
    }

    if((emax - emin) / emax < 1e-5)
    {
      msg("emax - emin = {} -> exit.", emax - emin);
      break;
    }

    float α = 0.5;
    p *= 1 + (α / emax - emin) * (Ea - emin);

    /*for(auto i = 0; i < n; i++)
    {
      float l;
      float eai = Ea(i);//std::abs(e(i));

      float alpha = 10;//0.4;
      //l = 1 + alpha * (eai - emin) / (emax - emin);

      //p(i) += alpha * (eai - emin) / (emax - emin);

      p(i) *= 1 + alpha * (eai - emin) / (emax - emin);
    }*/

    p /= p.sum();

    Ai = p.matrix().asDiagonal() * A;
    bi = (b.array() * p).matrix();

    /*for(auto i = 0; i < n; i++)
    {
      Ai.row(i) = A.row(i) * p(i);
      bi(i)     = b(i) * p(i);
    }*/

    x = lsq(Ai, bi);
  }
  msg("Fin lsi.");
  return x;
}
#endif



ArrayXf design_rif_eq2(int nc, IArrayXf D, IArrayXf W, bool debug)
{
  int n = D.rows();
  int m = (nc+1) / 2;

  msg("Design rif eq2 : nc = {}, n = {}, m = {}", nc, n, m);

  //tsd_assert((nc % 2) == 0); // Pour l'instant
  tsd_assert(W.rows() == n);

# if 0
  Eigen::MatrixXd A(n, m);
  //Eigen::VectorXf b(n);

  for(auto i = 0; i < n; i++)
  {
    double Omega = (i * π) / n;
    for(auto k = 0; k < m; k++)
    {
      if((k == m-1) && (nc & 1))
        A(i, k) = 1;
      else
      // Type I (nc impair)
      //if((nc % 2))
        A(i, k) = 2 * cos(Omega * (k - (nc-1)/2.0f));
      //else
        // Type II (nc pair)
        //A(i, k) = 2 * cos(Omega * (k - (nc-1)/2.0f));
    }
  }


  Eigen::VectorXd Dd = D.matrix().cast<double>();
  Eigen::VectorXd Wd = W.matrix().cast<double>();

  ArrayXf x = lsi(A, Dd, Wd, debug).array().cast<float>();
# endif


  Eigen::MatrixXf A(n, m);
  //Eigen::VectorXf b(n);

  for(auto i = 0; i < n; i++)
  {
    double Omega = (i * π) / n;
    for(auto k = 0; k < m; k++)
    {
      if((k == m-1) && (nc & 1))
        A(i, k) = 1;
      else
      // Type I (nc impair)
      //if((nc % 2))
        A(i, k) = 2 * std::cos(Omega * (k - (nc-1)/2.0));
      //else
        // Type II (nc pair)
        //A(i, k) = 2 * cos(Omega * (k - (nc-1)/2.0f));
    }
  }


  Eigen::VectorXf Dd = D.matrix();
  Eigen::VectorXf Wd = W.matrix();

  ArrayXf x = lsi(A, Dd, Wd, debug).array();

  ArrayXf h(nc);

  h.head(nc/2) = x.head(nc/2);
  h.tail(nc/2) = x.head(nc/2).reverse();
  // Filtre de type I, coefficient central
  if(nc & 1)
    h(nc/2) = x(nc/2);

  msg("ok.");

  return h;
}

ArrayXf design_rif_eq(int nc, const vector<SpecFreqIntervalle> &spec, bool debug, bool cheby_mode)
{
  float ε = 1e-10;
  auto l = 16 * nc;
  // Poids nul par défaut
  ArrayXf W = ArrayXf::Zero(l) + ε,
          D = ArrayXf::Ones(l);

  msg("rifeq : nb intervalles de fréq : {}", nc);

  for(auto i = 0u; i < spec.size(); i++)
  {
    auto &s = spec[i];
    auto ib = (int) round(2 * s.fb * l);
    auto ih = (int) round(2 * s.fh * l);

    ib = min(ib, l-1);
    ih = min(ih, l-1);

    msg("Intervalle : [{} ({}) - {} ({})], atten = {}, poids = {}", s.fb, ib, s.fh, ih, s.atten, s.poids);

    D.segment(ib, (ih - ib) + 1) = s.atten;
    W.segment(ib, (ih - ib) + 1) = s.poids + ε;

    if(i + 1 < spec.size())
    {
      auto ib2 = (unsigned int) round(2 * spec[i+1].fb * l);
      auto cnt = (ib2 - ih) + 1;
      D.segment(ih, cnt) = linspace(s.atten, spec[i+1].atten, cnt);
    }
  }
  if(cheby_mode)
    return design_rif_eq(nc, D, W, debug);
  else
    return design_rif_eq2(nc, D, W, debug);
}




}






