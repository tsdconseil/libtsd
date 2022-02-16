#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/telecom.hpp"
#include "tsd/figure.hpp"
#include <Eigen/QR>
#include <cctype>

using namespace tsd::vue;
using namespace std;

namespace tsd::telecom {


#if 0
//template<typename T, typename Tc>
//Eigen::Ref<Vecteur<T>> convol
//  (const Eigen::Ref<const Vecteur<Tc>> h, const Eigen::Ref<const Vecteur<T>> x)

template<typename T, typename Tc>
  auto convol(const Eigen::ArrayBase<Tc> &h, const Eigen::ArrayBase<T> &x)
{
  auto f = tsd::filtrage::filtre_rif<T, Tc>(h);
  return f->step(x);
}
#endif


auto convol(IArrayXf h, IArrayXf x)
{
  auto f = tsd::filtrage::filtre_rif<float, float>(h);
  return f->step(x);
}

auto convol(IArrayXf h, IArrayXcf x)
{
  auto f = tsd::filtrage::filtre_rif<float, cfloat>(h);
  return f->step(x);
}



struct EgaliseurRIF: FiltreGen<cfloat>
{
  sptr<FormeOnde> fo;
  int K = 1;
  float α = 0.01;
  int N1 = 11, N2 = 11, cnt = 0;

  ArrayXf h_avant, h_arriere;
  ArrayXcf wnd, wnd_deci;

  bool init_ok = false;

  enum Structure
  {
    FFE, DFE
  } structure = FFE;

  enum Errf
  {
    DEC, CMA
  } errf = DEC;


  EgaliseurRIF(sptr<FormeOnde> fo,
      const string &s_structure, const string &s_errf,
      int K, float α, int N1, int N2)
  {
    init_ok  = false;
    this->K  = K;
    this->α  = α;
    this->N1 = N1;
    this->N2 = N2;
    this->fo  = fo;

    tsd_assert_msg(fo, "Egaliseur RIF : la forme d'onde doit être spécifiée.");


    string s_errf2 = s_errf;
    for(auto &c: s_errf2)
      c = toupper(c);
    string s_structure2 = s_structure;
    for(auto &c: s_structure2)
      c = toupper(c);

    if((s_errf2 != "CMA") && (s_errf2 != "DEC"))
    {
      echec("Init égaliseur: la fonction d'erreur doit être \"DEC\" ou \"CMA\" (démandé : \"{}\").", s_errf);
      return;
    }

    if((s_structure2 != "FFE") && (s_structure2 != "DFE"))
    {
      echec("Init égaliseur: la struture doit être \"FFE\" ou \"DFE\" (démandé : \"{}\").", s_structure);
      return;
    }

    structure = s_structure2 == "FFE" ? FFE : DFE;
    errf      = s_errf2 == "DEC" ? DEC : CMA;

    // Identité par défaut (sortie = dernier échantillon reçu)
    h_avant    = vconcat(ArrayXf::Zero(N1-1), ArrayXf::Ones(1));
    wnd        = ArrayXf::Zero(N1);
    if(structure == DFE)
    {
      h_arriere  = ArrayXf::Zero(N2);
      wnd_deci   = ArrayXf::Zero(N2);
    }
    cnt        = 0;

    if((errf == CMA) && fo->infos.est_qam)
      msg_avert("Egaliseur RIF : la technique CMA n'est pas adaptée pour une modulation QAM.");

    init_ok = true;
  }

  void step(IArrayXcf x, ArrayXcf &y)
  {
    if(init_ok)
    {
      ArrayXf err;
      step_interne(x, y, err);
    }
  }
  void step_interne(IArrayXcf x, ArrayXcf &y, ArrayXf &err)
  {
    int n = x.rows();
    y.setZero(n);
    err.setZero(n);

    auto j = 0;
    for(auto i = 0; i < n; i++)
    {
      wnd.head(N1-1) = wnd.tail(N1-1).eval();
      wnd(N1-1)      = x(i);

      if(K > 1)
      {
        cnt = (cnt + 1) % K;
        if(cnt != 0)
          continue;
      }

      cfloat sortie      = (h_avant * wnd).sum();
      cfloat retroaction;

      if(structure == DFE)
      {
        retroaction = (h_arriere * wnd_deci).sum();
        sortie += retroaction;
      }


      cfloat décision;

      if((errf == DEC) || (structure == DFE))
        décision = fo->lis_symbole(fo->symbole_plus_proche(sortie));

      float e;
      if(errf == CMA)
      {
        e = 1 - carré(abs(sortie));
        h_avant   += α * e * (wnd * conj(sortie)).real();
        if(structure == DFE)
          h_arriere += α * e * (wnd_deci * conj(retroaction)).real();
      }
      else
      {
        auto ec = décision - sortie;
        e = abs(ec);
        h_avant   += α * (conj(ec) * wnd).real();
        if(structure == DFE)
          h_arriere += α * (conj(ec) * wnd_deci).real();
      }

      err(j) = e;
      y(j++) = sortie;

      if(structure == DFE)
      {
        wnd_deci.head(N2-1) = wnd_deci.tail(N2-1).eval();
        wnd_deci(N2-1)      = décision;
      }
    } // for i

    y.conservativeResize(j);
    err.conservativeResize(j);
  }

};


sptr<FiltreGen<cfloat>> égaliseur_rif_création(
    sptr<FormeOnde> forme_onde, const string &structure, const string &errf,
    int K, float α, int N1, int N2)
{
  return make_shared<EgaliseurRIF>(forme_onde, structure, errf, K, α, N1, N2);
}


Eigen::MatrixXf égaliseur_zfe_matrice(IArrayXf h, int n)
{
  int m = h.rows();
  // (2) Construit la matrice de filtrage : n + m - 1 lignes, n colonnes
  Eigen::MatrixXf F = Eigen::MatrixXf::Zero(n + m -1, n);

  // Ligne 0 à m-1
  for(auto i = 0; i < m; i++)
    F.block(i,0,1,i+1) = h.segment(0, i+1).reverse().transpose();

  // Ligne m à n-1
  for(auto i = m; i <= n-1; i++)
    F.block(i,i+1-m,1,m) = h.reverse().transpose();

  // Ligne n à n+m-2
  for(auto i = n; i < n+m-1; i++)
    for(auto j = 1; j <= m-1-(i-n); j++)
      F(i,n-j) = h(j+(i-n));

  return F;
}

ArrayXf égaliseur_zfe(IArrayXf h, int n)
{
  int m = h.rows();

  if(n < m)
  {
    msg_erreur("egaliseur_zfe(h[{}],n={}): n doit être supérieur ou égal à la longueur de la réponse du canal.", m, n);
    return ArrayXf();
  }

  Eigen::MatrixXf F = égaliseur_zfe_matrice(h, n);

  // (1) Calcule un délais
  // Suppose que h et g sont centrées, delai = m/2 + n/2
  auto d = (int) round((m + n)/2.0f);

  Eigen::VectorXf δ = Eigen::VectorXf::Zero(n + m - 1);
  δ(d) = 1.0f;

  // g = F \ δ;
  ArrayXf g = F.colPivHouseholderQr().solve(δ).array();

  //msg("F * g = {}", (F * g.matrix()).transpose());

  return g;
}






}


