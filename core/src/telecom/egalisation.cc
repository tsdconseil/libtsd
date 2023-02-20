#include "tsd/tsd.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/telecom.hpp"
#include "tsd/vue.hpp"
#include "tsd/eig-util.hpp"
#include <Eigen/QR>
#include <cctype>

using namespace tsd::vue;
using namespace std;

namespace tsd::telecom {


auto convol(const Vecf &h, const Vecf &x)
{
  soit f = tsd::filtrage::filtre_rif<float, float>(h);
  retourne f->step(x);
}

auto convol(const Vecf &h, const Veccf &x)
{
  soit f = tsd::filtrage::filtre_rif<float, cfloat>(h);
  retourne f->step(x);
}



struct EgaliseurRIF: FiltreGen<cfloat>
{
  sptr<FormeOnde> fo;
  entier K = 1;
  float α = 0.01;
  entier N1 = 11, N2 = 11, cnt = 0;

  Vecf h_avant, h_arriere;
  Veccf wnd, wnd_deci;

  bouléen init_ok = non;

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
      entier K, float α, entier N1, entier N2)
  {
    init_ok  = non;
    this->K  = K;
    this->α  = α;
    this->N1 = N1;
    this->N2 = N2;
    this->fo  = fo;

    tsd_assert_msg(fo, "Egaliseur RIF : la forme d'onde doit être spécifiée.");


    string s_errf2 = s_errf;
    pour(auto &c: s_errf2)
      c = toupper(c);
    string s_structure2 = s_structure;
    pour(auto &c: s_structure2)
      c = toupper(c);

    si((s_errf2 != "CMA") && (s_errf2 != "DEC"))
    {
      echec("Init égaliseur: la fonction d'erreur doit être \"DEC\" ou \"CMA\" (démandé : \"{}\").", s_errf);
      retourne;
    }

    si((s_structure2 != "FFE") && (s_structure2 != "DFE"))
    {
      echec("Init égaliseur: la struture doit être \"FFE\" ou \"DFE\" (démandé : \"{}\").", s_structure);
      retourne;
    }

    structure = s_structure2 == "FFE" ? FFE : DFE;
    errf      = s_errf2 == "DEC" ? DEC : CMA;

    // Identité par défaut (sortie = dernier échantillon reçu)
    h_avant    = vconcat(Vecf::zeros(N1-1), Vecf::ones(1));
    wnd        = Vecf::zeros(N1);
    si(structure == DFE)
    {
      h_arriere  = Vecf::zeros(N2);
      wnd_deci   = Vecf::zeros(N2);
    }
    cnt        = 0;

    si((errf == CMA) && fo->infos.est_qam)
      msg_avert("Egaliseur RIF : la technique CMA n'est pas adaptée pour une modulation QAM.");

    init_ok = oui;
  }

  void step(const Veccf &x, Veccf &y)
  {
    si(init_ok)
    {
      Vecf err;
      step_interne(x, y, err);
    }
  }
  void step_interne(const Veccf &x, Veccf &y, Vecf &err)
  {
    soit n = x.rows(), j = 0;
    y.setZero(n);
    err.setZero(n);

    pour(auto i = 0; i < n; i++)
    {
      wnd.head(N1-1) = wnd.tail(N1-1);
      wnd(N1-1)      = x(i);

      si(K > 1)
      {
        cnt = (cnt + 1) % K;
        si(cnt != 0)
          continue;
      }

      soit sortie      = (h_avant * wnd).somme();
      cfloat retroaction;

      si(structure == DFE)
      {
        retroaction = (h_arriere * wnd_deci).somme();
        sortie += retroaction;
      }


      cfloat décision;

      si((errf == DEC) || (structure == DFE))
        décision = fo->lis_symbole(fo->symbole_plus_proche(sortie));

      float e;
      si(errf == CMA)
      {
        e = 1 - carré(abs(sortie));
        h_avant   += α * e * real(wnd * conj(sortie));
        si(structure == DFE)
          h_arriere += α * e * real(wnd_deci * conj(retroaction));
      }
      sinon
      {
        soit ec = décision - sortie;
        e = abs(ec);
        h_avant   += α * real(wnd * conj(ec));
        si(structure == DFE)
          h_arriere += α * real(wnd_deci * conj(ec));
      }

      err(j) = e;
      y(j++) = sortie;

      si(structure == DFE)
      {
        wnd_deci.head(N2-1) = wnd_deci.tail(N2-1).eval();
        wnd_deci(N2-1)      = décision;
      }
    } // pour i

    y.conservativeResize(j);
    err.conservativeResize(j);
  }

};


sptr<FiltreGen<cfloat>> égaliseur_rif_création(
    sptr<FormeOnde> forme_onde, const string &structure, const string &errf,
    entier K, float α, entier N1, entier N2)
{
  retourne make_shared<EgaliseurRIF>(forme_onde, structure, errf, K, α, N1, N2);
}


Tabf égaliseur_zfe_matrice(const Vecf &h, entier n)
{
  soit m = h.rows();
  // (2) Construit la matrice de filtrage : n + m - 1 lignes, n colonnes
  soit F = Tabf::zeros(n + m -1, n);

  // Ligne 0 à m-1
  pour(auto i = 0; i < m; i++)
    //F.block(i,0,1,i+1) = h.segment(0, i+1).reverse().transpose();
    //F.set_row(i, h.head(i+1).reverse());
    for(auto k = 0; k < i+1; k++)
      F(i, k) = h(i-k);

  // Ligne m à n-1
  pour(auto i = m; i <= n-1; i++)
    for(auto k = 0; k < m; k++)
      F(i, i+1-m+k) = h(m-1-k);
    //F.set_row(i, h.reverse());
    //F.block(i,i+1-m,1,m) = h.reverse().transpose();

  // Ligne n à n+m-2
  pour(auto i = n; i < n+m-1; i++)
    pour(auto j = 1; j <= m-1-(i-n); j++)
      F(i,n-j) = h(j+(i-n));

  retourne F;
}

Vecf égaliseur_zfe(const Vecf &h, entier n)
{
  soit m = h.rows();

  si(n < m)
  {
    msg_erreur("egaliseur_zfe(h[{}],n={}): n doit être supérieur ou égal à la longueur de la réponse du canal.", m, n);
    retourne {};
  }

  soit F = égaliseur_zfe_matrice(h, n);

  // (1) Calcule un délais
  // Suppose que h et g sont centrées, delai = m/2 + n/2
  soit d = (entier) round((m + n)/2.0f);

  Eigen::VectorXf δ = Eigen::VectorXf::Zero(n + m - 1);
  δ(d) = 1.0f;

  // g = F \ δ;
  Eigen::VectorXf x = tab2etab(F).colPivHouseholderQr().solve(δ);
  retourne evec2vec(x);

  //msg("F * g = {}", (F * g.matrix()).transpose());
}






}


