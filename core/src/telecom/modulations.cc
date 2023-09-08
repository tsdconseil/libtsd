#include "tsd/tsd.hpp"
#include "tsd/telecom.hpp"
#include "tsd/telecom/bitstream.hpp"
#include "tsd/filtrage.hpp"
#include "tsd/vue.hpp"
#include <iostream>
#include <cassert>

using namespace tsd::vue;
using namespace tsd::filtrage;
using namespace std::complex_literals;

namespace tsd::telecom {


/*template<typename T>
  void clean(T &t)
{
    t = (t.abs() > 1e-8f).select(t, 0.0f);
}*/


static inline entier decode_QPSK(const cfloat &x)
{
  si(x.real() >= 0)
  {
    si(x.imag() >= 0)
      retourne 0;
    sinon
      retourne 3;
  }
  sinon si(x.imag() >= 0)
    retourne 1;
  sinon
    retourne 2;
}

static inline entier decode_BPSK(const cfloat &x)
{
  //retourne (x.real() >= 0) ? 0 : 1;
  retourne (x.real() >= 0) ? 1 : 0;
}

static Veccf psk_constellation(entier M)
{
  si(M == 2)
    retourne Veccf::valeurs({-1, 1});

  float dec = (M == 4) ? π_f / 4 : 0;

  retourne Veccf::int_expr(M, IMAP(
      std::polar(1.0f, dec + (2 * π_f * i) / M)));
}

static Vecf ask_constellation(entier M, float K1, float K2)
{
  retourne K1 + linspace(0, M-1, M) * (K2 / (M-1));
}

template<typename T>
  void clean(Vecteur<std::complex<T>> &t)
{
  soit tr = (t.real().abs() > 1e-7f).select(t.real(), 0.0f);
  soit ti = (t.imag().abs() > 1e-7f).select(t.imag(), 0.0f);
  t = tr + 1.0if * ti;
}


template<typename T>
 auto modulo(const Vecteur<T> &x, unsigned int M)
{
  retourne x - M * (x * (1.0f / M)).floor();
}



// Train binaire -> index
Veci symmap_binaire(const BitStream &x, entier k)
{
  entier n = x.lon();

  si(k <= 0)
    échec("symmap_binaire : k <= 0.");

  entier nsymbs = (n + k-1)/k;
  Veci y(nsymbs);
  pour(auto i = 0; i < nsymbs; i++)
  {
    entier symb = 0;
    entier e = 1;
    pour(auto j = 0; (j < k) && (i*k+j < n); j++)
    {
      symb = symb + e * x[i*k+j];
      e   *= 2;
    }
    y(i) = symb;
  }
  retourne y;
}

// Index -> train binaire
void symdemap_binaire(BitStream &bs, const Veci &x, entier k)
{
  entier n = x.rows();
  pour(auto i = 0; i < n; i++)
    pour(auto j = 0; j < k; j++)
      bs.push((((entier) x(i)) >> j) & 1);
}

Veccf FormeOnde::génère_symboles(const BitStream &bs)
{
  soit x = symmap_binaire(bs, infos.k);
  soit n = x.rows();
  Veccf y(n);
  pour(auto i = 0; i < n; i++)
  {
    cnt = i;
    y(i) = lis_symbole(x(i));
  }
  retourne y;
}

struct CtxGenImpl: FormeOnde::CtxGen
{
  FormeOnde *parent;
  sptr<FiltreGen<cfloat>> f;
  entier ncoefs, osf;

  CtxGenImpl(FormeOnde *parent, entier ncoefs, entier osf)
  {
    this->parent = parent;
    f = parent->filtre.filtre_mise_en_forme(ncoefs, osf);
    this->osf     = osf;
    this->ncoefs  = parent->filtre.get_coefs(ncoefs, osf).rows();
  }
  void reset()
  {

  }
  Veccf step(const BitStream &bs)
  {
    retourne {};
    /*soit &mconfig = config.format.modulation;

    soit osf = mconfig.fe / mconfig.fsymb;

    entier   d_ech           = (entier) ceil(mod->delais());
    entier   nbits_par_symb  = mconfig.wf->k;
    float nbits_par_echan = ((float) nbits_par_symb) / osf;
    entier   d_bit           = (entier) ceil(d_ech * nbits_par_echan);

    // Ici, non, il faut concaténer deux modulateurs
    retourne mod->step(bs);*/
  }
};


sptr<FormeOnde::CtxGen> FormeOnde::get_contexte_tx(entier ncoefs, entier osf)
{
  retourne std::make_shared<CtxGenImpl>(this, ncoefs, osf);
}


// TODO : stocker dans le (un) contexte les différents objet de filtrage
Veccf FormeOnde::génère_échantillons(const BitStream &bs, entier ncoefs,
    entier osf, float &retard)
{
  soit symbs = génère_symboles(bs);
  soit f = filtre.filtre_mise_en_forme(ncoefs, osf);
  ncoefs = filtre.get_coefs(ncoefs, osf).rows();
  entier nsymbs_flush = (ncoefs + osf - 1) / osf;//ncoefs;
  symbs = vconcat(symbs, Veccf::zeros(nsymbs_flush));

  // Temps vers le milieu du premier bit transmis
  retard = filtre_rif_ups_délais(ncoefs, osf);

  soit y = f->step(symbs);

  si(infos.est_fsk)
  {
    // df = 0.5 * h * fsymb
    // df sur un symbole = 0.5 * h / osf
    // 2 π / osf <=> h = 2
    // => θ = h * 2 * pi / osf / 2 = h * pi / osf

    float Ω_max = (π_f * infos.index) / osf;
    // h = 2 -> Omega_max = 2 * pi / osf

    msg("Ω max = {} degrés.", rad2deg(Ω_max));

    // h = 2 fd / fsymb = excursion / fsymb
    // => omega_max = 2%pi*fd = %pi * h * fsymb
    // => Omega_max = omega_max / fs
    //              = %pi * h / ovs;
    // Conversion phase -> IQ


    // normalisation entre [-fmax,fmax]
    soit vfreqs = real(y) * (Ω_max / abs(real(y)).valeur_max());

    // TODO : échantillon précédent !
    y = polar(cumsum(vfreqs));
  }

  retourne y;
}


struct CtxSansMemoire: FormeOnde::Ctx
{
  const FormeOnde *fo;
  CtxSansMemoire(const FormeOnde *fo)
  {
    this->fo = fo;
  }
  // Pas de contexte
  void reset(){}


  tuple<entier, cfloat> step(cfloat x)
  {
    entier i = fo->symbole_plus_proche(x);
    retourne {i, fo->lis_symbole(i)};
  }

  /*void avance(entier symbole){}
  cfloat lis_symbole(entier i)
  {
    retourne fo->lis_symbole(i);
  }
  entier symbole_plus_proche(cfloat x)
  {
    retourne fo->symbole_plus_proche(x);
  }*/
};

sptr<FormeOnde::Ctx> FormeOnde::get_ctx(entier OSF) const
{
  retourne std::make_shared<CtxSansMemoire>(this);
}

// Optimisable
void  FormeOnde::decode_symboles(BitStream &bs, const Veccf &x)
{
  soit n = x.rows();
  Veci x2(n);
  pour(auto i = 0; i < n; i++)
  {
    cnt = i;
    x2(i) = symbole_plus_proche(x(i));
  }
  symdemap_binaire(bs, x2, infos.k);
}


float FormeOnde::excursion() const
{
  retourne 1.0f;
}

entier FormeOnde::symbole_plus_proche(const cfloat &point) const
{
  soit bdist = 1e100;
  soit res = 0;
  pour(auto i = 0; i < infos.M; i++)
  {
    soit d2 = std::norm(point - lis_symbole(i));
    si(d2 < bdist)
    {
      bdist = d2;
      res = i;
    }
  }
  assertion((res >= 0) && (res < (entier) infos.M));
  retourne res;
}


struct FormeOndeASK: FormeOnde
{
  Veccf symbs;
  float K1 = -1, K2 = 2;

  entier symbole_plus_proche(const cfloat &x) const
  {
    retourne FormeOnde::symbole_plus_proche(x);
  }

  string desc_courte() const
  {
    retourne sformat("{}-ASK({},{})", infos.M, K1, K2);
  }

  string desc() const
  {
    retourne sformat("{}-ASK({},{})", infos.M, K1, K2);
  }


  FormeOndeASK(entier M, float K1, float K2, const SpecFiltreMiseEnForme &filtre)
  {
    this->filtre        = filtre;
    infos.M             = M;
    this->K1            = K1;
    this->K2            = K2;
    infos.est_lineaire  = oui;
    infos.est_ask       = oui;
    infos.k             = log2(M);
    symbs               = ask_constellation(M, K1, K2);
  }

  Veccf constellation() const
  {
    retourne symbs;
  }

  float ber(float EbN0)
  {
    // Equation 5.2.46 de Proakis
    float EbN0_lin = db2pow(EbN0);
    retourne ((infos.M-1.0f) / infos.M)
        * erfc(sqrt(3.0f * log2((float) infos.M) * EbN0_lin / (infos.M*infos.M - 1.0f)))/infos.k;
  }

  cfloat lis_symbole(unsigned int index) const
  {
    assertion((entier) index < infos.M);
    retourne symbs(index);
  }
};


struct FormeOndePSK: FormeOnde
{
  Veccf symbs;

  entier symbole_plus_proche(const cfloat &x) const
  {
    si(infos.M == 2)
      retourne x.real() > 0 ? 1 : 0;
    retourne FormeOnde::symbole_plus_proche(x);
  }

  string desc_courte() const
  {
    si(infos.M == 2)
      retourne "BPSK";
    sinon si(infos.M == 4)
      retourne "QPSK";
    sinon
      retourne sformat("{}PSK", infos.M);
  }

  string desc() const
  {
    retourne desc_courte() + sformat(", {}", filtre);
  }


  FormeOndePSK(unsigned int M, const SpecFiltreMiseEnForme &filtre)
  {
    this->filtre = filtre;
    infos.est_lineaire = oui;
    infos.M = M;
    infos.est_psk = oui;
    infos.k = std::log2(M);
    symbs = psk_constellation(M);
  }

  Veccf constellation() const
  {
    retourne symbs;
  }

  float ber(float EbN0)
  {
    float EbN0_lin = db2pow(EbN0);
    //eq 5.2.61 from proakis


    // Equation générale (valable sauf en BPSK)
    float ber = erfc(sqrt(infos.k*EbN0_lin)*sin(π/infos.M))/infos.k;

    // si k = 1, M = 2 (BPSK) :
    // sin(pi/2) = 1
    // -> erfc(std::sqrt(ebnoD))

    // si k = 2, M = 4 (QPSK) :
    // sin(pi/4) = sqrt(2)/2
    // erfc(std::sqrt(ebnoD) * sqrt(2) * sqrt(2) / 2) = erfc(sqrt(ebnoD)) / 2

    si(infos.M == 2)
      ber /= 2;
    retourne ber;
    //si(opt=='i')
      //ber=2*ber; // Differential modulation
    //end
  }

  cfloat lis_symbole(unsigned int index) const
  {
    assertion((entier) index < infos.M);
    retourne symbs(index);
  }
};


struct FormeOnde_π4QPSK: FormeOnde
{
  Veccf symbs[2];

  entier symbole_plus_proche(const cfloat &x) const
  {
    retourne FormeOnde::symbole_plus_proche(x);
  }

  string desc_courte() const
  {
    retourne "pi4-QPSK";
  }

  string desc() const
  {
    retourne sformat("pi4-QPSK, {}", filtre);
  }

  struct Ctxπ4QPSK: FormeOnde::Ctx
  {
    const FormeOnde_π4QPSK *parent;
    entier cnt = 0;
    const cfloat rot = std::polar(1.0f, -π_f/4);

    Ctxπ4QPSK(const FormeOnde_π4QPSK *parent)
    {
      this->parent = parent;
    }

    void reset()
    {
      cnt = 0;
    }

    /** Index = -1 si pas d'échantillon à sortir */
    tuple<entier, cfloat> step(cfloat x)
    {
      si(cnt++ & 1)
        x *= rot;
      entier r = decode_QPSK(x);
      retourne {r, parent->symbs[(cnt+1) & 1](r)};
    }
  };

  sptr<Ctx> get_ctx(entier OSF) const
  {
    retourne std::make_shared<Ctxπ4QPSK>(this);
  }


  FormeOnde_π4QPSK(const SpecFiltreMiseEnForme &filtre)
  {
    this->filtre = filtre;
    infos.est_lineaire  = oui;
    infos.M             = 4;
    infos.est_psk       = oui;
    infos.k             = 2;
    symbs[0] = psk_constellation(4);
    symbs[1] = symbs[0] * exp(ⅈ * π_f / 4.0f);
  }

  Veccf constellation() const
  {
    retourne vconcat(symbs[0], symbs[1]);
  }

  float ber(float EbN0)
  {
    // Ber idem QPSK
    float EbN0_lin = db2pow(EbN0);
    retourne erfc(sqrt(infos.k*EbN0_lin)*sin(π/infos.M))/infos.k;
  }

  cfloat lis_symbole(unsigned int index) const
  {
    assertion((entier) index < infos.M);
    retourne symbs[cnt & 1](index);
  }
};


struct FormeOndeQAM: FormeOnde
{
  Veccf symbs;

  string desc_courte() const
  {
    retourne sformat("QAM{}", infos.M);
  }

  string desc() const
  {
    retourne sformat("QAM{}, {}", infos.M, filtre);
  }

  FormeOndeQAM(unsigned int M, const SpecFiltreMiseEnForme &filtre)
  {
    this->filtre = filtre;
    infos.est_lineaire  = oui;
    infos.M             = M;
    infos.est_qam       = oui;
    infos.k             = log2(M);
    //nom = sformat("QAM{}", M);

    soit M2 = (entier) sqrt(M);

    si(M2 * M2 != M)
    {
      msg_erreur("FormeOnde QAM : M devrait être un carré (M = {}).", M);
      retourne;
    }

    symbs.resize(M);

    Tabcf map(M2, M2);
    pour(auto x = 0; x < M2; x++)
    {
      pour(auto y = 0; y < M2; y++)
      {
        soit xf = (float) x, yf = (float) y;

        xf /= (M2 - 1);
        yf /= (M2 - 1);

        // xf, yf compris entre 0 1
        // 0 -> -1
        // 1 ->  1
        map(x, y) = cfloat(xf * 2 - 1, yf * 2 - 1);
      }
    }
    symbs = map.reshape(M2 * M2);
  }

  Veccf constellation() const
  {
    retourne symbs;
  }

  /*ArrayXcf gene_symboles(const ArrayXf &bs)
  {
    retourne symmap(bs, k, 'q');
  }
  ArrayXf demodulation(const ArrayXcf &x)
  {
    retourne symdemap(x, k, 'q');
  }*/

  cfloat lis_symbole(unsigned int index) const
  {
    assertion((entier) index < infos.M);
    retourne symbs(index);
  }
  float ber(float EbN0)
  {
    float ebnoD = pow(10.0f, EbN0/10);
    // D'après proakis équation 5.2.79
    retourne (2.0f/infos.k)*(1-(1/sqrt(infos.M)))*erfc(sqrt(3.0f*infos.k*ebnoD/(2.0f*(infos.M-1))));
  }
};






struct FormeOndeFSK: FormeOnde
{
  Vecf symbs;

  struct CtxFSK: FormeOnde::Ctx
  {
    entier OSF; // Comment le configurer ?
    const FormeOndeFSK *fo;

    // Etat de phase initial du symbole en cours
    cfloat etat_phase = 0;
    // Nombre de pas depuis le début du symbole en cours
    entier    cnt = 0;

    Veccf symbs;

    // Tableau d'erreur pour chaque trajectoire de phase
    Vecf   err;
    Veccf  cands;

    CtxFSK(const FormeOndeFSK *fo, entier OSF)
    {
      this->OSF = OSF;
      this->fo = fo;

      entier M = fo->infos.M;

      symbs.resize(M);
      si(M != 2)
        échec("TODO: CtxFSK : M != 2");

      // si index = 2, il faut un OSF > 1
      symbs(0) = exp(- 2 * π_f * ⅈ * fo->infos.index / (2.0f * OSF));
      symbs(1) = exp(+ 2 * π_f * ⅈ * fo->infos.index / (2.0f * OSF));

      err.resize(M);
      cands.resize(M);
      reset();
    }

    void reset()
    {
      etat_phase = 0;
      cnt        = 0;
      err.setZero();
      cands.setConstant(1.0f);
    }

    tuple<entier, cfloat> step(cfloat x)
    {
      // pour chaque symbole, on fait M hypothèses
      pour(auto i = 0; i < (entier) fo->infos.M; i++)
      {
        cands(i) = (cands(i) * symbs(i)) / abs(cands(i));
        err(i)   += std::norm(cands(i) - x); // std::norm ABS2
      }
      cnt++;
      si(cnt == OSF)
      {
        entier i;
        cnt = 0;
        err.minCoeff(&i);
        err.setZero();
        cands.setConstant(cands(i));
        retourne {i, cands(i)};
      }
      sinon
      {
        retourne {-1, 0}; // Pas de symbole à sortir pour l'instant
      }
    }

    /*entier step(cfloat obs)
    {

      retourne 0;
    }*/



    /*void avance(entier i)
    {
      etat *= lis_symbole(i);
    }
    cfloat lis_symbole(entier i) const
    {
      retourne etat_phase * pow(symbs(i), cnt);
      // retourne fo->lis_symbole(i);
    }*/
  };


  /*sptr<Ctx> get_ctx(entier OSF) const
  {
    retourne std::make_shared<CtxFSK>(this, OSF);
  }*/

  string desc_courte() const
  {
    string nom = format("{}{}SK",
        filtre.type == SpecFiltreMiseEnForme::Type::GAUSSIEN ? "G" : "",
            infos.index == 0.5 ? "M" : "F");
    si(infos.M != 2)
      nom = sformat("{}{}", infos.M, nom);
    retourne nom;
  }

  string desc() const
  {
    retourne desc_courte() + sformat(", index={}", infos.index);
  }

  float excursion() const
  {
    retourne infos.index;
  }

  FormeOndeFSK(unsigned int M, float index, const SpecFiltreMiseEnForme &filtre)
  {
    infos.est_lineaire  = non;
    infos.M             = M;
    infos.index         = index;
    this->filtre        = filtre;
    infos.est_fsk       = oui;
    infos.k             = log2(M);

    symbs.resize(M);

    pour(auto x = 0u; x < M; x++)
    {
      float xf = (float) x;

      xf /= (M - 1);

      // xf, yf compris entre 0 1
      // 0 -> -1
      // 1 ->  1
      symbs(x) = xf * 2 - 1;
    }
  }

  Veccf constellation() const
  {
    retourne {};
  }
  float ber(float EbN0)
  {
    // Détecteur FSK cohérent et orthogonal (?)
    /*float F = pow(2.0f,k-1.0f)/(std::pow(2.0f,(float)k)-1);
    function [y]=fskcoherorth(q)
      A=sqrt(2*k*ebnoD(i)),
      B=q*q/2,
      C=0.5*erfc((-q-A)/sqrt(2)),
      y=(C^(M-1))/sqrt(2*(%pi))*exp(-B),
    endfunction
    Ps(i)=1-intg(-100,100,fskcoherorth);
    Pb(i)=F*Ps(i);*/
    // Détecteur FSK cohérent et non orthogonal
    float ebnoD = pow(10.0f, EbN0/10);
    cfloat Rho = 0.0f; // ?
    float A = ebnoD*(1-real(Rho));
    retourne 0.5 * erfc(sqrt(A/2));
  }
  // En supposant une discrimination polaire faite avant...
  cfloat lis_symbole(unsigned int index) const
  {
    //retourne 0.0f;
    assertion((entier) index < infos.M);
    retourne symbs(index);
  }
  float etat_phase = 0;
  cfloat lis_symbole(entier index) const
  {
    //retourne 0.0f;
    assertion(index < infos.M);
    retourne symbs(index);
  }
};


sptr<FormeOnde> forme_onde_fsk(unsigned int M, float index, const SpecFiltreMiseEnForme &filtre)
{
  retourne make_shared<FormeOndeFSK>(M, index, filtre);
}


sptr<FormeOnde> forme_onde_psk(unsigned int M, const SpecFiltreMiseEnForme &filtre)
{
  retourne make_shared<FormeOndePSK>(M, filtre);
}

sptr<FormeOnde> forme_onde_ask(entier M, float M1, float M2, const SpecFiltreMiseEnForme &filtre)
{
  retourne make_shared<FormeOndeASK>(M, M1, M2, filtre);
}

sptr<FormeOnde> forme_onde_qam(unsigned int M, const SpecFiltreMiseEnForme &filtre)
{
  retourne make_shared<FormeOndeQAM>(M, filtre);
}

sptr<FormeOnde> forme_onde_bpsk(const SpecFiltreMiseEnForme &filtre)
{
  retourne forme_onde_psk(2, filtre);
}

sptr<FormeOnde> forme_onde_qpsk(const SpecFiltreMiseEnForme &filtre)
{
  retourne forme_onde_psk(4, filtre);
}





sptr<FormeOnde> forme_onde_π4_qpsk(const SpecFiltreMiseEnForme &filtre)
{
  retourne std::make_shared<FormeOnde_π4QPSK>(filtre);
}




Vecf SpecFiltreMiseEnForme::get_coefs(entier ncoefs, entier osf) const
{
  assertion_msg(osf > 0, "SpecFiltreMiseEnForme : osf invalide ({})", osf);

  soit coefs = Vecf::zeros(osf);
  coefs(0) = 1;

  si(osf == 1)
  {
    retourne coefs; // Filtre neutre
    //échec("Spec filtre mise en forme::get_coefs() : osf = 1.");
  }

  si(ncoefs == 0)
  {
    ncoefs = 5 * osf + 1;
    si((ncoefs & 1) == 0)
      ncoefs++;
  }

  // Construction du filtre
  si(type == NRZ)
  {
    retourne Vecf::ones(osf) / osf; // Ntaps = osf = 4
  }
  sinon si(type == AUCUN)
  {
    retourne Vecf::ones(1);
  }
  sinon si(type == GAUSSIEN)
  {
    assertion(ncoefs > 0);
    //msg("gaussien : bt = {}, osf = {}, ncoefs = {}", BT, osf, ntaps);
    // Porte de largeur osf + filtre gaussien
    coefs = design_rif_gaussien_telecom(ncoefs, BT, osf);
  }
  sinon si(type == RCS)
  {
    // Ntaps = 3 * osf = 12
    assertion(ncoefs > 0);
    //msg("rcs : roff = {}, osf = {}, ncoefs = {}", α, osf, ncoefs);
    coefs = design_rif_rcs1(ncoefs, β, osf);
  }
  sinon
    échec("Type de filtre de mise en forme inconnu ({})", (entier) type);

  /*{
    Figure fg("Filtre de mise en forme");
    fg.plot(coefs, "b|r");
    stdo << fg;
  }*/

  //msg("Filtre de mise en forme : ncoefs = {} ({} prog, osf = {}).", coefs.rows(), ncoefs, osf);

  assertion(coefs.rows() > 0);
  si(coefs.hasNaN())
    échec("{} : Valeurs NaN dans les coefficients.", *this);

  retourne coefs;
}

sptr<FiltreGen<cfloat>> SpecFiltreMiseEnForme::filtre_adapté(entier ntaps, entier osf) const
{
  si(osf == 1)
    retourne filtre_id<cfloat>();
  soit h = get_coefs(ntaps, osf);
  retourne filtre_rif<float, cfloat>(h);
}

sptr<FiltreGen<cfloat>> SpecFiltreMiseEnForme::filtre_adapté_décimation(entier ntaps, entier osf) const
{
  si(osf == 1)
    retourne filtre_id<cfloat>();
  soit h = get_coefs(ntaps, osf);
  retourne filtre_rif_decim<float, cfloat>(h, osf);
}

sptr<FiltreGen<cfloat>> SpecFiltreMiseEnForme::filtre_mise_en_forme(entier ntaps, entier R) const
{
  si(R == 1)
    retourne filtre_id<cfloat>();

  soit h = get_coefs(ntaps, R);
  // Normalisation en énergie
  soit en = square(h).somme();

  // msg("filtre_mise_en_forme: en={}, R={}", en, R);

  // Facteur 1/R à cause de normalisation par R faite dans filtre_rif_ups (d'ailleurs à vérifier)
  h *= (sqrt(R) / sqrt(en)) / R;

  retourne filtre_rif_ups<float, cfloat>(h, R);
}



Vecf FormeOnde::ber(const Vecf &EbN0_dB)
{
  retourne Vecf::int_expr(EbN0_dB.rows(),
      IMAP(ber(EbN0_dB(i))));
}

std::ostream& operator<<(std::ostream &ss, const SpecFiltreMiseEnForme &t)
{
  //ss << "SpecFiltreMiseEnForme : ";
  si(t.type == SpecFiltreMiseEnForme::NRZ)
    ss << "NRZ";
  sinon si(t.type == SpecFiltreMiseEnForme::AUCUN)
    ss << "AUCUN";
  sinon si(t.type == SpecFiltreMiseEnForme::RCS)
    ss << sformat("RCS-dep={}", (entier) (t.β * 100.0f));
  sinon si(t.type == SpecFiltreMiseEnForme::GAUSSIEN)
    ss << sformat("Gaussien-BT={}", t.BT);
  sinon
    ss << "?";
  retourne ss;
}

std::ostream& operator<<(std::ostream &ss, const FormeOnde &t)
{
  ss << t.desc();
  retourne ss;
}

SpecFiltreMiseEnForme SpecFiltreMiseEnForme::aucun()
{
  retourne SpecFiltreMiseEnForme{.type = Type::AUCUN};
}

SpecFiltreMiseEnForme SpecFiltreMiseEnForme::gaussien(float BT)
{
  retourne SpecFiltreMiseEnForme{.type = Type::GAUSSIEN, .BT = BT};
}
SpecFiltreMiseEnForme SpecFiltreMiseEnForme::nrz()
{
  retourne SpecFiltreMiseEnForme{.type = Type::NRZ};
}
SpecFiltreMiseEnForme SpecFiltreMiseEnForme::rcs(float β)
{
  retourne SpecFiltreMiseEnForme{.type = Type::RCS, .β = β};
}

// Ajoute si besoin des zéros avant et après
static Vecf sympad(const Vecf &v, entier nc)
{
  assertion(v.rows() <= nc);
  soit r = v.clone();
  si(v.rows() < nc)
  {
    soit d = (nc - v.rows()) / 2;
    r = Vecf::zeros(d) | v | Vecf::zeros(nc - d - v.rows());
  }
  assertion(r.rows() == nc);
  retourne r;
}

SpecFiltreMiseEnForme::Analyse SpecFiltreMiseEnForme::analyse(entier ncoefs, entier osf) const
{
  SpecFiltreMiseEnForme::Analyse res;
  Figures f;

  // h, ncoefs : réponse réelle, avec ncoefs
  // ht, nc2   : approximation de la réponse théorique, avec beaucoup de coefficients

  soit h = sympad(get_coefs(ncoefs, osf), ncoefs);

  soit nc2  = 8*1024-1,
       osf2 = max(20, 2 * osf);
  soit ht = sympad(get_coefs(nc2, osf2), nc2);

  // Supprime les échantillons nuls, pour clarifier l'affichage
  soit i = 0;
  tantque((i < nc2/2) && (abs(ht(i)) < 1e-4) && ((nc2 - 2 * i) / ((1.0 * osf2) / osf) > ncoefs))// && (<abs(ht(i)) < h(0))<)
    i++;
  nc2 -= 2 * i;
  ht = sympad(get_coefs(nc2, osf2), nc2);


  soit m = ncoefs/2, m2 = nc2/2;
  // Intervalles temporels, en nombres de symboles
  soit t  = linspace(-(1.0f*m)/osf, (1.0f*m)/osf, ncoefs);
  soit t2 = linspace(-(1.0f*m2)/osf2, (1.0f*m2)/osf2, nc2);

  // Ratio entre les fréquences d'échantillonnage
  soit r = ((float) osf2) / osf;

  f.subplot();
  f.gcf().plot(t,  h,   "|bo", "Réel");
  f.gcf().plot(t2, ht * r, "-g",  "Théorique"); // ht * r : pour avoir la même amplitude
  f.gcf().titres(sformat("{} - {} coefs, OSF={}", *this, ncoefs, osf), "Symboles");
  f.subplot();
  soit [fr, mag]   = frmag(h);
  soit [fr2, mag2] = frmag(ht);
  f.gcf().plot(fr, mag,  "-b", "Réel");

  soit ns = (entier) (mag2.rows() / r);
  Vecf fr3(ns);
  pour(auto i = 0; i < ns; i++)
    fr3(i) = fr2((entier)floor(i * r));
  f.gcf().plot(fr3, mag2.head(ns), "-g", "Théorique");
  f.gcf().titre("Réponse fréquentielle");


  f.subplot();
  f.gcf().plot(fr,  mag2db(mag),           "-b", "Réel");
  f.gcf().plot(fr3, mag2db(mag2.head(ns)), "-g", "Théorique");
  f.gcf().titre("Réponse fréquentielle (log)");


  res.fig = f;
  retourne res;
}

}

