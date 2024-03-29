#include "tsd/telecom/lfsr.hpp"

namespace tsd::telecom {

/** @brief Maximum number of bits of the PRBS polynomial */
static const auto MAX_REGLEN = 16u;


/** From Xilinx xapp-052 (July 7, 1996), page 5 */
// Les polynômes sont stockés ainsi (n étant le degré) :
// Bit 0 : b0 * x^n (implicite : b0 = 1)
// Bit 1 : b1 * x^(n-1)
// Bit n : implicite (x^0) et en pratique non utilisé
//
// Structure du registre à décalage (avec les MSB en premier) :
// [r_n-1 r_n-2  ...  r_0]
// [b1    b1     ...  b_n]



// Exemples de polynômes primitif, encodés en binaire
// d'après "On Pseudo-Random and Orthogonal Binary Spreading Sequences", Mitra, 2007
// (bit de poids faible implicite = x^n,
//  bit de poids fort implicite = x^0)
const vector<uint32_t> pols_prim =
   {0,                                // n=0, n/a
    0,                                // n=1, 1 + x
    (1 << 1),                         // n=2, 1 + x + x^2
    (1 << 2),                         // n=3, 1 + x + x^3
    (1 << 3),                         // n=4, 1 + x + x^4
    (1 << 3),                         // n=5, 1 + x^2 + x^5
    (1 << 5),                         // n=6, 1 + x + x^6
    (1 << 6),                         // n=7, 1 + x + x^7
    (1 << 6) | (1 << 5) | (1 << 4),   // n=8, 1 + x^2 + x^3 + x^4 + x^8
    (1 << 5),                         // n=9
    (1 << 7),                         // n=10
    (1 << 9),                         // n=11
    (1 << 9) | (1 << 8) | (1 << 5),   // n=12, 1+x^3+x^4+x^7+x^12
    (1 << 9) | (1 << 10) | (1 << 12), // n=13, 1+x+x^3+x^4+x^13
    (1 << 5) | (1 << 3) | (1 << 1),
    (1 << 14),
    (1 << 5) | (1 << 3) | (1 << 2)
};

uint32_t polynome_primitif_binaire(entier n)
{
  si((n <= 0) || (n >= (entier) pols_prim.size()))
    échec("polynome_primitif : degré = {}, pas possible.", n);
  retourne pols_prim[n] | 1; // Ajoute X^n, implicite
}

Poly<entier> polynome_primitif(entier n)
{
  Poly<entier> res;
  soit z = Poly<entier>::z;
  soit p = polynome_primitif_binaire(n);
  pour(auto i = 0; i < n; i++)
  {
    si(p & (1 << i))
      res = res + (z ^ (n-i));
  }
  res = res + Poly<entier>::one(); // Ajoute c0 = 1, non utilisé dans les calculs
  retourne res;
}

BitStream code_mls(entier reglen)
{
  soit m = pow(2,reglen) - 1;
  soit pb = polynome_primitif_binaire(reglen);
  LFSRGenerateur gene;
  LFSRConfig config;
  config.sortie = LFSRConfig::POIDS_FAIBLE;
  config.reglen = reglen;
  config.pol    = pb;
  config.p0     = 1;
  gene.configure(config);
  BitStream bs;
  gene.step(bs, m);
  retourne bs;
}

LFSRGenerateur::LFSRGenerateur(unsigned int reglen)
{
  configure(reglen);
}

void LFSRGenerateur::configure(unsigned int reglen)
{
  assertion_msg(reglen <= MAX_REGLEN, "Invalid PRBS reg len: {}.", reglen);

  config.reglen = reglen;
  config.pol    = polynome_primitif_binaire(reglen);//best_pols[reglen] | 1;
  //pol    = (best_pols[reglen]<<1) | 1;
  reg    = 0xfa17;//0x0001;
}

void LFSRGenerateur::configure(const LFSRConfig &config)
{
  this->config = config;
  reg = config.p0;
}


void LFSRGenerateur::step(BitStream &bs, unsigned int nbits)
{
  pour(auto i = 0u; i < nbits; i++)
  {
    entier somme = (__builtin_popcount(reg & config.pol) & 1);
    entier sortie;

    si(config.sortie == LFSRConfig::POL)
    {
      soit tmp = reg & config.pol_sortie;
      soit s = 0u;
      pour(auto k = 0; k < config.reglen; k++)
        s = s ^ ((tmp >> k) & 1);
      sortie = s;
    }
    sinon si(config.sortie == LFSRConfig::POIDS_FAIBLE)
      sortie = reg & 1;
    sinon
      sortie = somme;

    bs.push(sortie);
    reg = (reg >> 1) | (somme << (config.reglen - 1));
  }
}


struct LFSRRecepteur::Impl
{
  uint16_t reglen = 0;
  uint32_t reg = 0, pol = 0;
  enum
  {
    PRBS_STATE_UNLOCKED = 0,
    PRBS_STATE_LOCKED
  } state = PRBS_STATE_UNLOCKED;

  uint32_t  nb_consecutive_bits_errors = 0,
            nb_consecutive_bits_ok = 0,
            nb_bits_a_ignorer = 0;
  uint64_t bit_counter = 0;
  bouléen est_verouille = non;
  // Number of bits decoded (right or bad) since the last change in locked state
  uint64_t nb_bits = 0;
  // Number of errors detected since the last change in locked state
  uint64_t nb_erreurs = 0;

  void step(const BitStream &bs)
  {
    soit nbits = bs.lon();
    pour(auto i = 0; i < nbits; i++)
    {
      uint8_t digital_bit = bs[i];
      //uint8_t digital_bit = (v[i] >> (7-j)) & 1;

      soit tmp = reg & pol;
      soit sum = 0u;
      pour(auto k = 0u; k < reglen; k++)
        sum = sum ^ ((tmp >> k) & 1);

      /* Unlocked state */
      si(state == Impl::PRBS_STATE_UNLOCKED)
      {
        /* Fill register with received data */
        reg = (reg >> 1) | (digital_bit << (reglen - 1));

        si(digital_bit == sum)
        {
          /* Bit received == processed: no error */
          nb_consecutive_bits_ok++;
        }
        sinon
        {
          /* Bit received != processed: error */
          nb_consecutive_bits_ok = 0;
        }
        si(nb_consecutive_bits_ok > 20)
        {
          //printf("prbs lock @%Ld.\n", bit_counter);
          /* 8 consecutive bits are good: lock the receiver */
          state = PRBS_STATE_LOCKED;
          nb_consecutive_bits_errors = 0;
          //nb_bits = 0;
          //nb_errors = 0;
          est_verouille = oui;
        }
      }
      /* Locked state */
      sinon
      {
        /* Fill register with computed data */
        reg = (reg >> 1) | (sum << (reglen - 1));
        nb_bits++;

        si(digital_bit == sum)
        {
          /* Bit received == processed: no error */
          nb_consecutive_bits_errors = 0;
        }
        sinon
        {
          /* Bit received != processed: error */
          nb_consecutive_bits_errors++;
          nb_erreurs++;
          //printf("Error @bit %Ld.\n", bit_counter);
        }
        si(nb_consecutive_bits_errors > 5)
        {
          //printf("prbs unlock @%Ld after %Ld bits (%Ld errors).\n",
              //bit_counter, nb_bits, nb_errors);
          /* 5 consecutive bits are bad: unlock the receiver */
          state = PRBS_STATE_UNLOCKED;
          nb_consecutive_bits_ok = 0;
        }
      }
      bit_counter++;
    } // pour(i) : bits
  }
};

LFSRRecepteur::LFSRRecepteur()
{
  impl = std::make_shared<Impl>();
  configure(16, 30);
}

entier LFSRRecepteur::configure(uint16_t reglen, entier nb_bits_to_ignore)
{
  assertion_msg(reglen <= MAX_REGLEN, "Invalid PRBS reg len: {}.", reglen);

  impl->nb_bits_a_ignorer = nb_bits_to_ignore;
  impl->reglen            = reglen;
  impl->pol               = polynome_primitif_binaire(reglen);
  impl->reg               = 0x0001;
  impl->bit_counter       = 0;
  reset();

  retourne 0;
}

void LFSRRecepteur::reset()
{
  impl->state = Impl::PRBS_STATE_UNLOCKED;
  impl->nb_consecutive_bits_ok = 0;
  impl->nb_bits = 0;
  impl->nb_erreurs = 0;
  impl->est_verouille = non;
}

void LFSRRecepteur::step(const BitStream &bs)
{
  impl->step(bs);
}

void LFSRRecepteur::affiche_etat() const
{
  bouléen is_locked;
  float ber;
  lis_etat(is_locked, ber);
  msg("PRBS récepteur : vérouillé = {}, ber = {:.1e}", is_locked ? "oui" : "non", ber);
  msg("  nb bits = {}, nb erreurs = {}", (entier) impl->nb_bits, (entier) impl->nb_erreurs);
}

void LFSRRecepteur::lis_etat(bouléen &is_locked, float &ber) const
{
  is_locked = (impl->state == Impl::PRBS_STATE_LOCKED);
  si((impl->nb_bits == 0) || (!is_locked) || (impl->nb_bits < 10))
    ber = 0.5;
  sinon
    ber = ((float) impl->nb_erreurs) / ((float) impl->nb_bits);
}





}

