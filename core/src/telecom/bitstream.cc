#include "tsd/telecom/bitstream.hpp"
#include <cstdio>
#include <random>

namespace tsd::telecom {


BitStream randstream(entier n)
{
  retourne BitStream::rand(n);
}

BitStream BitStream::zéros(entier nbits)
{
  BitStream r;
  // TODO : plus efficace
  pour(auto i = 0; i < nbits; i++)
    r.push(0);
  retourne r;
}

BitStream BitStream::altern(entier nbits)
{
  BitStream r;

  pour(auto i = 0; i < nbits + 1; i += 2)
  {
    r.push(1);
    r.push(0);
  }

  si(nbits & 1)
    r.push(1);

  retourne r;
}

BitStream BitStream::uns(entier nbits)
{
  BitStream r;
  // TODO : plus efficace
  pour(auto i = 0; i < nbits; i++)
    r.push(1);
  retourne r;
}

BitStream BitStream::rand(entier nbits)
{
  BitStream r;
  r.resize(nbits);
  //std::mt19937 gene;
  std::uniform_int_distribution<> dis(0, 255);
  pour(auto j = 0u; j < r.buffer.size(); j++)
    r.buffer[j] = dis(generateur_aleatoire);
  retourne r;
}

void BitStream::set(unsigned int index, bouléen valeur)
{
  soit octet = index / 8;
  soit bit   = index & 7;

  si(octet >= buffer.size())
    échec("Ecriture bitstream : index invalide ({}, nb bits = {}).", index, buffer.size() * 8);

  si(valeur)
    buffer[octet] |= (1 << bit);
  sinon
    buffer[octet] &= ~(1 << bit);
}

bouléen BitStream::operator [](unsigned int index) const
{
  soit octet = index / 8;
  soit bit   = index & 7;

  si(octet >= buffer.size())
  {
    msg_erreur("Lecture bitstream : index invalide ({}, nb bits = {}, buffer size = {} bits).",
        index, nbits, buffer.size() * 8);
    retourne non;
  }

  retourne buffer[octet] & (1 << bit) ? oui : non;
}

BitStream::BitStream(cstring s)
{
  pour(auto i = 0u; i < s.size(); i++)
  {
    si(s[i] == '1')
      push(oui);
    sinon si(s[i] == '0')
      push(non);
    sinon
      échec("Construction bitstream : chaine invalide \"{}\"", s);
  }
}

entier BitStream::dst_Hamming(const BitStream &bs) const
{
  soit n = lon();
  si(bs.lon() != n)
    échec("Distance de Hamming entre bitstream : dimensions différentes ({} vs {} bits).", n, bs.lon());

  soit cnt = 0;
  pour(auto i = 0; i < n; i++)
    si((*this)[i] != bs[i])
      cnt++;

  retourne cnt;
}

BitStream::BitStream(entier n)
{
  resize(n);
}

BitStream::BitStream(const Vecf &x)
{
  soit nbits = x.rows();
  buffer.reserve((nbits + 7) / 8);
  pour(auto i = 0; i < nbits; i++)
    push(x(i) != 0.0f);
}

Veci BitStream::iarray() const
{
  soit n = lon();
  Veci res(n);
  pour(auto i = 0; i < n; i++)
    res(i) = (*this)[i] ? 1 : 0;
  retourne res;
}

Vecf BitStream::array() const
{
  soit n = lon();
  Vecf res(n);
  pour(auto i = 0; i < n; i++)
    res(i) = (*this)[i] ? 1.0f : 0.0f;
  retourne res;
}

entier BitStream::lon() const
{
  retourne nbits;
}

bouléen BitStream::operator ==(const BitStream &t2)
{
  retourne buffer == t2.buffer;
}

BitStream operator +(const BitStream &t1, const BitStream &t2)
{
  BitStream res;
  pour(auto i = 0; i < t1.lon(); i++)
    res.push(t1[i]);
  pour(auto i = 0; i < t2.lon(); i++)
    res.push(t2[i]);
  retourne res;
}

std::ostream& operator<<(std::ostream &ss, const BitStream &t)
{
  pour(auto i = 0; i < t.lon(); i++)
    ss << ((entier) t[i]);
  retourne ss;
}

void BitStream::restart()
{
  pos = 0;
}

void BitStream::clear()
{
  resize(0);
}

unsigned char *BitStream::get_ptr()
{
  retourne &(buffer[0]);
}

void BitStream::pad(entier nzeros)
{
  pour(auto i = 0; i < nzeros; i++)
    push(0);
}

void BitStream::pad_mult(entier m)
{
  soit nbits = lon();
  soit reste = nbits % m;

  // Exemples :
  //  - k = 1 => RAS
  //  - k = 2, nbits = 11 => push 1 bit
  //
  si(reste != 0)
    pad(m - reste);
}

void BitStream::operator +=(const BitStream &t)
{
  *this = *this + t;
}

void BitStream::resize(entier nbits)
{
  this->nbits = nbits;
  buffer.resize((nbits + 7) / 8, 0);
  pos = 0;
}

void BitStream::push_u32(uint32_t i)
{
  pour(auto j = 0; j < 32; j++)
    push(i & (1u << j));
}

uint32_t BitStream::pop_u32()
{
  uint32_t res = 0;

  pour(auto j = 0; j < 32; j++)
      res |= ((uint32_t) pop() & 1) << j;

  retourne res;
}

void BitStream::push(bouléen b)
{
  nbits++;
  si(nbits > (entier) (8 * buffer.size()))
    buffer.resize((nbits + 7) / 8, 0);

  set(nbits-1, b);
}

bouléen BitStream::pop()
{
  retourne (*this)[pos++];
} 

bouléen BitStream::eof()
{
  retourne pos >= nbits;
}


}
