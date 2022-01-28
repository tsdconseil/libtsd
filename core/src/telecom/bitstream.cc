#include "tsd/telecom/bitstream.hpp"
#include <cstdio>
#include <random>

namespace tsd::telecom {


BitStream randstream(int n)
{
  return BitStream::rand(n);
}

BitStream BitStream::zéros(int nbits)
{
  BitStream r;
  // TODO : plus efficace
  for(auto i = 0; i < nbits; i++)
    r.push(0);
  return r;
}

BitStream BitStream::altern(int nbits)
{
  BitStream r;

  for(auto i = 0; i < nbits + 1; i += 2)
  {
    r.push(1);
    r.push(0);
  }

  if(nbits & 1)
    r.push(1);

  return r;
}

BitStream BitStream::uns(int nbits)
{
  BitStream r;
  // TODO : plus efficace
  for(auto i = 0; i < nbits; i++)
    r.push(1);
  return r;
}

BitStream BitStream::rand(int nbits)
{
  BitStream r;
  r.resize(nbits);
  //std::default_random_engine gene;
  std::mt19937 gene;
  std::uniform_int_distribution<> dis(0, 255);
  for(auto j = 0u; j < r.buffer.size(); j++)
  {
    r.buffer[j] = dis(gene);
    //msg("rand[{}] = {}", j, r.buffer[j]);
  }
  //r.pos = nbits;
  //for(auto j = 0u; j < n; j++)
    //r.set(j, dis(generator));
  return r;
}

void BitStream::set(unsigned int index, bool valeur)
{
  auto byte = index / 8;
  auto bit  = index & 7;

  if(byte >= buffer.size())
    echec("Ecriture bitstream : index invalide ({}, nb bits = {}).", index, buffer.size() * 8);

  if(valeur)
    buffer[byte] |= (1 << bit);
  else
    buffer[byte] &= ~(1 << bit);
}

bool BitStream::operator [](unsigned int index) const
{
  auto byte = index / 8;
  auto bit  = index & 7;


  if(byte >= buffer.size())
  {
    msg_erreur("Lecture bitstream : index invalide ({}, nb bits = {}, buffer size = {} bits).",
        index, nbits, buffer.size() * 8);
    return false;
  }

  return buffer[byte] & (1 << bit) ? true : false;
}

BitStream::BitStream(const std::string &s)
{
  for(auto i = 0u; i < s.size(); i++)
  {
    if(s[i] == '1')
      push(true);
    else if(s[i] == '0')
      push(false);
    else
      echec("Construction bitstream : chaine invalide \"{}\"", s);
  }
}

int BitStream::dst_Hamming(const BitStream &bs) const
{
  auto n = lon();
  if(bs.lon() != n)
    echec("Distance de Hamming entre bitstream : dimensions différentes ({} vs {} bits).", n, bs.lon());

  int cnt = 0;
  for(auto i = 0; i < n; i++)
    if((*this)[i] != bs[i])
      cnt++;

  return cnt;
}

BitStream::BitStream(int n)
{
  resize(n);
}

BitStream::BitStream(IArrayXf x)
{
  auto nbits = x.rows();
  buffer.reserve((nbits + 7) / 8);
  for(auto i = 0; i < nbits; i++)
    push(x(i) != 0.0f);
}

Eigen::ArrayXi BitStream::iarray() const
{
  auto n = lon();
  Eigen::ArrayXi res(n);
  for(auto i = 0; i < n; i++)
    res(i) = (*this)[i] ? 1 : 0;
  return res;
}

ArrayXf BitStream::array() const
{
  auto n = lon();
  ArrayXf res(n);
  for(auto i = 0; i < n; i++)
    res(i) = (*this)[i] ? 1.0f : 0.0f;
  return res;
}

int BitStream::lon() const
{
  return nbits;//buffer.size() * 8;
}

bool BitStream::operator ==(const BitStream &t2)
{
  return buffer == t2.buffer;
}

BitStream operator +(const BitStream &t1, const BitStream &t2)
{
  BitStream res;
  for(auto i = 0; i < t1.lon(); i++)
    res.push(t1[i]);
  for(auto i = 0; i < t2.lon(); i++)
    res.push(t2[i]);
  return res;
}

std::ostream& operator<<(std::ostream &ss, const BitStream &t)
{
  for(auto i = 0; i < t.lon(); i++)
    ss << ((int) t[i]);
  return ss;
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
  return &(buffer[0]);
}

void BitStream::pad(int nzeros)
{
  for(auto i = 0; i < nzeros; i++)
    push(0);
}

void BitStream::pad_mult(int m)
{
  int nbits = lon();
  int reste = nbits % m;

  // Exemples :
  //  - k = 1 => RAS
  //  - k = 2, nbits = 11 => push 1 bit
  //
  if(reste != 0)
    pad(m - reste);
}

void BitStream::operator +=(const BitStream &t)
{
  *this = *this + t;
}

void BitStream::resize(int nbits)
{
  this->nbits = nbits;
  buffer.resize((nbits + 7) / 8, 0);
  pos = 0;
}

void BitStream::push_u32(uint32_t i)
{
  for(auto j = 0; j < 32; j++)
    push(i & (1u << j));
}

uint32_t BitStream::pop_u32()
{
  uint32_t res = 0;

  for(auto j = 0; j < 32; j++)
      res |= ((uint32_t) pop() & 1) << j;

  return res;
}

void BitStream::push(bool b)
{
  //if(pos >= (int) (8 * buffer.size()))
    //buffer.resize(buffer.size() + 1, 0);

  nbits++;
  if(nbits > (int) (8 * buffer.size()))
    buffer.resize((nbits + 7) / 8, 0);

  set(nbits-1, b);

  //tsd_assert(!eof());

  /*auto byte = nbits / 8;
  auto bit  = nbits & 7;

  if(b)
    buffer[byte] |= (1 << bit);
  else
    buffer[byte] &= ~(1 << bit);

  //pos++;
  nbits++;*/
}

bool BitStream::pop()
{
  return (*this)[pos++];
} 

bool BitStream::eof()
{
  return pos >= nbits;//(int) (8 * buffer.size());
}


}
