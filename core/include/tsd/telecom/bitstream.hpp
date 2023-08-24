#ifndef TSD_BITSTREAM_H
#define TSD_BITSTREAM_H

#include "tsd/tsd.hpp"

namespace tsd::telecom {

 /** @addtogroup telecom
  *  @{
  */


// En écriture : on pousse des bits un par un, ré-allocation automatique au fur et à mesure du remplissage
// En lecture: on lis un par un, jusqu'à la fin

/** @brief Chaine binaire
 *
 * Cette classe facilite la manipulation de séquences binaires.
 */
class BitStream
{
public:

  /** @brief Construction d'une chaine vide */
  BitStream(){}

  /** @brief Constructeur, à partir d'un vecteur de flottants.
   *
   *  Les éléments non nuls de x sont interprétés comme des 1.
   */
  BitStream(const Vecf &x);

  /** @brief Constructeur, d'après la dimension (nombre de bits) */
  BitStream(entier n);

  /** @brief Constructeur, d'après une chaine de caractères de type "0100111001..." */
  BitStream(cstring s);

  entier dst_Hamming(const BitStream &bs) const;

  /** @brief Conversion vers vecteur de flottants (valeurs = 0 ou 1). */
  Vecf array() const;

  /** @brief Conversion vers vecteur d'entiers. */
  Veci iarray() const;

  /** @brief Lecture d'un bit d'index donné */
  bouléen operator [](unsigned int index) const;

  /** @brief Ecriture d'un bit d'index donné */
  void set(unsigned int index, bouléen valeur);

  /** @brief Train de bit pseudo-aléatoire */
  static BitStream rand(entier nbits);

  /** @brief Zéros */
  static BitStream zéros(entier nbits);

  /** @brief Uns */
  static BitStream uns(entier nbits);

  /** @brief Alternance de 0 et de 1 */
  static BitStream altern(entier nbits);

  void restart();

  /** @brief Ré-allocation de la dimension.
   *
   * @warning Remet à zéro tous les bits. */
  void resize(entier nbits);

  /** @brief Insertion d'un nouveau bit en fin de séquence. */
  void push(bouléen b);

  /** @brief Insertion d'un entier 32 bits en fin de séquence. */
  void push_u32(uint32_t i);

  /** @brief Retire un entier 32 bits du début de la séquence. */
  uint32_t pop_u32();

  /** @brief Lecture et dépilement d'un bit en début de séquence. */
  bouléen pop();

  /** @brief Fin de séquence atteinte (en lecture) ? */
  bouléen eof();

  /** @brief Nombre de bits */
  entier lon() const;

  /** @brief Suppression de tous les bits.
   *
   * La séquence est vide après l'appel à cette méthode.
   */
  void clear();

  /** @brief Concaténation avec une autre chaîne binaire */
  void operator +=(const BitStream &t);

  /** @brief Comparaison */
  bouléen operator ==(const BitStream &t2);

  /** @cond undoc */
  unsigned char *get_ptr();
  /** @endcond */

  /** @brief Ajoute des zéros à la fin */
  void pad(entier nzeros);

  /** @brief Ajoute (si besoin) des zéros de manière à ce que le nombre de bits soit un multiple de m */
  void pad_mult(entier m);

private:
  vector<unsigned char> buffer;
  entier pos = 0, nbits = 0;
};



/** @cond undoc */
extern std::ostream& operator<<(std::ostream &ss, const BitStream &t);
/** @endcond */

/** @brief Concaténation de deux chaines binaires. */
extern BitStream operator +(const BitStream &t1, const BitStream &t2);


/** @brief Chaîne binaire pseudo-aléatoire */
extern BitStream randstream(entier n);

/** @} */

}


ostream_formater(tsd::telecom::BitStream)


#endif

