#if !defined(FR_HPP) || defined(__CDT_PARSER__)
#define FR_HPP




#define soit      auto
#define Soient    auto
#define soient    auto
#define Si        if
#define si        if
#define sinon     else
#define retourne  return
#define retourne  return
#define Pour      for
#define Tantque   while
#define pour      for
#define tantque   while
#define ou        ||
#define et        &&

template<typename T>
  using fonction = std::function<T>;

static const auto
  non = false,
  oui = true;

using bouléen = bool;
using entier  = int;


#endif

