#if !defined(COMMUN_HPP) || defined(__CDT_PARSER__)
#define COMMUN_HPP


/**
 *   Alias et mise à disposition de composants essentiels du C++
 *   (pour un code plus compact et lisible)
 */

#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <functional>

using std::string;
using std::vector;
using std::tuple;

using cstring = const string &;


template<typename T>
  using sptr = std::shared_ptr<T>;

template<typename T>
  using wptr = std::weak_ptr<T>;


using std::make_shared;


template <class T, class U>
sptr<T> dyncast(const sptr<U> &r) noexcept
{
    return std::dynamic_pointer_cast<T>(r);
}

#include "fr.hpp"

typedef fonction<void(const char *fn, entier ligne, entier niveau, cstring str)> logger_t;


// Peut-être implémenté :
//  - Dans libcutil
//  - Dans libtsd
//  - Dans libcutil+libtsd

inline logger_t &get_logger()
{
  // Variable globale
  static logger_t logger;
  retourne logger;
}


inline void log_msg_impl(const char *fn, entier ligne, entier niveau, cstring str)
{
  si(get_logger())
    get_logger()(fn, ligne, niveau, str);
}

#ifndef FMT_THROW
#define FMT_THROW(AA) log_msg_impl("", 0, 4, "Fmt error.")
#endif

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/color.h>


#ifndef FMT_RUNTIME
#if (FMT_VERSION >= 80000)
# define FMT_RUNTIME(AA) fmt::runtime(AA)
#else
# define FMT_RUNTIME(AA) AA
#endif
#endif

#if OLD_GCC
template<typename ... Ts>
  string sformat(string fmt, Ts &&... args)
{
  return fmt::format(fmt, std::forward<Ts>(args)...);
}
#else
template<typename ... Ts>
  string sformat(fmt::format_string<Ts...> fmt, Ts &&... args)
{
  return fmt::format(fmt, std::forward<Ts>(args)...);
}
#endif






#if OLD_GCC
template<typename ... Ts>
  void log_msg(const char *fn, const entier ligne, entier niveau, string fmt, Ts &&... args)
{
  string s = fmt::format(fmt, args...);
  log_msg_impl(fn, ligne, niveau, s);
}
#else
template<typename ... Ts>
  void log_msg(const char *fn, const entier ligne, entier niveau, fmt::format_string<Ts...> fmt, Ts &&... args)
{
  string s = fmt::format(fmt, std::forward<Ts>(args)...);
  log_msg_impl(fn, ligne, niveau, s);
}
#endif

#ifndef msg
#define msg_verb(...)     log_msg(__FILE__, __LINE__, 0, __VA_ARGS__)
#define msg(...)          log_msg(__FILE__, __LINE__, 1, __VA_ARGS__)
#define msg_majeur(...)   log_msg(__FILE__, __LINE__, 2, __VA_ARGS__)
#define msg_avert(...)    log_msg(__FILE__, __LINE__, 3, __VA_ARGS__)
#define msg_erreur(...)   log_msg(__FILE__, __LINE__, 4, __VA_ARGS__)
#endif


/** Génération d'une exception */
#if OLD_GCC
template<typename ... Ts>
  [[noreturn]] void failure(string fmt, Ts &&... args)
{
  soit s = sformat(fmt, std::forward<Ts>(args)...);
  msg_erreur("{}", s);
  throw s;
}
#else
template<typename ... Ts>
  [[noreturn]] void failure(fmt::format_string<Ts...> fmt, Ts &&... args)
{
  soit s = sformat(fmt, std::forward<Ts>(args)...);
  msg_erreur("{}", s);
  throw s;
}
#endif

#if !OLD_GCC
#define échec failure
#endif


/** @brief Display an error message in the log before asserting */
#define assertion_msg(AA, ...)  do{si(!(AA)) {log_msg(__FILE__, __LINE__, 5, __VA_ARGS__);}} while(0)
#define assertion(AA)           do{si(!(AA)) {log_msg(__FILE__, __LINE__, 5, "{}", "Echec assertion : " #AA ".");}} while(0)

#ifdef DEBUG_MODE
# define assertion_safe assertion_msg
#else
# define assertion_safe(...)
#endif

/** Motif de conception PIMPL (Pointer to Implementaiton) */
#define _PIMPL_ private: struct Impl; sptr<Impl> impl;


#endif
