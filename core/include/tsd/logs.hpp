#ifndef CORE2_INCLUDE_TSD_LOGS_HPP_
#define CORE2_INCLUDE_TSD_LOGS_HPP_


//#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <fmt/ostream.h>
//#include <fmt/color.h>
#include <functional>
#include <cassert>

#include "tsd/fr.hpp"

#ifndef FMT_RUNTIME
#if (FMT_VERSION >= 80000)
# define FMT_RUNTIME(AA) fmt::runtime(AA)
#else
# define FMT_RUNTIME(AA) AA
#endif
#endif


namespace tsd {
extern void msg_impl2(const char *fn, const entier ligne, entier niveau, const char *fonction, const std::string &str);


typedef std::function<void(const char *fn, const entier ligne, entier niveau, const char *fonction, const std::string &str)> logger_t;

template<typename ... Ts>
  void msg_impl(const char *fn, const entier ligne, entier niveau, const char *fonction, const std::string &format_str, Ts &&... args)
{
  auto s = fmt::format(FMT_RUNTIME(format_str), args...);
  msg_impl2(fn, ligne, niveau, fonction, s);
}

extern void set_logger(logger_t logger);
extern void reset_logger();

}

#ifndef msg

#define msg_verb(...)     tsd::msg_impl(__FILE__, __LINE__, 0, __PRETTY_FUNCTION__, __VA_ARGS__)
#define msg(...)          tsd::msg_impl(__FILE__, __LINE__, 1, __PRETTY_FUNCTION__, __VA_ARGS__)
#define msg_majeur(...)   tsd::msg_impl(__FILE__, __LINE__, 2, __PRETTY_FUNCTION__, __VA_ARGS__)
#define msg_avert(...)    tsd::msg_impl(__FILE__, __LINE__, 3, __PRETTY_FUNCTION__, __VA_ARGS__)
#define msg_erreur(...)   tsd::msg_impl(__FILE__, __LINE__, 4, __PRETTY_FUNCTION__, __VA_ARGS__)


#endif

#ifndef echec
#define echec(...)        tsd::msg_impl(__FILE__, __LINE__, 5, __PRETTY_FUNCTION__, __VA_ARGS__)
#endif

/** @brief Display an error message in the log before asserting */
#define tsd_assert_msg(AA, ...)  if(!(AA)) {tsd::msg_impl(__FILE__, __LINE__, 5, __PRETTY_FUNCTION__, __VA_ARGS__);}
#define tsd_assert(AA)           if(!(AA)) {tsd::msg_impl(__FILE__, __LINE__, 5, __PRETTY_FUNCTION__, "{}", "Echec assertion : " #AA ".");}

//#define dsp_assertion(AA, ...) {printf("\033[1;37;41mEchec assertion : " #AA ".\nFichier: %s, ligne: %d.\033[0m", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); fflush(0); {*((char *) 0) = 5;} assert(AA);}

/** @brief Returns -1 if the condition is not true (and add an error message in the log) */
//#define tsd_check(AA, ...) if(!(AA)){dsp_assertion(AA, __VA_ARGS__);}

#ifdef DEBUG_MODE
# define tsd_assert_safe tsd_assert_msg
#else
# define tsd_assert_safe(...)
#endif



#endif /* CORE2_INCLUDE_TSD_LOGS_HPP_ */
