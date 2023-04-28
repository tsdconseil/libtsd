#ifndef COMMUN_HPP
#define COMMUN_HPP


#include <string>
#include <vector>
#include <cstdint>
#include <memory>

using std::string;
using std::vector;
using std::tuple;

using cstring = const string &;


template<typename T>
  using sptr = std::shared_ptr<T>;

template<typename T>
  using wptr = std::weak_ptr<T>;

#include <fmt/core.h>
#include <fmt/ostream.h>



#ifndef FMT_RUNTIME
#if (FMT_VERSION >= 80000)
# define FMT_RUNTIME(AA) fmt::runtime(AA)
#else
# define FMT_RUNTIME(AA) AA
#endif
#endif

template<typename ... Ts>
  string sformat(fmt::format_string<Ts...> fmt, Ts &&... args)
{
  return fmt::format(fmt, std::forward<Ts>(args)...);
}

#include "fr.hpp"


#define _PIMPL_ private: struct Impl; sptr<Impl> impl;


#endif
