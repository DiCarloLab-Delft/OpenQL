/** \file
 * Provides utilities for working with strings that the STL fails to
 * satisfactorily provide.
 */

#pragma once

#include <string>
#include <sstream>
#include <typeinfo>
#include "utils/num.h"

namespace ql {
namespace utils {

/**
 * Shorthand for std::string.
 */
using Str = std::string;

/**
 * Shorthand for std::ostringstream.
 */
using StrStrm = std::ostringstream;

// Anonymous namespace to make the SFINAE and overload resolution of the
// debug_str_internal() function "private". Use debug_str().
namespace {

template <typename T, typename = decltype(std::declval<std::ostream&>() << std::declval<T>())>
Str try_to_string_internal(T val) {
    StrStrm ss{};
    ss << val;
    return ss.str();
}

template <typename T, typename... X>
Str try_to_string_internal(T val, X... vals) {
    std::ostringstream ss{};
    ss << "<unknown value of type " << typeid(T).name() << ">";
    return ss.str();
}

} // anonymous namespace

/**
 * If a << overload for streams exists for type T, return the result of said
 * overload as a string. If no overload exists, "<unknown value of type ...>"
 * is returned, where ... is derived using typeid::name().
 */
template <typename T>
Str try_to_string(T val) {
    return try_to_string_internal(val);
}

/**
 * Convert the given value to a string using its stream operator<< overload.
 */
template <typename T>
Str to_string(T arg) {
    StrStrm ss;
    ss << arg;
    return ss.str();
}

UInt parse_uint(const Str &str);
Int parse_int(const Str &str);
Real parse_real(const Str &str);

Str to_lower(Str str);
Str replace_all(Str str, const Str &from, const Str &to);

} // namespace utils
} // namespace ql
