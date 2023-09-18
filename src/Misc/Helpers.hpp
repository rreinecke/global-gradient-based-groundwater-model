/*
 * Copyright (c) <2016>, <Robert Reinecke>
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef GLOBAL_FLOW_HELPERS_HPP
#define GLOBAL_FLOW_HELPERS_HPP

#include <boost/timer/timer.hpp>

#include "colors.hpp"

#include <boost/preprocessor.hpp>

#define X_DEFINE_ENUM_WITH_STRING_CONVERSIONS_TOSTRING_CASE(r, data, elem)    \
    case elem : return BOOST_PP_STRINGIZE(elem);

#define DEFINE_ENUM_WITH_STRING_CONVERSIONS(name, enumerators)                \
    enum name {                                                               \
        BOOST_PP_SEQ_ENUM(enumerators)                                        \
    };                                                                        \
                                                                              \
    inline const char* ToString(name v)                                       \
    {                                                                         \
        switch (v)                                                            \
        {                                                                     \
            BOOST_PP_SEQ_FOR_EACH(                                            \
                X_DEFINE_ENUM_WITH_STRING_CONVERSIONS_TOSTRING_CASE,          \
                name,                                                         \
                enumerators                                                   \
            )                                                                 \
            default: return "[Unknown " BOOST_PP_STRINGIZE(name) "]";         \
        }                                                                     \
    }

#include "../Logging/Logging.hpp"

class NANInSolutionException : public std::exception {
    virtual const char *what() const throw() { return "NaN value in result"; }
};

class InfInSolutionException : public std::exception {
    virtual const char *what() const throw() { return "Inf value in result"; }
};

inline void NANChecker(const double &value, std::string message) {
    if (std::isnan(value)) {
        LOG(GlobalFlow::critical) << "NAN value! :((" << message << "\n";
        throw new NANInSolutionException();
    }
    if (std::isinf(value)) {
        LOG(GlobalFlow::critical) << "INF value! :((" << message << "\n";
        throw new InfInSolutionException();
    }
}

template<typename T>
T const pi = std::acos(-T(1));

/*
 * Helper function for rounding double values up
 */
inline double
roundValue(double valueToRound)
{
    return ceil(valueToRound * 100) / 100;
}

/** http://stackoverflow.com/questions/15181579 **/
template<class T>
struct Is
{
    T d_;
    bool
    in(T a)
    {
        return a == d_;
    }
    template<class Arg, class... Args>
    bool
    in(Arg a, Args... args)
    {
        return in(a) || in(args...);
    }
};

template<class T>
Is<T>
is(T d)
{
    return Is<T>{d};
}

class Position
{
public:
    Position(double lat, double lon) : lat(lat), lon(lon) {}

    const double lat{0};
    const double lon{0};
};


template<typename C>
class Singleton {
public:
    static C *instance() {
        if (!_instance) {
            _instance = new C();
        }
        return _instance;
    }

    virtual
    ~Singleton() {
        _instance = 0;
    }

private:
    static C *_instance;
protected:
    Singleton() {}
};

template<typename C> C *Singleton<C>::_instance = 0;

//Usage:
/*class LoggerInterface : public Singleton<LoggerInterface> {
    friend class Singleton<LoggerInterface>;
private:
public:
    virtual ~LoggerInterface() {};
};*/

#endif