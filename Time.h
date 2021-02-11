#pragma once
#include <time.h> 

namespace SiriusFM {
inline double IntervalYearFrac(time_t t) {
    constexpr double SecY = 365.25*86400.;
    return (double) t/SecY;
}

inline double YearFrac(time_t t) {
    return 1970.+ IntervalYearFrac(t);
}

};
