#ifndef MUSCL_LIMITER_HH
#define MUSCL_LIMITER_HH

#include <cmath>

struct gminmod {
  static double limit(double qL, double qC, double qR) {
    const double GMMLTheta{1.5};
    const double sC{0.5 * (qR - qL)};
    const double sL{GMMLTheta * (qC - qL)};
    const double sR{GMMLTheta * (qR - qC)};
    return (sL * sR < 0.0
              ? 0.0
              : (sC > 0.0 ? 1.0 : -1.0) *
                  std::min(std::abs(sC), std::min(std::abs(sL), std::abs(sR))));
  } // operator()
}; // gminmod

#endif // MUSCL_LIMITER_HH
