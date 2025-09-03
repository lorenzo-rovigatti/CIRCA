#pragma once
namespace circa {
struct FE_AC_Linear {
    double Mc = 0.1, gcoef = 1.0;
    inline double dfdc(double c, double driver) const { return Mc * (c - gcoef * driver); }
};
}  // namespace circa
