#pragma once
namespace circa {
struct FE_CH_Landau {
    double eps = 0.8, kappa = 1.0;

    inline double bulk(double u) const {
        return -(eps*0.5)*u*u + 0.25*u*u*u*u;
    }

    inline double interfacial(double u) const {
        return 0.5 * kappa;
    }

    inline double mu(double u, double lap_u) const { 
        return -eps * u + u * u * u - kappa * lap_u; 
    }
};
}  // namespace circa
