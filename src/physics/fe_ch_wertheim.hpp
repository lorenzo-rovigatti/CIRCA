#pragma once
#include "../util/toml.hpp"

namespace circa {

struct FE_CH_Wertheim {
    double B2, delta, two_valence_delta;
    int valence;

    FE_CH_Wertheim(const toml::table &fe_tbl) {
        B2 = *value_or_die<double>(fe_tbl, "B2");
        delta = *value_or_die<double>(fe_tbl, "delta");
        valence = *value_or_die<double>(fe_tbl, "valence");
        two_valence_delta = 2.0 * valence * delta;
    }

    inline double X(double rho) const {
        return (-1.0 + std::sqrt(1.0 + 2.0 * two_valence_delta * rho)) / (two_valence_delta * rho);
    }

    inline double bulk(double rho) const {
        double f_ref = rho * std::log(rho) - rho + B2 * rho * rho;
        double f_bond = (rho > 0.) ? valence * rho * (std::log(X(rho)) + 0.5 * (1. - X(rho))) : 0.0;

        return f_ref + f_bond;
    }

    inline double mu(double rho) const {
        double der_f_ref = std::log(rho) + 2 * B2 * rho;
        double der_f_bond = (rho > 0.) ? valence * std::log(X(rho)) : 0.0;
        
        return der_f_ref + der_f_bond;
    }

    inline double dmu_drho(double rho) const {
        double my_X = X(rho);
        double d2f_ref = 1.0 / rho + 2.0 * B2;
        double d2f_bond = (rho > 0.) ? valence * (my_X - 1.0) / ((2.0 - my_X) * rho) : 0.0;

        return d2f_ref + d2f_bond;
    }
};

}  // namespace circa
