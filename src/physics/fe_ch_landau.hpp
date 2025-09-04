#pragma once
#include "../util/toml.hpp"

namespace circa {
struct FE_CH_Landau {
    double eps;

    FE_CH_Landau(const toml::table &fe_tbl) {
        eps = *value_or_die<double>(fe_tbl, "eps");
    }

    inline double bulk(double u) const {
        return -(eps*0.5)*u*u + 0.25*u*u*u*u;
    }

    inline double mu(double u) const { 
        return -eps * u + u * u * u;
    }
};
}  // namespace circa
