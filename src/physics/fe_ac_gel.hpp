#pragma once
#include "../util/toml.hpp"

namespace circa {
struct FE_AC_Gel {
    double critical_OP;
    double M_c;
    double p_gel;
    bool rescale_OP;

    FE_AC_Gel(const toml::table &fe_tbl) {
        critical_OP = *value_or_die<double>(fe_tbl, "critical_OP");
        M_c = *value_or_die<double>(fe_tbl, "M_c");
        p_gel = *value_or_die<double>(fe_tbl, "p_gel");
        rescale_OP = value_or<bool>(fe_tbl, "rescale_OP", true);
    }

    inline double dfdc(double c, double driver) const {
        double phi = (rescale_OP) ? (driver + 1.0) / 2.0 : driver;
        double g = (p_gel * phi - critical_OP) / (1.0 - critical_OP);
        return M_c * (c * c - g * c);
    }
};
}  // namespace circa
