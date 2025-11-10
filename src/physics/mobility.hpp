#pragma once
#include <cmath>
#include <string>

#include "../core/field_store.hpp"

namespace circa {

template <int D>
struct MobConst {
    double M0 = 1.0;

    inline double operator()(int /*i*/, const FieldStore<D>& /*S*/) const { 
        return M0;
    }
};

template <int D>
struct MobExpOfField {
    std::string field;
    double c0;

    inline double operator()(int i, const FieldStore<D>& S) const { 
        return std::exp(-S.get(field).a[i] / c0);
    }
};

// Mobility model based on Wertheim theory
// First, a "tag" type to hold configuration
template<int D>
struct MobWertheimAuto {
    double D0 = 1.0;
    std::string field = "c";
};

// And then the actual mobility model that uses a FE_CH_Wertheim instance
template<int D, class FE>
struct MobWertheimBound {
    MobWertheimAuto<D> cfg;
    const FE* fe = nullptr;

    inline double operator()(int i, const FieldStore<D>& S) const {
        const double rho = S.get(cfg.field).a[i];
        const double dmu_drho = fe->dmu_drho(rho);
        const double X = fe->X(rho);
        return cfg.D0 * std::pow(X, fe->valence) / dmu_drho;
    }   
};


}  // namespace circa
