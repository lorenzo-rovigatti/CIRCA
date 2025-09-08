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

}  // namespace circa
