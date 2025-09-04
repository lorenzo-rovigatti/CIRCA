#pragma once
#include <vector>

#include "../core/field_store.hpp"

namespace circa {

// Diagonal, constant M_i
template <int D>
struct MobilityDiagConst {
    std::vector<double> M;  // size N
    inline double M_i(int i_species, int /*idx*/, const FieldStore<D>& /*S*/) const {
        return M[i_species];
    }
};

// Full constant matrix M_{iβ}
template <int D>
struct MobilityFullConst {
    std::vector<std::vector<double>> M;  // size N x N
    inline double M_ibeta(int i, int b, int /*idx*/, const FieldStore<D>& /*S*/) const {
        return M[i][b];
    }
};

}  // namespace circa
