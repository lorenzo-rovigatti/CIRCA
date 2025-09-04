#pragma once
#include <array>

namespace circa {

template <int D>
struct Grid {
    std::array<int, D> n{};
    std::array<double, D> L{};
    std::array<double, D> dx{};
    double dV;
    int size = 0;

    Grid() = default;
    Grid(const std::array<int, D>& n_, const std::array<double, D>& L_) : n(n_), L(L_) {
        size = 1;
        dV = 1.0;
        for (int d = 0; d < D; ++d) {
            dx[d] = L[d] / n[d];
            dV *= dx[d];
            size *= n[d];
        }
    }
};

template <int D>
inline int flat(const std::array<int, D>& I, const std::array<int, D>& n) {
    int s = 0, m = 1;
    for (int d = 0; d < D; ++d) {
        s += I[d] * m;
        m *= n[d];
    }
    return s;
}
template <int D>
inline std::array<int, D> unflat(int lin, const std::array<int, D>& n) {
    std::array<int, D> I{};
    for (int d = 0; d < D; ++d) {
        I[d] = lin % n[d];
        lin /= n[d];
    }
    return I;
}

}  // namespace circa
