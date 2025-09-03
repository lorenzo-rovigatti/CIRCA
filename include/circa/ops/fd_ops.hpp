#pragma once
#include "../core/grid.hpp"
#include "deriv_ops.hpp"

namespace circa {

template <int D>
struct FDOps : DerivOps<D> {
    Field<D> laplacian(const Field<D>& f) const override {
        Field<D> out(f.g);
        for (int i = 0; i < f.g.size; ++i) {
            auto I = unflat<D>(i, f.g.n);
            double acc = 0.0;
            for (int d = 0; d < D; ++d) {
                auto Ip = I, Im = I;
                Ip[d] = (I[d] + 1 == f.g.n[d]) ? 0 : I[d] + 1;
                Im[d] = (I[d] == 0) ? f.g.n[d] - 1 : I[d] - 1;
                acc += (f.a[flat<D>(Ip, f.g.n)] - 2.0 * f.a[i] + f.a[flat<D>(Im, f.g.n)]) / (f.g.dx[d] * f.g.dx[d]);
            }
            out.a[i] = acc;
        }
        return out;
    }
    std::array<Field<D>, D> gradient(const Field<D>& f) const override {
        std::array<Field<D>, D> g{Field<D>(f.g)};
        for (int d = 1; d < D; ++d) g[d] = Field<D>(f.g);
        for (int i = 0; i < f.g.size; ++i) {
            auto I = unflat<D>(i, f.g.n);
            for (int d = 0; d < D; ++d) {
                auto Ip = I, Im = I;
                Ip[d] = (I[d] + 1 == f.g.n[d]) ? 0 : I[d] + 1;
                Im[d] = (I[d] == 0) ? f.g.n[d] - 1 : I[d] - 1;
                g[d].a[i] = (f.a[flat<D>(Ip, f.g.n)] - f.a[flat<D>(Im, f.g.n)]) / (2.0 * f.g.dx[d]);
            }
        }
        return g;
    }
    Field<D> divergence(const std::array<Field<D>, D>& v) const override {
        Field<D> out(v[0].g);
        for (int i = 0; i < out.g.size; ++i) {
            auto I = unflat<D>(i, out.g.n);
            double acc = 0.0;
            for (int d = 0; d < D; ++d) {
                auto Ip = I, Im = I;
                Ip[d] = (I[d] + 1 == out.g.n[d]) ? 0 : I[d] + 1;
                Im[d] = (I[d] == 0) ? out.g.n[d] - 1 : I[d] - 1;
                acc += (v[d].a[flat<D>(Ip, out.g.n)] - v[d].a[flat<D>(Im, out.g.n)]) / (2.0 * out.g.dx[d]);
            }
            out.a[i] = acc;
        }
        return out;
    }
};

}  // namespace circa
