#pragma once
#include "integrator.hpp"

namespace circa {

struct RK4Options {
    bool enable_filter = false;
};

template <int D>
struct RK4 : IIntegrator<D> {
    RK4Options opt;
    explicit RK4(const BuildFn<D>& build, FieldStore<D>& S0, RK4Options o) : IIntegrator<D>(build, S0), opt(std::move(o)) {}

    void step(FieldStore<D>& S, double dt) override {
        FieldStore<D> k1(S.g), k2(S.g), k3(S.g), k4(S.g);
        k1.zero();
        this->sys_.set_state(&S, &k1);
        this->sys_.rhs();

        auto S2 = plus_scaled(S, k1, 1.0, 0.5 * dt);
        k2.zero();
        this->sys_.set_state(&S2, &k1);
        this->sys_.rhs();

        auto S3 = plus_scaled(S, k2, 1.0, 0.5 * dt);
        k3.zero();
        this->sys_.set_state(&S3, &k1);
        this->sys_.rhs();

        auto S4 = plus_scaled(S, k3, 1.0, dt);
        k4.zero();
        this->sys_.set_state(&S4, &k1);
        this->sys_.rhs();

        auto tmp = plus_scaled(k1, k2, 1.0, 2.0);
        auto tmp2 = plus_scaled(tmp, k3, 1.0, 2.0);
        auto sum = plus_scaled(tmp2, k4, 1.0, 1.0);
        axpy(S, sum, dt / 6.0);
    }
};

}  // namespace circa
