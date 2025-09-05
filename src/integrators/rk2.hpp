#pragma once
#include "integrator.hpp"
#include "../util/config.hpp"

namespace circa {

template <int D>
struct RK2 : public IIntegrator<D> {
    explicit RK2(const BuildSysFn<D>& build, FieldStore<D>& S0, const cfg::GeneralConfig<D> &config) : IIntegrator<D>(build, S0) {
        
    }

    void step(FieldStore<D>& S, double dt) override {
        FieldStore<D> k1(S.g);
        this->sys_.set_state(&S, &k1);
        this->sys_.rhs();

        FieldStore<D> S_tmp = plus_scaled(S, k1, 1.0, dt);

        FieldStore<D> k2(S.g);
        k2.zero();
        this->sys_.set_state(&S_tmp, &k2);
        this->sys_.rhs();

        FieldStore<D> sum = plus_scaled(k1, k2, 1.0, 1.0);
        axpy(S, sum, 0.5 * dt);
    }
};

}  // namespace circa
