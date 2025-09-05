#pragma once
#include "integrator.hpp"
#include "../util/config.hpp"

namespace circa {

template <int D>
struct Euler : public IIntegrator<D> {
    explicit Euler(const BuildSysFn<D>& build, FieldStore<D>& S0, const cfg::GeneralConfig<D> &config) : IIntegrator<D>(build, S0) {}

    void step(FieldStore<D>& S, double dt) override {
        FieldStore<D> k1(S.g);
        k1.zero();
        this->sys_.set_state(&S, &k1);
        this->sys_.rhs();

        axpy(S, k1, dt);
    }
};

}  // namespace circa
