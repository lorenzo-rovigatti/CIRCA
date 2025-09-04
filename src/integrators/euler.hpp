#pragma once
#include "../core/field_store.hpp"
#include "integrator.hpp"

namespace circa {

struct EulerOptions {
    bool conservative_mass_fix = false;
    std::string mass_field = "phi";
};

template <int D>
struct Euler : public IIntegrator<D> {
    EulerOptions opt;
    explicit Euler(const BuildFn<D>& build, FieldStore<D>& S0, EulerOptions o) : IIntegrator<D>(build, S0), opt(std::move(o)) {}

    void step(FieldStore<D>& S, double dt) override {
        FieldStore<D> k1(S.g);
        k1.zero();
        this->sys_.set_state(&S, &k1);
        this->sys_.rhs();

        axpy(S, k1, dt);

        if (opt.conservative_mass_fix) {
            try {
                auto& f = S.map.at(opt.mass_field);
                double m = mean(f);
                for (int i = 0; i < f.g.size; ++i) f.a[i] -= m;
            } 
            catch (...) {
            }
        }
    }
};

}  // namespace circa
