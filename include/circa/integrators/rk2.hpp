#pragma once
#include "integrator.hpp"

namespace circa {

struct RK2Options {
    bool conservative_mass_fix = false;
    std::string mass_field = "phi";
};

template <int D>
struct RK2 : public IIntegrator<D> {
    RK2Options opt;

    explicit RK2(const BuildFn<D>& build, FieldStore<D>& S0, RK2Options o) : IIntegrator<D>(build, S0), opt(std::move(o)) {
        
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
