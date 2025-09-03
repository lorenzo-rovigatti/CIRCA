#pragma once
#include "system.hpp"

namespace circa {

template <int D>
struct Diagnostics {
    // total mass of a scalar field (integral ∫ u dV)
    static double total_mass(const Field<D>& f) {
        double sum = 0.0;
        for(double v : f.a) {
            sum += v;
        }
        return sum * f.g.dV;
    }
    static double mean_value(const Field<D>& f) { return mean(f); }  // already had this

    // ask each term that implements IEnergy<D>
    static double total_free_energy(const System<D>& sys) {
        double E = 0.0;
        for (const auto& up : sys.terms) {
            if (auto e = dynamic_cast<const IEnergy<D>*>(up.get())) E += e->energy();
        }
        return E;
    }
};

}  // namespace circa
