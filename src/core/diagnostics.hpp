#pragma once
#include "system.hpp"

namespace circa {

template <int D>
struct Diagnostics {
    static double total_mass(const Field<D>& f) {
        double sum = 0.0;
        for(double v : f.a) {
            sum += v;
        }
        return sum * f.g.dV;
    }
    
    static double total_free_energy(const System<D>& sys) {
        double FE = 0.0;
        for(const auto& up : sys.terms) {
            if(auto e = dynamic_cast<const IEnergy<D>*>(up.get())) {
                FE += e->energy();
            }
        }
        return FE;
    }
};

}  // namespace circa
