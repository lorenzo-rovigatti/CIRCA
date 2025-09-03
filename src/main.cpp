#include <iostream>
#include <random>

#include "circa/core/diagnostics.hpp"
#include "circa/core/field_store.hpp"
#include "circa/core/grid.hpp"
#include "circa/core/system.hpp"
#include "circa/io/log.hpp"
#include "circa/io/plain.hpp"
#include "circa/integrators/registry.hpp"
#include "circa/ops/fd_ops.hpp"
#include "circa/physics/fe_ac_linear.hpp"
#include "circa/physics/fe_ch_landau.hpp"
#include "circa/physics/mobility.hpp"
#include "circa/terms/ac_term.hpp"
#include "circa/terms/ch_term.hpp"

using namespace circa;

int main() {
    circa::log::init();

    std::array<int, DIM> grid_cells;
    grid_cells.fill(128);
    std::array<double, DIM> grid_size;
    grid_size.fill(256.0);

    Grid<DIM> grid(grid_cells, grid_size);
    FieldStore<DIM> S(grid);
    FieldStore<DIM> scratch(grid); // used later as a placeholder

    for(auto n : {"phi", "c"}) {
        S.ensure(n);
    }

    std::mt19937 rng(42);
    std::normal_distribution<double> nphi(0.0, 0.05);
    for(int i = 0; i < grid.size; ++i) {
        S.map["phi"].a[i] = nphi(rng);
        S.map["c"].a[i] = 0.0;
    }

    FDOps<DIM> fd;

    auto build_system = [&](FieldStore<DIM>& Sin, FieldStore<DIM>& dS) -> System<DIM> {
        System<DIM> sys;
        sys.add(std::make_unique<CHTerm<DIM, FE_CH_Landau, MobConst<DIM>, FDOps<DIM>>>(
            Sin, dS, fd, "phi", FE_CH_Landau{0.8, 1.0}, MobConst<DIM>{1.0}));
        // sys.add(std::make_unique<ACTerm<DIM, FE_AC_Linear, FDOps<DIM>>>(Sin, dS, fd, "c", "phi", FE_AC_Linear{0.1, 1.0}));
        return sys;
    };

    ParsedInput cfg;
    cfg.integrator = "euler";
    cfg.mass_fix = false;
    auto registry = make_integrator_registry<DIM>();
    auto it = registry.find(cfg.integrator);
    if(it == registry.end()) {
        CIRCA_CRITICAL("Unknown integrator");
        exit(0);
    }
    auto stepper = it->second(cfg, build_system, S);

    double dt = 0.01, t = 0.0;
    int steps = 1000, out_every = 10000, conf_every = 10000;

    circa::io::dump_all_fields_plain<DIM>(S, "init", 0, 0.0, false);

    for(int s = 0; s <= steps; s++) {
        if(s % out_every == 0) {
            const auto& phi = S.get("phi");
            double m_avg = mean(phi);
            double m_tot = circa::Diagnostics<DIM>::total_mass(phi);

            scratch.zero();
            auto sys_now = build_system(S, scratch);
            double FE = circa::Diagnostics<DIM>::total_free_energy(sys_now);
            std::cout << t << " " << FE << " " << m_avg << " " << s << "\n";
        }
        if(s % conf_every == 0) {
            circa::io::dump_all_fields_plain<DIM>(S, "last", s, t, false);
        }
        stepper->step(S, dt);
        t += dt;
    }
    CIRCA_INFO("END OF SIMULATION");
    return 0;
}
