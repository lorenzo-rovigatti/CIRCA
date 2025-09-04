#include <iostream>
#include <random>

#include "core/diagnostics.hpp"
#include "core/field_store.hpp"
#include "core/grid.hpp"
#include "core/system.hpp"
#include "integrators/registry.hpp"
#include "io/log.hpp"
#include "io/plain.hpp"
#include "ops/fd_ops.hpp"
#include "physics/fe_ac_linear.hpp"
#include "physics/fe_ch_landau.hpp"
#include "physics/mobility.hpp"
#include "terms/ac_term.hpp"
#include "terms/ch_term.hpp"
#include "util/config.hpp"

using namespace circa;

int main(int argc, char *argv[]) {
    if(argc < 2) {
		std::cerr << fmt::format("Usage is {} configuration_file", argv[0]) << std::endl;
		return 0;
	}

    circa::log::init();

    circa::cfg::GeneralConfig<DIM> config; // this instance contains the TOML table storing all the "raw" options,
                                           // which is passed around and therefore have to remain alive over the course 
                                           // of the simulation.
    
    try {
        config = circa::cfg::load<DIM>(argv[1]);
    }
    catch (std::runtime_error &e) {
        CIRCA_CRITICAL(e.what());
        exit(1);
    }

    Grid<DIM> grid(config.grid.n, config.grid.L);
    FieldStore<DIM> S(grid);
    FieldStore<DIM> scratch(grid); // used later as a placeholder

    for(const auto &n : config.fields.names) {
        S.ensure(n);
    }

    std::mt19937 rng(42);
    std::normal_distribution<double> nphi(0.0, 0.05);
    for(int i = 0; i < grid.size; ++i) {
        S.map["phi"].a[i] = nphi(rng);
        S.map["c"].a[i] = 0.0;
    }

    ParsedInput cfg;
    cfg.integrator = "euler";
    cfg.mass_fix = false;
    auto registry = make_integrator_registry<DIM>();
    auto it = registry.find(config.integrator.name);
    if(it == registry.end()) {
        CIRCA_CRITICAL("Unknown integrator");
        exit(0);
    }
    auto stepper = it->second(cfg, config.build_system_fn, S);

    circa::io::dump_all_fields_plain<DIM>(S, "init", 0, 0.0, false);

    for(int64_t s = 0; s <= config.time.steps; s++) {
        double t = s * config.time.dt;
        if(s % config.out.output_every == 0) {
            const auto& phi = S.get("phi");
            double m_avg = mean(phi);
            double m_tot = circa::Diagnostics<DIM>::total_mass(phi);

            auto sys_now = config.build_system_fn(S, scratch);
            double FE = circa::Diagnostics<DIM>::total_free_energy(sys_now);
            std::cout << fmt::format("{:.5} {:.8} {:.5} {:L}", t, FE, m_avg, s) << std::endl;
        }
        if(s % config.out.conf_every == 0) {
            circa::io::dump_all_fields_plain<DIM>(S, "last", s, t, false);
        }
        stepper->step(S, config.time.dt);
    }
    CIRCA_INFO("END OF SIMULATION");
    return 0;
}
