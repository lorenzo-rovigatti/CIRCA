#include "config.hpp"

#include "../core/field_store.hpp"
#include "../core/grid.hpp"
#include "../core/system.hpp"
#include "../io/log.hpp"
#include "../ops/fd_ops.hpp"
#include "../physics/fe_ac_linear.hpp"
#include "../physics/fe_ch_landau.hpp"
#include "../physics/mobility.hpp"
#include "../terms/ac_term.hpp"
#include "../terms/ch_term.hpp"

#include <string_view>

namespace circa::cfg {

template <int D>
struct TermSpec {
    std::string id;
    std::string kind;        // "CH" | "AC" | ...
    std::string target;      // field to update
    std::string ops_type;    // "fd" | "spectral" | ...
    const toml::table* tbl;  // pointer into root TOML (kept alive by caller)
};

// -------------------- builder from TOML --------------------
//
// Expected TOML for each [[terms]]:
//
// [[terms]]
// id="ch_phi"; kind="CH"; target="phi"; enabled=true
//   [terms.ops] type="fd"
//   [terms.free_energy] type="landau"; eps=0.8; kappa=1.0
//   [terms.mobility] type="const"; M0=1.0
//
// [[terms]]
// id="ac_c"; kind="AC"; target="c"; enabled=true
//   [terms.ops] type="fd"
//   [terms.free_energy] type="linear"; Mc=0.1; gcoef=1.0
//   [terms.coupling] driver="phi"
//
// ------------------------------------------------------------

template <int D>
inline std::vector<TermSpec<D>> parse_term_specs(const toml::table& root) {
    std::vector<TermSpec<D>> specs;
    auto terms = root["terms"].as_array();
    if(!terms) {
        CIRCA_CRITICAL("[[terms]] missing or not an array");
        throw std::runtime_error("");
    }

    for(auto& node : *terms) {
        auto t = node.as_table();
        if(!t) continue;

        bool enabled = t->operator[]("enabled").template value<bool>().value_or(true);
        if(!enabled) {
            continue;
        }

        TermSpec<D> s;
        s.id = t->operator[]("id").template value<std::string>().value_or(std::string{});
        s.kind = t->operator[]("kind").template value<std::string>().value_or(std::string{});
        s.target = t->operator[]("target").template value<std::string>().value_or(std::string{});
        s.tbl = t;
        if(s.kind.empty() || s.target.empty())
            throw std::runtime_error("term missing 'kind' or 'target'");

        if(auto ops = t->operator[]("ops").as_table()) {
            s.ops_type = ops->operator[]("type").template value<std::string>().value_or("fd");
        }
        else {
            s.ops_type = "fd";
        }

        specs.push_back(std::move(s));
    }
    return specs;
}

// Concrete "ops" resolver (FD for now; extend later)
template <int D>
inline const DerivOps<D>& resolve_ops(const std::string& ops_type, const toml::table* /*ops_tbl*/) {
    if(ops_type == "fd") {
        static FDOps<D> fd;
        return fd;
    }
    throw std::runtime_error("Unknown ops.type: " + ops_type);
}

// -------------------- term builders --------------------

// CH term factory for specific FE & MOB types
template <int D, class FE, class MOB>
inline std::unique_ptr<ITerm<D>> make_CH_term(FieldStore<D>& S, FieldStore<D>& dS,
                                              const DerivOps<D>& ops_any,
                                              const std::string& target,
                                              const FE& fe, const MOB& mob, double k) {
    auto fd = dynamic_cast<const FDOps<D>*>(&ops_any);
    if(!fd) {
        throw std::runtime_error("CHTerm currently requires FDOps backend");
    }
    return std::make_unique<CHTerm<D, FE, MOB, FDOps<D>>>(S, dS, *fd, target, fe, mob, k);
}

// AC term factory for specific FE type
template <int D, class FE>
inline std::unique_ptr<ITerm<D>> make_AC_term(FieldStore<D>& S, FieldStore<D>& dS,
                                              const DerivOps<D>& ops_any,
                                              const std::string& target,
                                              const std::string& driver,
                                              const FE& fe) {
    auto fd = dynamic_cast<const FDOps<D>*>(&ops_any);
    if(!fd) {
        throw std::runtime_error("ACTerm currently requires FDOps backend");
    }
    return std::make_unique<ACTerm<D, FE, FDOps<D>>>(S, dS, *fd, target, driver, fe);
}

// Build one term instance from a TermSpec + its TOML subtree
template <int D>
inline std::unique_ptr<ITerm<D>> build_one_term(FieldStore<D>& S, FieldStore<D>& dS, const TermSpec<D>& spec) {
    // Resolve ops
    const toml::table *ops_tbl = as_table_ptr(spec.tbl->operator[]("ops"));
    const DerivOps<D>& ops = resolve_ops<D>(spec.ops_type, ops_tbl);

    // Branch on term kind
    if(spec.kind == "CH") {
        const toml::table* fe_tbl = as_table_ptr(spec.tbl->operator[]("free_energy"));
        const toml::table* mob_tbl = as_table_ptr(spec.tbl->operator[]("mobility"));
        if(!fe_tbl) {
            throw std::runtime_error(spec.id + ": [free_energy] missing");
        }
        // FE selection
        const std::string fe_type = value_or<std::string>(fe_tbl, "type", "");
        if(fe_type.empty()) {
            throw std::runtime_error(spec.id + ": free_energy.type missing");
        }

        // Mobility selection (default const)
        std::string mob_type = value_or<std::string>(mob_tbl, "type", "const");

        // --- FE_CH_Landau ---
        if(fe_type == "landau") {
            double k = *value_or_die<double>(*spec.tbl, "kappa");
            FE_CH_Landau fe(*fe_tbl);

            if(mob_type == "const") {
                MobConst<D> m;
                m.M0 = value_or<double>(mob_tbl, "M0", m.M0);
                return make_CH_term<D, FE_CH_Landau, MobConst<D>>(S, dS, ops, spec.target, fe, m, k);
            }
            else if(mob_type == "exp_of_field") {
                MobExpOfField<D> m;
                m.field = value_or<std::string>(mob_tbl, "field", "c");
                m.c0 = value_or<double>(mob_tbl, "c0", 1.0);
                return make_CH_term<D, FE_CH_Landau, MobExpOfField<D>>(S, dS, ops, spec.target, fe, m, k);
            }
            else {
                throw std::runtime_error(spec.id + ": unknown mobility.type: " + mob_type);
            }
        }

        throw std::runtime_error(spec.id + ": unknown CH free_energy.type: " + fe_type);
    } 
    else if(spec.kind == "AC") {
        const toml::table* fe_tbl = as_table_ptr(spec.tbl->operator[]("free_energy"));
        const toml::table* c_tbl = as_table_ptr(spec.tbl->operator[]("coupling"));
        if(!fe_tbl) throw std::runtime_error(spec.id + ": [free_energy] missing");

        std::string driver = value_or<std::string>(c_tbl, "driver", "phi");
        const std::string fe_type = value_or<std::string>(fe_tbl, "type", "");

        if(fe_type == "linear") {
            FE_AC_Linear fe;
            fe.Mc = value_or<double>(fe_tbl, "Mc", fe.Mc);
            fe.gcoef = value_or<double>(fe_tbl, "gcoef", fe.gcoef);
            return make_AC_term<D, FE_AC_Linear>(S, dS, ops, spec.target, driver, fe);
        }

        throw std::runtime_error(spec.id + ": unknown AC free_energy.type: " + fe_type);
    }

    throw std::runtime_error(spec.id + ": unknown term kind: " + spec.kind);
}

template <int D> GeneralConfig<D> load(const std::string& path) {
    GeneralConfig<D> config;

    config.raw_table = toml::parse_file(path.c_str());
    if(!config.raw_table) {
        throw std::runtime_error(fmt::format("Parsing failed with error '{}'", config.raw_table.error().description()));
	}

    config.seed = value_or<uint64_t>(config.raw_table, "seed", std::time(NULL));

    // grid
    auto gsec = config.raw_table["grid"];
    if(!gsec) {
        throw std::runtime_error("[grid] section missing");
    }

    if(gsec["n"].is_array()) {
        config.grid.n = array_from_toml<int, D>(*gsec["n"].as_array(), "grid.n");
    }
    else {
        config.grid.n.fill(*gsec["n"].template value<int>());
    }
    if(gsec["L"].is_array()) {
        config.grid.L = array_from_toml<double, D>(*gsec["L"].as_array(), "grid.L");
    }
    else {
        config.grid.L.fill(*gsec["L"].template value<double>());
    }

    // time
    if(auto t = config.raw_table["time"]) {
        config.time.dt = t["dt"].value_or(config.time.dt);
        config.time.steps = t["steps"].value_or(config.time.steps);
    }

    // output
    if(auto o = config.raw_table["output"]) {
        config.out.txt_append = o["txt_append"].value_or(config.out.txt_append);
        config.out.output_every = o["output_every"].value_or(config.out.output_every);
        config.out.conf_every = o["conf_every"].value_or(config.out.conf_every);
    }

    // integrator
    if(auto isec = config.raw_table["integrator"]) {
        config.integrator.name = isec["name"].value_or(config.integrator.name);
    }

    // fields
    auto fields = config.raw_table["fields"].as_array();
    if(!fields) {
        CIRCA_CRITICAL("[[fields]] missing or not an array");
        throw std::runtime_error("");
    }

    for(auto& node : *fields) {
        auto t = node.as_table();
        if(!t) {
            continue;
        }
        std::string name = *value_or_die<std::string>(*t, "name");
        std::string init_opt = *value_or_die<std::string>(*t, "initialisation");
        FieldInitialisation init_strat;
        if(init_opt == "constant") {
            init_strat.strategy = init_strat.CONSTANT;
            init_strat.average = *value_or_die<double>(*t, "average");
        }
        else if(init_opt == "random") {
            init_strat.strategy = init_strat.RANDOM;
            init_strat.average = *value_or_die<double>(*t, "average");
            init_strat.random_stddev = *value_or_die<double>(*t, "random_stddev");
        }
        else if(init_opt == "from_file") {
            init_strat.strategy = init_strat.READ_FROM_FILE;
            init_strat.filename = *value_or_die<std::string>(*t, "filename");
        }
        else {
            CIRCA_CRITICAL("Field '{}', the specified initialisation strategy '{}' is invalid", name, init_opt);
            throw std::runtime_error("");
        }

        config.fields.names.push_back(name);
        config.fields.init_strategies.push_back(init_strat);
    }

    auto specs = parse_term_specs<D>(config.raw_table);
    config.build_system_fn = [specs](FieldStore<D>& S_in, FieldStore<D>& dSdt_out) -> System<D> {
        System<D> sys;
        for(const auto& spec : specs) {
            auto term = build_one_term<D>(S_in, dSdt_out, spec);
            sys.add(std::move(term));
        }
        return sys;
    };

    return config;
}

template GeneralConfig<1> load(const std::string& path);
template GeneralConfig<2> load(const std::string& path);
template GeneralConfig<3> load(const std::string& path);

}  // namespace circa::cfg
