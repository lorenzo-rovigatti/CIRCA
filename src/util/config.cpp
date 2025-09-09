#include "config.hpp"

#include "../core/field_store.hpp"
#include "../core/grid.hpp"
#include "../core/system.hpp"
#include "../io/log.hpp"
#include "../ops/fd_ops.hpp"
#include "../physics/fe_ac_gel.hpp"
#include "../physics/fe_ch_landau.hpp"
#include "../physics/fe_ch_wertheim.hpp"
#include "../physics/mobility.hpp"
#include "../terms/ac_term.hpp"
#include "../terms/ch_term.hpp"

#include <variant>
#include <string_view>

namespace circa::cfg {

using FE_CH_Any = std::variant<
    FE_CH_Landau,
    FE_CH_Wertheim
>;

FE_CH_Any parse_ch_fe_any(const toml::table& fe_tbl){
    const auto type = fe_tbl["type"].template value<std::string>().value_or("");
    if(type == "landau"){
        FE_CH_Landau fe(fe_tbl);
        return fe;
    }
    else if(type == "wertheim"){ 
        FE_CH_Wertheim fe(fe_tbl);
        return fe;
    }
    throw std::runtime_error("unknown CH free_energy.type: " + type);
}

template<int D>
using MobAny = std::variant<
    MobConst<D>,
    MobExpOfField<D>
>;

template<int D>
MobAny<D> parse_mob_any(const toml::table* mob_tbl){
    const std::string type = value_or<std::string>(mob_tbl, "type", "const");
    if(type == "const"){
        MobConst<D> m;
        m.M0 = value_or<double>(mob_tbl, "M0", m.M0);
        return m;
    }

    if(type == "exp_of_field"){
        MobExpOfField<D> m;
        m.field = value_or<std::string>(mob_tbl, "field", "c");
        m.c0    = value_or<double>(mob_tbl, "c0", 1.0);
        return m;
    }
    throw std::runtime_error("unknown mobility.type: " + type);
}

using FE_AC_Any = std::variant<
    FE_AC_Gel
>;

FE_AC_Any parse_ac_fe_any(const toml::table& fe_tbl){
    const std::string fe_type = value_or<std::string>(fe_tbl, "type", "");
    if(fe_type == "gel") {
        FE_AC_Gel fe(fe_tbl);
        return fe;
    }
    throw std::runtime_error("unknown AC free_energy.type: " + fe_type);
}

template <int D>
struct TermSpec {
    std::string id;
    std::string kind;        // "CH" | "AC" | ...
    std::string target;      // field to update
    std::string ops_type;    // "fd" | "spectral" | ...
    const toml::table* tbl;  // pointer into root TOML (kept alive by caller)
};

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
    if (spec.kind == "CH") {
        const toml::table* fe_tbl  = as_table_ptr(spec.tbl->operator[]("free_energy"));
        const toml::table* mob_tbl = as_table_ptr(spec.tbl->operator[]("mobility"));
        if(!fe_tbl) {
            throw std::runtime_error(spec.id + ": [free_energy] missing");
        }

        double k = *value_or_die<double>(spec.tbl, "kappa");

        auto fe_any  = parse_ch_fe_any(*fe_tbl);
        auto mob_any = parse_mob_any<D>(mob_tbl);

        return std::visit(
            [&](auto&& fe, auto&& mob) -> std::unique_ptr<ITerm<D>> {
                using FE  = std::decay_t<decltype(fe)>;
                using MOB = std::decay_t<decltype(mob)>;
                auto fd = dynamic_cast<const FDOps<D>*>(&ops);
                if(!fd) {
                    throw std::runtime_error("CHTerm requires FDOps backend");
                }
                return std::make_unique<CHTerm<D, FE, MOB, FDOps<D>>>(S, dS, *fd, spec.target, fe, mob, k);
            },
            fe_any, mob_any
        );
    }
    else if(spec.kind == "AC") {
        const toml::table *fe_tbl = as_table_ptr(spec.tbl->operator[]("free_energy"));
        const toml::table *c_tbl = as_table_ptr(spec.tbl->operator[]("coupling"));
        if(!fe_tbl) {
            throw std::runtime_error(spec.id + ": [free_energy] missing");
        }

        auto fe_any = parse_ac_fe_any(*fe_tbl);
        std::string driver = value_or<std::string>(c_tbl, "driver", "phi");

        return std::visit(
            [&](auto&& fe) -> std::unique_ptr<ITerm<D>> {
                using FE = std::decay_t<decltype(fe)>;
                auto fd = dynamic_cast<const FDOps<D>*>(&ops);
                if(!fd) {
                    throw std::runtime_error("ACTerm requires FDOps backend");
                }
                return std::make_unique<ACTerm<D, FE, FDOps<D>>>(S, dS, *fd, spec.target, driver, fe);
            },
            fe_any
        );
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

    // time
    if(auto t = config.raw_table["time"]) {
        config.time.dt = t["dt"].value_or(config.time.dt);
        config.time.steps = t["steps"].value_or(config.time.steps);
    }

    // output
    if(auto o = config.raw_table["output"]) {
        config.out.output_append = o["output_append"].value_or(config.out.output_append);
        config.out.output_filename = o["output_filename"].value_or(config.out.output_filename);
        config.out.output_every = *value_or_die<int>(*o.as_table(), "output_every");
        config.out.conf_every = *value_or_die<int>(*o.as_table(), "conf_every");

        if(auto arr = o["mass_fields"].as_array()) {
            config.out.mass_fields.reserve(arr->size());
            for(auto&& v : *arr) {
                if(auto s = v.template value<std::string>()) {
                    config.out.mass_fields.push_back(*s);
                }
                else {
                    throw std::runtime_error("mass_fields must be strings");
                }
            }
        }
        else {
            auto field_str = *value_or_die<std::string>(*o.as_table(), "mass_fields");
            config.out.mass_fields.push_back(field_str);
        }

        // check that the "mass fields" exist
        for(auto &s : config.out.mass_fields) {
            if(std::find(config.fields.names.begin(), config.fields.names.end(), s) == config.fields.names.end()) {
                throw std::runtime_error("mass_fields refers to unknown field: " + s);
            }
        }
    }

    // integrator
    if(auto isec = config.raw_table["integrator"]) {
        config.integrator.name = isec["name"].value_or(config.integrator.name);
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
