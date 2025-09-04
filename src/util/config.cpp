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

// -------------------- helpers --------------------
inline const toml::table* as_table_ptr(const toml::node_view<const toml::node>& nv) {
    if(!nv) return nullptr;
    return nv.as_table();
}

template <typename T>
inline T value_or_t(const toml::table& t, std::string_view key, T def) {
    if(auto v = t[key].template value<T>()) return *v;
    return def;
}

template <typename T>
inline T value_or_t(const toml::table* tp, std::string_view key, T def) {
    if(!tp) return def;
    return value_or_t<T>(*tp, key, def);
}

template <int D, typename T>
inline std::array<T, D> array_or(const toml::array* a, const std::array<T, D>& def) {
    if(!a || a->size() != D) return def;
    std::array<T, D> out{};
    for(int i = 0; i < D; ++i) {
        out[i] = (*a)[size_t(i)].template value<T>().value_or(def[i]);
    }
    return out;
}

// -------------------- parsed term spec --------------------

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
    if(!terms) throw std::runtime_error("[[terms]] missing or not an array");

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
                                              const FE& fe, const MOB& mob) {
    auto fd = dynamic_cast<const FDOps<D>*>(&ops_any);
    if(!fd) {
        throw std::runtime_error("CHTerm currently requires FDOps backend");
    }
    return std::make_unique<CHTerm<D, FE, MOB, FDOps<D>>>(S, dS, *fd, target, fe, mob);
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
    const toml::table* ops_tbl = as_table_ptr(spec.tbl->operator[]("ops"));
    const DerivOps<D>& ops = resolve_ops<D>(spec.ops_type, ops_tbl);

    // Branch on term kind
    if(spec.kind == "CH") {
        const toml::table* fe_tbl = as_table_ptr(spec.tbl->operator[]("free_energy"));
        const toml::table* mob_tbl = as_table_ptr(spec.tbl->operator[]("mobility"));
        if(!fe_tbl) throw std::runtime_error(spec.id + ": [free_energy] missing");
        // FE selection
        const std::string fe_type = value_or_t<std::string>(fe_tbl, "type", "");
        if(fe_type.empty()) throw std::runtime_error(spec.id + ": free_energy.type missing");

        // Mobility selection (default const)
        std::string mob_type = value_or_t<std::string>(mob_tbl, "type", "const");

        // --- FE_CH_Landau ---
        if(fe_type == "landau") {
            FE_CH_Landau fe;
            fe.eps = value_or_t<double>(fe_tbl, "eps", fe.eps);
            fe.kappa = value_or_t<double>(fe_tbl, "kappa", fe.kappa);

            if(mob_type == "const") {
                MobConst<D> m;
                m.M0 = value_or_t<double>(mob_tbl, "M0", m.M0);
                return make_CH_term<D, FE_CH_Landau, MobConst<D>>(S, dS, ops, spec.target, fe, m);
            } 
            else if(mob_type == "exp_of_field") {
                MobExpOfField<D> m;
                m.field = value_or_t<std::string>(mob_tbl, "field", "c");
                m.c0 = value_or_t<double>(mob_tbl, "c0", 1.0);
                return make_CH_term<D, FE_CH_Landau, MobExpOfField<D>>(S, dS, ops, spec.target, fe, m);
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

        std::string driver = value_or_t<std::string>(c_tbl, "driver", "phi");
        const std::string fe_type = value_or_t<std::string>(fe_tbl, "type", "");

        if(fe_type == "linear") {
            FE_AC_Linear fe;
            fe.Mc = value_or_t<double>(fe_tbl, "Mc", fe.Mc);
            fe.gcoef = value_or_t<double>(fe_tbl, "gcoef", fe.gcoef);
            return make_AC_term<D, FE_AC_Linear>(S, dS, ops, spec.target, driver, fe);
        }

        throw std::runtime_error(spec.id + ": unknown AC free_energy.type: " + fe_type);
    }

    throw std::runtime_error(spec.id + ": unknown term kind: " + spec.kind);
}

template <typename T, size_t N>
static std::array<T, N> array_from_toml(const toml::array& a, const char* key) {
    if(a.size() != N) {
        throw std::runtime_error(std::string("Expected ") + key + " to have " + std::to_string(N) + " elements");
    }
    std::array<T, N> out{};
    for(size_t i = 0; i < N; ++i) {
        out[i] = a[i].value<T>().value_or(T{});
    }
    return out;
}

template <int D> GeneralConfig<D> load(const std::string& path) {
    GeneralConfig<D> config;

    config.raw_table = toml::parse_file(path.c_str());
    if(!config.raw_table) {
        throw std::runtime_error(fmt::format("Parsing failed with error '{}'", config.raw_table.error().description()));
	}

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
    if(auto fsec = config.raw_table["fields"]) {
        if(auto names = fsec["names"].as_array()) {
            config.fields.names.clear();
            for(auto& v : *names) {
                if(auto s = v.template value<std::string>()) {
                    config.fields.names.push_back(*s);
                }
            }
        }
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
