#pragma once

#include <array>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#define TOML_EXCEPTIONS 0
#define TOML_ENABLE_FORMATTERS 0
// the two define's that follow are require to work around a known toml++ bug (see https://github.com/marzer/tomlplusplus/issues/213)
#define TOML_RETURN_BOOL_FROM_FOR_EACH_BROKEN 1
#define TOML_RETURN_BOOL_FROM_FOR_EACH_BROKEN_ACKNOWLEDGED 1
#include <toml++/toml.hpp>

namespace circa::cfg {

// ---------- user-facing config structs ----------
template <int D>
struct GridCfg {
    std::array<int, D> n{};
    std::array<double, D> L{};
};

struct TimeCfg {
    double dt = 1e-3;
    int steps = 1000;
    int out_every = 100;
};

struct OutputCfg {
    std::string vtk_dir = "out_vtk";
    std::string txt_dir = "out_txt";
    bool txt_append = false;
    int conf_every = 1000;
    int output_every = 100;
};

struct IntegratorCfg {
    std::string name = "rk2";
};

struct FieldsCfg {
    std::vector<std::string> names;
};

struct FE_CH_Landau_Cfg {
    double eps = 0.8, kappa = 1.0;
};

struct MobConst_Cfg {
    double M0 = 1.0;
};

struct FE_AC_Linear_Cfg {
    double Mc = 0.1, gcoef = 1.0;
};

template <int D>
struct PhysicsCfg {
    // CH term
    std::string ch_target = "phi";
    FE_CH_Landau_Cfg ch_fe{};
    MobConst_Cfg ch_mob{};
    // AC term
    std::string ac_field = "c";
    std::string ac_driver = "phi";
    FE_AC_Linear_Cfg ac_fe{};
};

template <int D>
struct GeneralConfig {
    toml::parse_result raw_table;

    GridCfg<D> grid{};
    TimeCfg time{};
    OutputCfg out{};
    IntegratorCfg integrator{};
    FieldsCfg fields{};
    PhysicsCfg<D> phys{};
};

// ---------- helpers ----------
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

// ---------- loader ----------
template <int D>
inline GeneralConfig<D> load(const std::string& path) {
    GeneralConfig<D> config;

    config.raw_table = toml::parse_file(path.c_str());
    if(!config.raw_table) {
        throw std::runtime_error(fmt::format("Parsing failed with error '{}'", config.raw_table.error().description()));
	}

    // grid
    auto gsec = config.raw_table["grid"];
    if(!gsec) throw std::runtime_error("[grid] section missing");
    config.grid.n = array_from_toml<int, D>(*gsec["n"].as_array(), "grid.n");
    config.grid.L = array_from_toml<double, D>(*gsec["L"].as_array(), "grid.L");

    // time
    if(auto t = config.raw_table["time"]) {
        config.time.dt = t["dt"].value_or(config.time.dt);
        config.time.steps = t["steps"].value_or(config.time.steps);
        config.time.out_every = t["out_every"].value_or(config.time.out_every);
    }

    // output
    if(auto o = config.raw_table["output"]) {
        config.out.vtk_dir = o["vtk_dir"].value_or(config.out.vtk_dir);
        config.out.txt_dir = o["txt_dir"].value_or(config.out.txt_dir);
        config.out.txt_append = o["txt_append"].value_or(config.out.txt_append);
        config.out.output_every = o["output_every"].value_or(config.out.output_every);
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

    // physics - CH
    if(auto p = config.raw_table["physics"]) {
        if(auto ch = p["ch"]) {
            config.phys.ch_target = ch["target"].value_or(config.phys.ch_target);
            if(auto fe = ch["fe"]) {
                config.phys.ch_fe.eps = fe["eps"].value_or(config.phys.ch_fe.eps);
                config.phys.ch_fe.kappa = fe["kappa"].value_or(config.phys.ch_fe.kappa);
            }
            if(auto mob = ch["mobility"]) {
                // type is present for future variants; we ignore if not "const"
                config.phys.ch_mob.M0 = mob["M0"].value_or(config.phys.ch_mob.M0);
            }
        }
        // physics - AC
        if(auto ac = p["ac"]) {
            config.phys.ac_field = ac["field"].value_or(config.phys.ac_field);
            config.phys.ac_driver = ac["driver"].value_or(config.phys.ac_driver);
            if(auto fe = ac["fe"]) {
                config.phys.ac_fe.Mc = fe["Mc"].value_or(config.phys.ac_fe.Mc);
                config.phys.ac_fe.gcoef = fe["gcoef"].value_or(config.phys.ac_fe.gcoef);
            }
        }
    }

    return config;
}

}  // namespace circa::cfg
