#pragma once

#include "../core/system.hpp"
#include "toml.hpp"

#include <array>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

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
};

struct OutputCfg {
    bool output_append = false;
    std::string output_filename = "energy.dat";
    int output_every;
    int conf_every;
};

struct IntegratorCfg {
    std::string name = "euler";
};

struct FieldInitialisation {
    enum Strategy {
        READ_FROM_FILE,
        RANDOM,
        CONSTANT
    };
    Strategy strategy;

    std::string filename;
    double average;
    double random_stddev;
};

struct FieldsCfg {
    std::vector<std::string> names;
    std::vector<FieldInitialisation> init_strategies;
};

template <int D>
struct GeneralConfig {
    toml::parse_result raw_table;

    uint64_t seed;
    GridCfg<D> grid{};
    TimeCfg time{};
    OutputCfg out{};
    IntegratorCfg integrator{};
    FieldsCfg fields{};
    BuildSysFn<D> build_system_fn;
};

template <int D> GeneralConfig<D> load(const std::string& path);

}  // namespace circa::cfg
