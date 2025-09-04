#pragma once

#include "../core/system.hpp"

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
};

struct OutputCfg {
    bool txt_append = false;
    int conf_every = 1000;
    int output_every = 100;
};

struct IntegratorCfg {
    std::string name = "euler";
};

struct FieldsCfg {
    std::vector<std::string> names;
};

template <int D>
struct GeneralConfig {
    toml::parse_result raw_table;

    GridCfg<D> grid{};
    TimeCfg time{};
    OutputCfg out{};
    IntegratorCfg integrator{};
    FieldsCfg fields{};
    BuildSysFn<D> build_system_fn;
};

template <int D> GeneralConfig<D> load(const std::string& path);

}  // namespace circa::cfg
