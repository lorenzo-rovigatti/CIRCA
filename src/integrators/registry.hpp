#pragma once
#include <memory>
#include <string>
#include <unordered_map>

#include "../util/config.hpp"
#include "euler.hpp"
#include "rk2.hpp"
#include "rk4.hpp"

namespace circa {

template<int D>
using IntegratorFactory = std::function<std::unique_ptr<IIntegrator<D>>(const cfg::GeneralConfig<D>&, const BuildSysFn<D>&, FieldStore<D>& /*S0*/)>;

template<int D>
std::unordered_map<std::string, IntegratorFactory<D>> make_integrator_registry() {
    std::unordered_map<std::string, IntegratorFactory<D>> R;

    R["euler"] = [](const cfg::GeneralConfig<D>& cfg, const BuildSysFn<D>& build, FieldStore<D>& S0){
        return std::make_unique<Euler<D>>(build, S0, cfg);
    };

    R["rk2"] = [](const cfg::GeneralConfig<D>& cfg, const BuildSysFn<D>& build, FieldStore<D>& S0){
        return std::make_unique<RK2<D>>(build, S0, cfg);
    };

    R["rk4"] = [](const cfg::GeneralConfig<D>& cfg, const BuildSysFn<D>& build, FieldStore<D>& S0) {
        return std::make_unique<RK4<D>>(build, S0, cfg);
    };

    return R;
}

}  // namespace circa
