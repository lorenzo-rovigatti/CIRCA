#pragma once
#include <memory>
#include <string>
#include <unordered_map>

#include "euler.hpp"
#include "integrator.hpp"
#include "rk2.hpp"
#include "rk4.hpp"

namespace circa {

struct ParsedInput {
    std::string integrator = "rk2";  // "euler" | "rk2" | "rk4"
    bool mass_fix = false;
    std::string mass_field = "phi";
};

template<int D>
using IntegratorFactory = std::function<std::unique_ptr<IIntegrator<D>>(const ParsedInput&, const BuildSysFn<D>&, FieldStore<D>& /*S0*/)>;

template<int D>
std::unordered_map<std::string, IntegratorFactory<D>> make_integrator_registry() {
    std::unordered_map<std::string, IntegratorFactory<D>> R;

    R["euler"] = [](const ParsedInput& cfg, const BuildSysFn<D>& build, FieldStore<D>& S0){
        EulerOptions opt; opt.conservative_mass_fix = cfg.mass_fix;
        opt.mass_field = cfg.mass_field;
        return std::make_unique<Euler<D>>(build, S0, opt);
    };

    R["rk2"] = [](const ParsedInput& cfg, const BuildSysFn<D>& build, FieldStore<D>& S0){
        RK2Options opt; opt.conservative_mass_fix = cfg.mass_fix;
        opt.mass_field = cfg.mass_field;
        return std::make_unique<RK2<D>>(build, S0, opt);
    };

    R["rk4"] = [](const ParsedInput& cfg, const BuildSysFn<D>& build, FieldStore<D>& S0) {
        RK4Options opt;
        return std::make_unique<RK4<D>>(build, S0, opt);
    };

    return R;
}

}  // namespace circa
