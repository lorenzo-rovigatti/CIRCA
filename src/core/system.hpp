#pragma once
#include <memory>
#include <vector>

#include "../core/field_store.hpp"

namespace circa {

template <int D>
struct ITerm {
    virtual ~ITerm() = default;
    virtual void add_rhs() = 0;
    virtual void set_state(FieldStore<D>* S_in, FieldStore<D>* dSdt_out) = 0;
};

template <int D>
struct IEnergy {
    virtual ~IEnergy() = default;
    // Return total free-energy contribution (integrated over space)
    virtual double energy() const = 0;
};

template <int D>
struct System {
    std::vector<std::unique_ptr<ITerm<D>>> terms;
    void add(std::unique_ptr<ITerm<D>> t) { terms.emplace_back(std::move(t)); }
    void rhs() {
        for (auto& t : terms) t->add_rhs();
    }
    void set_state(FieldStore<D>* S_in, FieldStore<D>* dSdt_out) {
        for (auto& t : terms) t->set_state(S_in, dSdt_out);
    }
};

template <int D>
using BuildSysFn = std::function<System<D>(FieldStore<D>& /*S_in*/, FieldStore<D>& /*dSdt_out*/)>;

}  // namespace circa
