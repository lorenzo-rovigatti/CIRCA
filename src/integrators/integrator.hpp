#pragma once
#include <functional>

#include "../core/system.hpp"

namespace circa {

template <int D>
using BuildFn = std::function<System<D>(FieldStore<D>& /*S_in*/, FieldStore<D>& /*dSdt_out*/)>;

template <int D>
struct IIntegrator {
    System<D> sys_;

    IIntegrator(const BuildFn<D>& build, FieldStore<D>& S0) {
        // temporary dummy dSdt to satisfy constructor signatures
        FieldStore<D> dummy(S0.g);
        sys_ = build(S0, dummy); // moves terms in; they point to (&S0,&dummy) for now
    }

    virtual ~IIntegrator() = default;
    virtual void step(FieldStore<D>& S, double dt) = 0;
};

}  // namespace circa
