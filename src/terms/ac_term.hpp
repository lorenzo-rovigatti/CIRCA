#pragma once
#include "../core/system.hpp"

namespace circa {

template <int D, class FE, class Ops>
struct ACTerm : ITerm<D> {
    FieldStore<D>* S = nullptr;
    FieldStore<D>* dSdt = nullptr;
    const Ops& ops;
    std::string c_name, driver_name;
    FE fe;

    ACTerm(FieldStore<D>& S0, FieldStore<D>& dS0, const Ops& ops_,
           std::string cfield, std::string driver, FE fe_)
        : S(&S0), dSdt(&dS0), ops(ops_), c_name(std::move(cfield)), driver_name(std::move(driver)), fe(fe_) {}

    void set_state(FieldStore<D>* Sin, FieldStore<D>* dSout) override {
        S = Sin;
        dSdt = dSout;
    }

    void add_rhs() override {
        const Field<D>& c = S->get(c_name);
        const Field<D>* drv = driver_name.empty() ? nullptr : S->maybe(driver_name);
        Field<D>& out = dSdt->ensure(c_name);
        if (out.empty()) out = Field<D>(c.g);
        for (int i = 0; i < c.g.size; ++i) {
            double driver = drv ? drv->a[i] : 0.0;
            out.a[i] += -fe.dfdc(c.a[i], driver);
        }
    }
};

}  // namespace circa
