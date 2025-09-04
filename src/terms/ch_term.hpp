#pragma once
#include "../core/system.hpp"
#include "../ops/deriv_ops.hpp"

namespace circa {

template <int D, class FE, class M, class Ops>
struct CHTerm : ITerm<D>, IEnergy<D> {
    FieldStore<D>* S = nullptr;
    FieldStore<D>* dSdt = nullptr;
    const Ops& ops;
    std::string target;
    FE fe;
    M Mfun;
    double kappa;

    CHTerm(FieldStore<D>& S0, FieldStore<D>& dS0, const Ops& ops_, std::string tgt, FE fe_, M m_, double k)
        : S(&S0), dSdt(&dS0), ops(ops_), target(std::move(tgt)), fe(fe_), Mfun(m_), kappa(k) {}

    void set_state(FieldStore<D>* Sin, FieldStore<D>* dSout) override {
        S = Sin;
        dSdt = dSout;
    }

    void add_rhs() override {
        const Field<D>& u = S->get(target);
        Field<D> lap_u = ops.laplacian(u);

        Field<D> mu(u.g);
        for(int i = 0; i < u.g.size; ++i) {
            mu.a[i] = fe.mu(u.a[i]) - kappa * lap_u.a[i];
        }

        auto grad_mu = ops.gradient(mu);
        std::array<Field<D>, D> flux{Field<D>(u.g)};
        for(int d = 1; d < D; ++d) {
            flux[d] = Field<D>(u.g);
        }
        
        for(int i = 0; i < u.g.size; ++i) {
            double Mv = Mfun(i, *S);
            for (int d = 0; d < D; ++d) flux[d].a[i] = Mv * grad_mu[d].a[i];
        }

        Field<D> dudt = ops.divergence(flux);
        Field<D>& out = dSdt->ensure(target);
        if(out.empty()) out = Field<D>(u.g);
        for(int i = 0; i < u.g.size; ++i) {
            out.a[i] += dudt.a[i];
        }
    }

    double energy() const override {
        const Field<D>& u = S->get(target);

        // gradient
        auto gu = ops.gradient(u);

        double E = 0.0;
        for (int i = 0; i < u.g.size; ++i) {
            double grad2 = 0.0;
            for (int d = 0; d < D; ++d) {
                grad2 += gu[d].a[i] * gu[d].a[i];
            }
            double e_bulk = fe.bulk(u.a[i]);
            double e_grad = 0.5 * kappa * grad2;
            E += (e_bulk + e_grad) * u.g.dV;
        }
        return E;
    }
};

}  // namespace circa
