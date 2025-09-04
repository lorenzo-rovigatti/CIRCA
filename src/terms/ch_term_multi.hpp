#pragma once
#include <string>
#include <type_traits>
#include <vector>

#include "../core/system.hpp"
#include "../ops/deriv_ops.hpp"

namespace circa {

// ---- detection idiom for C++17 ----
template <typename T, int D>
struct has_M_i {
   private:
    template <typename U>
    static auto test(int) -> decltype(std::declval<U>().M_i(0, 0, std::declval<const FieldStore<D>&>()), std::true_type{});
    template <typename>
    static std::false_type test(...);

   public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

// MOB must offer either:
//   double M_i(int i_species, int idx_site, const FieldStore<D>& S)          // diagonal case
// or
//   double M_ibeta(int i_species, int beta_species, int idx_site, const FieldStore<D>& S) // full matrix

template <int D, class FE, class MOB, class Ops>
struct CHMultiTerm : ITerm<D> {
    FieldStore<D>* S = nullptr;
    FieldStore<D>* dSdt = nullptr;
    const Ops& ops;
    std::vector<std::string> target;  // names of the N species
    FE fe;
    MOB mob;

    CHMultiTerm(FieldStore<D>& S0, FieldStore<D>& dS0, const Ops& ops_,
                std::vector<std::string> targets, FE fe_, MOB mob_)
        : S(&S0), dSdt(&dS0), ops(ops_), target(std::move(targets)), fe(fe_), mob(mob_) {}

    void set_state(FieldStore<D>* Sin, FieldStore<D>* dSout) override {
        S = Sin;
        dSdt = dSout;
    }

    void add_rhs() override {
        const int N = (int)target.size();
        // gather φ_i and ∇²φ_i
        std::vector<const Field<D>*> phi(N);
        std::vector<Field<D>> lap(N);
        for (int a = 0; a < N; ++a) {
            phi[a] = &S->get(target[a]);
            lap[a] = ops.laplacian(*phi[a]);
        }

        // μ_i
        std::vector<Field<D>> mu = fe.template compute_mu<D>(*S, phi, lap);

        // ∇μ_i
        std::vector<std::array<Field<D>, D>> grad_mu;
        grad_mu.resize(N);
        for (int i = 0; i < N; ++i) {
            for (int d = 0; d < D; ++d) grad_mu[i][d] = Field<D>(phi[i]->g);
            grad_mu[i] = ops.gradient(mu[i]);
        }

        // For each species i: J_i = -sum_beta M_{iβ} ∇μ_β   (diagonal => only β=i)
        for (int i = 0; i < N; ++i) {
            std::array<Field<D>, D> flux;
            for (int d = 0; d < D; ++d) flux[d] = Field<D>(phi[i]->g);

            if constexpr (has_M_i<MOB, D>::value) {
                // diagonal mobility
                for (int p = 0; p < phi[i]->g.size; ++p) {
                    const double Mi = mob.M_i(i, p, *S);
                    for (int d = 0; d < D; ++d)
                        flux[d].a[p] = -Mi * grad_mu[i][d].a[p];
                }
            } else {
                // full matrix mobility
                for (int p = 0; p < phi[i]->g.size; ++p) {
                    for (int d = 0; d < D; ++d) {
                        double acc = 0.0;
                        for (int b = 0; b < N; ++b) {
                            const double Mib = mob.M_ibeta(i, b, p, *S);
                            acc += -Mib * grad_mu[b][d].a[p];
                        }
                        flux[d].a[p] = acc;
                    }
                }
            }

            // dφ_i/dt = -∇·J_i
            Field<D> dphi_dt = ops.divergence(flux);
            Field<D>& out = dSdt->ensure(target[i]);
            if (out.empty()) out = Field<D>(phi[i]->g);
            for (int p = 0; p < phi[i]->g.size; ++p) out.a[p] += dphi_dt.a[p];
        }
    }
};

}  // namespace circa
