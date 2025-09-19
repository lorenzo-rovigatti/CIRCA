#pragma once
#include <cassert>
#include <vector>

#include "../core/field.hpp"
#include "../core/field_store.hpp"

namespace circa {

struct FE_CH_MultiQuad {
    std::vector<double> a, b, kappa;       // size N
    std::vector<std::vector<double>> chi;  // NxN symmetric

    template <int D>
    std::vector<Field<D>> mu(const std::vector<const Field<D>*>& phi) const {
        const int N = (int)phi.size();

        assert((int)a.size() == N && (int)b.size() == N && (int)kappa.size() == N && (int)chi.size() == N);
        for(int i = 0; i < N; ++i) {
            assert((int)chi[i].size() == N);
        }

        std::vector<Field<D>> mu_values;
        mu_values.reserve(N);
        for(int i = 0; i < N; i++) {
            mu_values.emplace_back(phi[i]->g);
        }

        const int size = phi[0]->g.size;
        for(int p = 0; p < size; p++) {
            // cache φ_j(p)
            std::vector<double> ph(N);
            for(int j = 0; j < N; j++) {
                ph[j] = phi[j]->a[p];
            }

            for(int i = 0; i < N; i++) {
                double bulk = a[i] * ph[i] + b[i] * ph[i] * ph[i] * ph[i];  // a_i φ_i + b_i φ_i^3
                double coup = 0.0;
                for (int j = 0; j < N; j++) {
                    if (j != i) {
                        coup += chi[i][j] * ph[j];
                    }
                }
                mu_values[i].a[p] = bulk + coup;
            }
        }
        return mu_values;
    }

    template <int D>
    double bulk(const std::vector<const Field<D>*>& phi, int p) const {
        const int N = (int)phi.size();
        double s = 0.0;
        for(int i = 0; i < N; i++) {
            const double x = phi[i]->a[p];
            s += 0.5 * a[i] * x * x + 0.25 * b[i] * x * x * x * x;
        }
        for(int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; ++j) {
                s += chi[i][j] * phi[i]->a[p] * phi[j]->a[p];
            }
        }
        return s;
    }
};

}  // namespace circa
