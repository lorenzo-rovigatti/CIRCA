#pragma once
#include <string>
#include <unordered_map>

#include "field.hpp"

namespace circa {

template <int D>
struct FieldStore {
    Grid<D> g;
    std::unordered_map<std::string, Field<D>> map;
    explicit FieldStore(const Grid<D>& gg) : g(gg) {}

    Field<D>& ensure(const std::string& name) {
        auto it = map.find(name);
        if(it == map.end()) {
            it = map.emplace(name, Field<D>(g)).first;
        }
        return it->second;
    }

    const Field<D>& get(const std::string& name) const {
        auto it = map.find(name);
        if(it == map.end()) throw std::runtime_error("Missing field: " + name);
        return it->second;
    }

    const Field<D>* maybe(const std::string& name) const {
        auto it = map.find(name);
        return it == map.end() ? nullptr : &it->second;
    }
    
    void zero() {
        for(auto& kv : map) {
            kv.second.fill(0.0);
        }
    }
};

template <int D>
inline void axpy(Field<D>& y, const Field<D>& x, double a) {
    for (int i = 0; i < y.g.size; ++i) y.a[i] += a * x.a[i];
}
template <int D>
inline void axpy(FieldStore<D>& y, const FieldStore<D>& x, double a) {
    for (const auto& kv : x.map) {
        auto& yf = y.ensure(kv.first);
        if (yf.empty()) yf = Field<D>(y.g);
        axpy(yf, kv.second, a);
    }
}
template <int D>
inline FieldStore<D> plus_scaled(const FieldStore<D>& X, const FieldStore<D>& Y, double aX, double aY) {
    FieldStore<D> Z(X.g);
    for (const auto& kv : X.map) Z.ensure(kv.first);
    for (const auto& kv : Y.map) Z.ensure(kv.first);
    for (auto& kv : Z.map) {
        const Field<D>* xf = X.maybe(kv.first);
        const Field<D>* yf = Y.maybe(kv.first);
        for (int i = 0; i < Z.g.size; ++i) {
            double xv = xf ? xf->a[i] : 0.0;
            double yv = yf ? yf->a[i] : 0.0;
            kv.second.a[i] = aX * xv + aY * yv;
        }
    }
    return Z;
}

}  // namespace circa
