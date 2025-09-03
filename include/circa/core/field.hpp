#pragma once
#include <algorithm>
#include <vector>

#include "grid.hpp"

namespace circa {

template <int D>
struct Field {
    Grid<D> g;
    std::vector<double> a;
    Field() = default;
    explicit Field(const Grid<D>& gg) : g(gg), a(gg.size, 0.0) {}
    double& at(int i) { return a[i]; }
    double at(int i) const { return a[i]; }
    bool empty() const { return a.empty(); }
    void fill(double v) { std::fill(a.begin(), a.end(), v); }
};

template <int D>
inline double mean(const Field<D>& f) {
    double s = 0;
    for (double v : f.a) s += v;
    return f.g.size ? s / f.g.size : 0.0;
}
template <int D>
inline double var(const Field<D>& f) {
    double m = mean(f), s = 0;
    for (double v : f.a) {
        double d = v - m;
        s += d * d;
    }
    return f.g.size ? s / f.g.size : 0.0;
}

}  // namespace circa
