#pragma once
#include <array>

#include "../core/field.hpp"

namespace circa {

template <int D>
struct DerivOps {
    virtual ~DerivOps() = default;
    virtual Field<D> laplacian(const Field<D>& f) const = 0;
    virtual std::array<Field<D>, D> gradient(const Field<D>& f) const = 0;
    virtual Field<D> divergence(const std::array<Field<D>, D>& v) const = 0;
};

}  // namespace circa
