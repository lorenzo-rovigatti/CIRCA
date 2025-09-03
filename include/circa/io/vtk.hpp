#pragma once
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

#include "../core/field_store.hpp"
#include "../core/grid.hpp"

namespace circa::io {

namespace detail {
template <int D>
inline int dim_or_one(const std::array<int, D>& n, int idx) {
    return (idx < D) ? n[idx] : 1;
}
template <int D>
inline double dx_or_one(const std::array<double, D>& dx, int idx) {
    return (idx < D) ? dx[idx] : 1.0;
}
}  // namespace detail

// Write a single scalar field to VTK (STRUCTURED_POINTS, ASCII).
// Works for D=1/2/3; for D<3 we set missing dims to 1 and nz=1 (2D) or ny=nz=1 (1D).
template <int D>
inline void write_vtk_scalar(const Field<D>& f, const std::string& filename, const std::string& scalar_name) {
    const int nx = detail::dim_or_one(f.g.n, 0);
    const int ny = detail::dim_or_one(f.g.n, 1);
    const int nz = detail::dim_or_one(f.g.n, 2);
    const double sx = detail::dx_or_one(f.g.dx, 0);
    const double sy = detail::dx_or_one(f.g.dx, 1);
    const double sz = detail::dx_or_one(f.g.dx, 2);

    std::ofstream os(filename);
    if(!os) {
        CIRCA_CRITICAL("Cannot open {} for writing", filename);
        exit(1);
    }

    os << "# vtk DataFile Version 3.0\n";
    os << "CIRCA scalar output\n";
    os << "ASCII\n";
    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    os << "ORIGIN 0 0 0\n";
    os << "SPACING " << std::setprecision(16) << sx << " " << sy << " " << sz << "\n";
    os << "POINT_DATA " << (static_cast<size_t>(nx) * ny * nz) << "\n";
    os << "SCALARS " << scalar_name << " double\n";
    os << "LOOKUP_TABLE default\n";

    // VTK expects x fastest, then y, then z. Our flat() does x-fastest too,
    // so we can index with flat({i,j,k}).
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                std::array<int, 3> I3{i, j, k};
                // Build an index array of size D: missing dims are 0.
                std::array<int, D> ID{};
                if constexpr (D >= 1) ID[0] = I3[0];
                if constexpr (D >= 2) ID[1] = I3[1];
                if constexpr (D >= 3) ID[2] = I3[2];
                const int lin = flat<D>(ID, f.g.n);
                os << std::setprecision(16) << f.a[lin] << "\n";
            }
        }
    }
}

// Convenience: write all fields in a FieldStore to files like <dir>/<step>_<name>.vtk
template <int D>
inline void dump_all_fields_vtk(const FieldStore<D>& S, const std::string& out_dir, int step) {
    std::filesystem::create_directories(out_dir);
    for (const auto& kv : S.map) {
        const std::string& name = kv.first;
        const auto fname = out_dir + "/" + std::to_string(step) + "_" + name + ".vtk";
        write_vtk_scalar(kv.second, fname, name);
    }
}

}  // namespace circa::io
