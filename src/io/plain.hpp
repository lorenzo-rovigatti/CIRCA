#pragma once
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <type_traits>

#include "../core/field_store.hpp"
#include "../core/grid.hpp"

namespace circa::io {

// Write a single scalar field to a plain-text file.
// Header format:
//   # step = XXX, t = XXX, size = Nx[,Ny], dx = dx[,dy]
// Data:
//   D=1: single column (Nx lines)
//   D=2: Ny rows, Nx columns (matrix)
// If append=true, header+data are appended to the file.
template <int D>
inline void write_plain_scalar(const Field<D>& f,
                               const std::string& filename,
                               int step, double t,
                               bool append = false) {
    if(DIM > 2) {
        return;
    }

    std::ofstream os(filename, append ? (std::ios::out | std::ios::app)
                                      : (std::ios::out | std::ios::trunc));
    if(!os) {
        CIRCA_CRITICAL("Cannot open {} for writing", filename);
        exit(1);
    }

    const int nx = f.g.n[0];
    const double dx = f.g.dx[0];

    os << std::setprecision(16);
    if constexpr (D == 1) {
        os << "# step = " << step
           << ", t = " << t
           << ", size = " << nx
           << ", dx = " << dx
           << "\n";
        for (int i = 0; i < nx; ++i) {
            std::array<int, 1> I{i};
            const int lin = flat<1>(I, f.g.n);
            os << f.a[lin] << "\n";
        }
    } 
    else if constexpr (D == 2) {
        const int ny = f.g.n[1];
        const double dy = f.g.dx[1];
        os << "# step = " << step
           << ", t = " << t
           << ", size = " << nx << " " << ny
           << ", dx = " << dx << " " << dy
           << "\n";
        // Row-major print: y as rows, x as columns
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                std::array<int, 2> I{i, j};
                const int lin = flat<2>(I, f.g.n);
                os << f.a[lin];
                if (i + 1 < nx) os << " ";
            }
            os << "\n";
        }
    }
    // optional blank line separator when appending multiple snapshots:
    os << "\n";
}

template <int D>
inline void dump_all_fields_plain(const FieldStore<D>& S,
                                  const std::string& prefix,
                                  int step, double t,
                                  bool append = false) {
    if(DIM > 2) {
        return;
    }

    for(const auto& kv : S.map) {
        const std::string& name = kv.first;
        const auto& f = kv.second;

        std::string fname = prefix + "_" + name + ".dat";
        write_plain_scalar<D>(f, fname, step, t, append);
    }
}

}  // namespace circa::io
