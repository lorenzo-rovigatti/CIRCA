#pragma once
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <type_traits>

#include "../core/field_store.hpp"
#include "../core/grid.hpp"
#include "../util/strings.h"

namespace circa::io {

template <int D>
inline uint64_t init_field_from_plain(const std::string& filename, Field<D>& f) {
    std::ifstream is(filename);
    if(!is) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    // Read header line (starts with '#')
    if(!std::getline(is, line)) {
        throw std::runtime_error("File is empty: " + filename);
    }
    if(!util::starts_with(line, "#")) {
        throw std::runtime_error("Expected header starting with '#': " + filename);
    }

    util::trim(line);
    auto spl = util::split(line, ",");
    uint64_t initial_time;
    if(spl.size() == 3) {
        // parse the initial time
        auto time_spl = util::split(spl[0], "=");
        initial_time = std::stoll(time_spl[1]);
        CIRCA_INFO("Initial time step (as parsed from '{}'): {}", filename, initial_time);

        // parse the grid size
        auto size_spl = util::split(spl[2], "=");
        auto size_str = size_spl[1];
        // remove the x that can be used to separate dimensions
        for(auto& c : size_str) {
            if(c == 'x') c = ' ';
        }
        std::istringstream ss(size_str);
        std::vector<int> sizevals;
        int v;
        while (ss >> v) sizevals.push_back(v);
        if((int)sizevals.size() != D) {
            throw std::runtime_error(fmt::format("Dimension mismatch when reading '{}'", filename));
        }
        for(int d = 0; d < D; ++d) {
            if(sizevals[d] != f.g.n[d]) {
                throw std::runtime_error(fmt::format("Grid size mismatch in {}: size along the dimension {} is {}, should be {}", filename, d, sizevals[d], f.g.n[d]));
            }
        }
    }

    // Now read data
    if constexpr (D == 1) {
        for(int i = 0; i < f.g.n[0]; ++i) {
            if (!(is >> f.a[i])) {
                throw std::runtime_error("Unexpected EOF in " + filename);
            }
        }
    } 
    else if constexpr (D == 2) {
        for(int j = 0; j < f.g.n[1]; ++j) {
            for(int i = 0; i < f.g.n[0]; ++i) {
                int idx = j * f.g.n[0] + i;
                if(!(is >> f.a[idx])) {
                    throw std::runtime_error("Unexpected EOF in " + filename);
                }
            }
        }
    } 
    else {
        static_assert(D <= 2, "Plaintext format only defined for D=1,2");
    }

    return initial_time;
}

// Write a single scalar field to a plain-text file.
// Header format:
//   # step = XXX, t = XXX, size = Nx[,Ny], dx = dx[,dy]
// Data:
//   D=1: single column (Nx lines)
//   D=2: Ny rows, Nx columns (matrix)
// If append=true, header+data are appended to the file.
template <int D>
inline void write_field_to_plain(const Field<D>& f,
                               const std::string& filename,
                               int step, double t,
                               bool append = false) {
    if constexpr (DIM > 2) {
        return;
    }

    std::ofstream os(filename, append ? (std::ios::out | std::ios::app)
                                      : (std::ios::out | std::ios::trunc));
    if (!os) {
        CIRCA_CRITICAL("Cannot open {} for writing", filename);
        throw std::runtime_error("");
    }

    const int nx = f.g.n[0];
    const double dx = f.g.dx[0];

    os << std::setprecision(16);
    if constexpr (D == 1) {
        os << fmt::format("# step = {}, t = {}, size = {}, dx = {}", step, t, nx, dx) << std::endl;

        for(int i = 0; i < nx; ++i) {
            std::array<int, 1> I{i};
            const int lin = flat<1>(I, f.g.n);
            os << f.a[lin] << std::endl;
        }
    } 
    else if constexpr (D == 2) {
        const int ny = f.g.n[1];
        const double dy = f.g.dx[1];

        os << fmt::format("# step = {}, t = {}, size = {} {}, dx = {} {}", step, t, nx, ny, dx, dy) << std::endl;
        // Row-major print: y as rows, x as columns
        for(int j = 0; j < ny; ++j) {
            for(int i = 0; i < nx; ++i) {
                std::array<int, 2> I{i, j};
                const int lin = flat<2>(I, f.g.n);
                os << f.a[lin];
                if (i + 1 < nx) os << " ";
            }
            os << std::endl;
        }
    }
    os << std::endl;
}

template <int D>
inline void dump_all_fields_plain(const FieldStore<D>& S,
                                  const std::string& prefix,
                                  int step, double t,
                                  bool append = false) {
    if constexpr (DIM > 2) {
        return;
    }

    for(const auto& kv : S.map) {
        const std::string& name = kv.first;
        const auto& f = kv.second;

        std::string fname = prefix + "_" + name + ".dat";
        write_field_to_plain<D>(f, fname, step, t, append);
    }
}

}  // namespace circa::io
