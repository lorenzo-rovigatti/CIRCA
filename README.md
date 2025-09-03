# CIRCA
## Coupled Integrators for Cahn–Hilliard / Allen–Cahn

CIRCA is a modular C++17 framework for phase-field PDEs with a focus on Cahn–Hilliard (CH) and Allen–Cahn (AC) equations, and their coupling. It is the heir to [an earlier attempt](https://github.com/lorenzo-rovigatti/cahn-hilliard).

Features
- Composable terms (free energies, mobilities)
- Swappable derivative operators (FD now; spectral later)
- Runtime-selectable integrators: Euler, RK2, RK4
- Clean class-based registry for integrators
- Minimal example (CPU/FD) builds out of the box with CMake

## Build
```
mkdir build && cd build
cmake ..
cmake --build . -j
./circa
```
