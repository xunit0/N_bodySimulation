# N-Body Simulation in C++

This project is a simple **N-body gravitational simulator** written in C++.  
It computes the gravitational forces between particles (planets, moons, stars, or random bodies) and integrates their motion forward in time. Output is written in a `.tsv` format so it can be visualized with a provided Python plotting script.

---

## Features
- Particle state includes: **mass, position, velocity, force**.
- Built-in initialization modes:
  - `sem` → Sun–Earth–Moon toy system
  - Random bodies (user chooses `N`)
  - Load from `.tsv` file (same format as output)
- Gravitational force calculation with softening factor to prevent singularities.
- Euler integration of motion (`v += a*dt; x += v*dt`).
- Output in `.tsv` format for plotting with Python.
- Centering & scaling output (optional) so solar system fits nicely on charts.

---

## 🛠️ Build Instructions

### Prerequisites
- C++17 or newer (e.g., `g++`, `clang++`)
- `make` (for using the Makefile)
- Python 3 with `matplotlib` (for plotting)

### Build
Clone and build:
```bash
git clone https://github.com/xunit0/N-bodySimulation.git
cd N-bodySimulation
make run
