# 3D Quasi-Magnetostatic (QMS) Field Solver

A high-performance 3D Quasi-Magnetostatic (QMS) field solver suite coupling a Fortran 2018 / Intel MKL PARDISO computational engine with Python domain setup, automated job control, Matplotlib spectral analysis, and ParaView VTK post-processing pipelines.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Governing Physics & Mathematical Formulation](#governing-physics--mathematical-formulation)
  - [Quasi-Magnetostatic Approximations](#quasi-magnetostatic-approximations)
  - [Coulomb Gauge Vector Potential ($\mathbf{A}-\phi$)](#coulomb-gauge-vector-potential-a-\phi)
  - [Discretization & Linear System Formulations](#discretization--linear-system-formulations)
- [Prerequisites & System Requirements](#prerequisites--system-requirements)
- [Compilation & Build Guide](#compilation--build-guide)
- [Simulation Workflow & Input Structure](#simulation-workflow--input-structure)
  - [1. Geometry & Domain Setup (`master.py`)](#1-geometry--domain-setup-masterpy)
  - [2. Input Specification (`inputs.txt`)](#2-input-specification-inputstxt)
  - [3. Binary Input Buffers](#3-binary-input-buffers)
- [Execution & Solver Output](#execution--solver-output)
- [Post-Processing & Visualization](#post-processing--visualization)
  - [Matplotlib Rendering & FFT Analysis (`plot.py`)](#matplotlib-rendering--fft-analysis-plotpy)
  - [ParaView Static Field Export (`paraview_still_frame.py`)](#paraview-static-field-export-paraview_still_framepy)
  - [ParaView Time-Series Visualization (`paraview_time.py`)](#paraview-time-series-visualization-paraview_timepy)
- [License & Citation](#license--citation)

---

## Overview

The **3D QMS Field Solver** is designed for high-resolution electromagnetic modeling where propagation delay effects (retardation) are negligible relative to inductive and diffusive mechanisms, such as in low-to-mid frequency magnetic devices, induction heating systems, magnetic sensors, and complex conductor structures.

Key features include:
* **Fortran 2018 Engine (`QMS_3D.f90`)**: Fully parallelized spatial stencil discretizations integrated with **Intel MKL PARDISO** direct sparse solvers.
* **Coulomb Gauge Formulation**: Enforces $
abla \cdot \mathbf{A} = 0$ to decouple potential fields and eliminate spurious modes.
* **Automated Python Pipeline (`master.py`)**: Direct generation of complex 3D material property tensors ($\mu_r$, $\sigma$) and source current distributions ($\mathbf{J}_{	ext{ext}}$).
* **Flexible Visualization Suite**:
  * Matplotlib scripts (`plot.py`) for spatial slice heatmaps, signal histories, and Fast Fourier Transform (FFT) frequency spectrum extraction.
  * ParaView integration scripts (`paraview_still_frame.py`, `paraview_time.py`) exporting native VTK ImageData (`.vti`) formats for volume rendering and field line streamlines.

---

## Repository Structure

| File / Module | Language / Tool | Description |
| :--- | :--- | :--- |
| `QMS_3D.f90` | Fortran 2018 | Core 3D Quasi-Magnetostatic solver engine using Intel MKL PARDISO. |
| `compile.sh` | Bash Script | Build script compiling `QMS_3D.f90` with `ifx` and MKL linking. |
| `master.py` | Python 3 | Master driver script for parameter definition, grid construction, binary generation, and Fortran execution. |
| `plot.py` | Python 3 | Post-processor for 2D field slice plotting, time-domain animations (MP4/GIF), and FFT spectral plots. |
| `paraview_still_frame.py` | Python 3 | VTK exporter for static single-frame snapshot datasets of magnetic fields ($\mathbf{B}, \mathbf{H}$). |
| `paraview_time.py` | Python 3 | Multi-frame VTK exporter for time-dependent transient magnetic fields. |
| `LICENSE.txt` | Text | MIT Open Source License (CU Boulder Regents & Daniel Richardson). |

---

## Governing Physics & Mathematical Formulation

### Quasi-Magnetostatic Approximations

In the quasi-magnetostatic (QMS) regime, the displacement current density $ rac{\partial \mathbf{D}}{\partial t}$ in Maxwell's equations is neglected ($ rac{\partial \mathbf{D}}{\partial t} \ll \mathbf{J}$):

$$
abla 	imes \mathbf{H} = \mathbf{J}$$

$$
abla 	imes \mathbf{E} = - rac{\partial \mathbf{B}}{\partial t}$$

$$
abla \cdot \mathbf{B} = 0$$

$$
abla \cdot \mathbf{D} = 
ho_v$$

Constitutive relations couple the macroscopic fields:

$$\mathbf{B} = \mu \mathbf{H} = \mu_0 \mu_r \mathbf{H}$$

$$\mathbf{J} = \sigma \mathbf{E} + \mathbf{J}_{	ext{ext}}$$

where $\mu_0$ is the vacuum permeability ($4\pi 	imes 10^{-7} 	ext{ H/m}$), $\mu_r(\mathbf{r})$ is the relative magnetic permeability tensor/scalar, $\sigma(\mathbf{r})$ is electrical conductivity ($	ext{S/m}$), and $\mathbf{J}_{	ext{ext}}$ represents externally imposed source excitation currents.

---

### Coulomb Gauge Vector Potential ($\mathbf{A}-\phi$)

Representing the magnetic flux density in terms of the magnetic vector potential $\mathbf{A}$:

$$\mathbf{B} = 
abla 	imes \mathbf{A}$$

Substituting into Faraday's Law yields:

$$
abla 	imes \left( \mathbf{E} +  rac{\partial \mathbf{A}}{\partial t} 
ight) = 0 \implies \mathbf{E} = -
abla \phi -  rac{\partial \mathbf{A}}{\partial t}$$

where $\phi$ is the scalar electric potential. Substituting into Amp√®re's Law with conductivity $\sigma$:

$$
abla 	imes \left(  rac{1}{\mu} 
abla 	imes \mathbf{A} 
ight) = \mathbf{J}_{	ext{ext}} - \sigma 
abla \phi - \sigma  rac{\partial \mathbf{A}}{\partial t}$$

Using the vector identity $
abla 	imes (
abla 	imes \mathbf{A}) = 
abla (
abla \cdot \mathbf{A}) - 
abla^2 \mathbf{A}$ and imposing the **Coulomb Gauge Condition**:

$$
abla \cdot \mathbf{A} = 0$$

For a region with piecewise uniform magnetic properties $\mu$, the vector potential governing PDE simplifies to:

$$
abla^2 \mathbf{A} - \mu \sigma  rac{\partial \mathbf{A}}{\partial t} = -\mu \mathbf{J}_{	ext{ext}} + \mu \sigma 
abla \phi$$

---

### Discretization & Linear System Formulations

The spatial domain $(N_x 	imes N_y 	imes N_z)$ is discretized on a structured 3D Cartesian mesh with grid spacings $\Delta x, \Delta y, \Delta z$. Second-order central finite differences discretize the spatial Laplacian $
abla^2 \mathbf{A}$:

$$
abla^2 A_k  ig|_{i,j,k}  pprox  rac{A_{k, i+1,j,k} - 2A_{k, i,j,k} + A_{k, i-1,j,k}}{\Delta x^2} +  rac{A_{k, i,j+1,k} - 2A_{k, i,j,k} + A_{k, i,j-1,k}}{\Delta y^2} +  rac{A_{k, i,j,k+1} - 2A_{k, i,j,k} + A_{k, i,j,k-1}}{\Delta z^2}$$

Time-stepping applies implicit Euler or Crank-Nicolson schemes to ensure unconditional stability in highly conductive domains ($\sigma \gg 0$).

The sparse system matrix $\mathbf{K} \mathbf{u} = \mathbf{f}$ is passed directly to **Intel MKL PARDISO**, utilizing parallelized direct factorizations ($L U$ or $L D L^T$) to efficiently solve multi-million degree-of-freedom field problems.

---

## Prerequisites & System Requirements

### Hardware Requirements
* **CPU**: x86_64 Processor supporting AVX2 or AVX-512 instruction sets. Multi-core execution is highly recommended for PARDISO threads.
* **RAM**: 16 GB minimum (32 GB+ recommended for spatial grids larger than $128 	imes 128 	imes 128$).

### Software Dependencies
* **Intel OneAPI Toolkits**:
  * Intel Fortran Compiler (`ifx` or `ifort`).
  * Intel Math Kernel Library (Intel MKL).
* **Python Environment (Python 3.8+)**:
  * `numpy`
  * `matplotlib`
  * `pyvista` or `vtk` (for VTK data conversion)
* **ParaView** (v5.8+ recommended for interactive 3D visualization).

---

## Compilation & Build Guide

The Fortran backend source code `QMS_3D.f90` is compiled using the provided shell script `compile.sh`.

### `compile.sh` Script Mechanics:
```bash
#!/bin/bash
# Compilation script for QMS 3D Fortran Solver using Intel Fortran Compiler and MKL

ifx -O3 -qopenmp -qmkl=parallel QMS_3D.f90 -o QMS_3D
```

### Build Instructions:
1. Ensure Intel OneAPI environment variables are sourced:
   ```bash
   source /opt/intel/oneapi/setvars.sh
   ```
2. Grant execution permissions and run `compile.sh`:
   ```bash
   chmod +x compile.sh
   ./compile.sh
   ```
3. Upon successful compilation, an executable binary named `QMS_3D` will be generated in the root directory.

---

## Simulation Workflow & Input Structure

The workflow follows a 3-step pipeline: setup via Python, binary execution via Fortran, and visualization via Python/ParaView.

```
‚îå‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îê       ‚îå‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îê       ‚îå‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îê
‚îÇ   master.py     ‚îÇ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ>‚îÇ inputs.txt      ‚îÇ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ>‚îÇ   QMS_3D        ‚îÇ
‚îÇ (Geometry/Grid) ‚îÇ       ‚îÇ Binary Datasets ‚îÇ       ‚îÇ (Fortran Solver)‚îÇ
‚îî‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îò       ‚îî‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îò       ‚îî‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚î¨‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îò
                                                             ‚îÇ
                                                             ‚ñº
‚îå‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îê       ‚îå‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îê       ‚îå‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îê
‚îÇ ParaView (.vti) ‚îÇ<‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÇ plot.py / VTK   ‚îÇ<‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÇ Raw Output      ‚îÇ
‚îÇ 3D Rendering    ‚îÇ       ‚îÇ Post-Processors ‚îÇ       ‚îÇ Field Binaries  ‚îÇ
‚îî‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îò       ‚îî‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îò       ‚îî‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îÄ‚îò
```

### 1. Geometry & Domain Setup (`master.py`)

`master.py` configures global simulation parameters, defines geometry via 3D matrix masks, generates physical coefficient arrays ($\mu_r, \sigma$), writes the parameter config `inputs.txt`, and generates raw unformatted binary inputs.

Key variables defined within `master.py`:
* **Grid dimensions**: `Nx`, `Ny`, `Nz`
* **Spatial steps**: `dx`, `dy`, `dz`
* **Time parameters**: `dt`, `Nt` (Number of timesteps)
* **Frequency / Excitation parameters**: `freq`, excitation current pulse parameters
* **Material Maps**: Permeability matrix `mu_r[Nx, Ny, Nz]`, Conductivity matrix `sigma[Nx, Ny, Nz]`, External Current density matrices `Jx`, `Jy`, `Jz`.

### 2. Input Specification (`inputs.txt`)

The parameter control file `inputs.txt` created by `master.py` contains basic numerical inputs read directly by `QMS_3D.f90`:

```text
128 128 128      ! Nx, Ny, Nz (Grid Dimensions)
0.001 0.001 0.001 ! dx, dy, dz (Grid Spacing in meters)
1e-6 500         ! dt (Time step in seconds), Nt (Total time steps)
1.0              ! Relaxation / Solver Convergence Criteria
```

### 3. Binary Input Buffers

For performance, 3D spatial field structures are exported by `master.py` as unformatted IEEE 754 floating-point binary buffers:
* `mu_r.bin`: Relative permeability distribution matrix.
* `sigma.bin`: Electrical conductivity map ($	ext{S/m}$).
* `Jx_ext.bin`, `Jy_ext.bin`, `Jz_ext.bin`: Time-independent or spatial envelope matrices of external excitation current density ($	ext{A/m}^2$).

---

## Execution & Solver Output

Execute the main simulation driver script:

```bash
python3 master.py
```

`master.py` performs the following steps:
1. Builds geometry and coefficient matrices.
2. Writes `inputs.txt` and `.bin` material files.
3. Invokes the compiled binary `./QMS_3D`.
4. The Fortran engine reads inputs, initializes MKL PARDISO matrices, performs time-stepping, and streams output binary files at pre-configured intervals:
   * `Bx_time.bin`, `By_time.bin`, `Bz_time.bin`: Magnetic flux density arrays across timesteps.
   * `Hx_time.bin`, `Hy_time.bin`, `Hz_time.bin`: Magnetic field intensity arrays.
   * `A_time.bin`: Computed 3D magnetic vector potential distributions.

---

## Post-Processing & Visualization

### Matplotlib Rendering & FFT Analysis (`plot.py`)

`plot.py` handles 2D slice visualization, animation creation, and temporal signal extraction.

* **Spatial Slice Heatmaps**: Plots cross-sectional scalar cuts ($x$-$y$, $y$-$z$, $x$-$z$ planes) of $|\mathbf{B}|$ or $|\mathbf{H}|$.
* **Time-Series Tracking**: Monitors magnetic field values at specified virtual probe locations $(x_p, y_p, z_p)$ over time.
* **FFT Frequency Spectrum**: Computes Discrete Fourier Transforms (DFT/FFT) of transient probe signals to identify spectral harmonics.
* **Animation Export**: Generates `.mp4` or `.gif` video files showing time-varying field distributions.

Run matplotlib post-processing:
```bash
python3 plot.py
```

---

### ParaView Static Field Export (`paraview_still_frame.py`)

To convert single frame snapshot binary output into VTK ImageData (`.vti`) for high-fidelity 3D volume rendering:

```bash
python3 paraview_still_frame.py
```

* Exports vector field quantities $\mathbf{B} = (B_x, B_y, B_z)$ and scalar field magnitude $|\mathbf{B}|$.
* Output file: `QMS_still_frame.vti` (Loadable directly in ParaView).

---

### ParaView Time-Series Visualization (`paraview_time.py`)

To render transient multi-timestep animations in ParaView:

```bash
python3 paraview_time.py
```

* Converts time-series binary output datasets into a sequence of VTK ImageData files (`QMS_frame_0000.vti`, `QMS_frame_0001.vti`, ...).
* Generates a ParaView metadata collection file (`QMS_transient.pvd`) referencing the series, enabling time-step playback, streamlines, vector glyph animations, and iso-surface extractions in ParaView.

---

## License & Citation

This project is licensed under the **MIT License**.

```text
Copyright (c) Regents of the University of Colorado & Daniel Richardson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
