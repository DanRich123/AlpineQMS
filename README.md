# 3D Quasi-Magnetostatic (QMS) Field Solver

A high-performance 3D Quasi-Magnetostatic (QMS) field solver suite coupling a Fortran 2018 / Intel MKL PARDISO computational engine with Python domain setup, automated job control, Matplotlib spectral analysis, and ParaView VTK (`pyevtk`) post-processing pipelines.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Governing Physics & Mathematical Formulation](#governing-physics--mathematical-formulation)
  - [Quasi-Magnetostatic Approximations](#quasi-magnetostatic-approximations)
  - [Unreduced Vector Potential Equation & Gauge Penalty](#unreduced-vector-potential-equation--gauge-penalty)
  - [Steady-State Initialization & Advective Motion](#steady-state-initialization--advective-motion)
  - [Spatial Staggering (Yee Grid) & Discretization](#spatial-staggering-yee-grid--discretization)
  - [Time Integration](#time-integration)
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

---

## Overview

The **3D QMS Field Solver** is designed for high-resolution electromagnetic modeling where propagation delay effects (retardation) are negligible relative to inductive and diffusive mechanisms, such as in low-to-mid frequency magnetic devices, induction heating systems, magnetic sensors, and complex conductor structures.

Key features include:
* **Fortran 2018 Engine (`QMS_3D.f90`)**: Fully parallelized spatial stencil discretizations integrated with **Intel MKL PARDISO** direct sparse solvers.
* **Gauge Penalty Formulation**: Retains the unreduced formulation for $\mathbf{A}$ and incorporates a numerical gauge penalty term to maintain stability without analytical gauge reduction.
* **Coulomb Gauge Enforcement**: Coulomb gauge is enforced directly however via background field and gradient fields and Yee grid placement for the fields.
* **Yee Grid Spatial Staggering**: Enforces strict geometric duality by positioning vector potential $\mathbf{A}$ along edges, scalar potential $V$ at nodes, and magnetic fields $(\mathbf{B}, \mathbf{H})$ at cell faces.
* **Automated Python Pipeline (`master.py`)**: Direct generation of complex 3D material property tensors ($\mu_r$, $\sigma$).
* **Flexible Visualization Suite**:
  * Matplotlib scripts (`plot.py`) for spatial slice heatmaps, signal histories, and Fast Fourier Transform (FFT) frequency spectrum extraction.
  * ParaView integration scripts using `pyevtk` (`paraview_still_frame.py`, `paraview_time.py`) exporting native VTK ImageData (`.vti`) formats for volume rendering and field line streamlines.

---

## Repository Structure

| File / Module | Language / Tool | Description |
| :--- | :--- | :--- |
| `QMS_3D.f90` | Fortran 2018 | Core 3D Quasi-Magnetostatic solver engine using Intel MKL PARDISO. |
| `compile.sh` | Bash Script | Build script compiling `QMS_3D.f90` with `ifx` and MKL linking. |
| `master.py` | Python 3 | Master driver script for parameter definition, grid construction, binary generation, and Fortran execution. |
| `plot.py` | Python 3 | Post-processor for 2D field slice plotting, time-domain animations (MP4/GIF), and FFT spectral plots. |
| `paraview_still_frame.py` | Python 3 | VTK exporter using `pyevtk` for static single-frame snapshot datasets of magnetic fields ($\mathbf{B}, \mathbf{H}$). |
| `paraview_time.py` | Python 3 | Multi-frame VTK exporter using `pyevtk` for time-dependent transient magnetic fields. |
| `LICENSE.txt` | Text | MIT Open Source License (CU Boulder Regents & Daniel Richardson). |

---

## Governing Physics & Mathematical Formulation

### Quasi-Magnetostatic Approximations

In the quasi-magnetostatic (QMS) regime, the displacement current density $\frac{\partial \mathbf{D}}{\partial t}$ in Maxwell's equations is neglected ($\frac{\partial \mathbf{D}}{\partial t} \ll \mathbf{J}$):

$$\nabla \times \mathbf{H} = \mathbf{J}$$

$$\nabla \times \mathbf{E} = -\frac{\partial \mathbf{B}}{\partial t}$$

$$\nabla \cdot \mathbf{B} = 0$$

Constitutive relations couple the macroscopic fields:

$$\mathbf{B} = \mu \mathbf{H} = \mu_0 \mu_r \mathbf{H}$$

$$\mathbf{J} = \sigma \mathbf{E} + \mathbf{J}_{\text{ext}}$$

> **Note**: In the current version, external excitation current density is set to zero ($$\mathbf{J}_{\text{ext}} = 0$$). Future releases will introduce full support for external excitation sources.

---

### Unreduced Vector Potential Equation & Gauge Penalty

The electric field $\mathbf{E}$ is defined in terms of the scalar potential $V$ and vector potential $\mathbf{A}$:

$$\mathbf{E} = -\nabla V - \frac{\partial \mathbf{A}}{\partial t}$$

Substituting into Ampère's Law with conductivity $\sigma$:

$$\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A} \right) + \sigma \frac{\partial \mathbf{A}}{\partial t} + \sigma \nabla V = 0$$

Rather than setting $$\nabla \cdot \mathbf{A} = 0$$ explicitly to eliminate terms analytically, the solver retains the full vector operator form for $$\mathbf{A}$$ and introduces a **gauge penalty parameter** ($$\gamma \nabla (\nabla \cdot \mathbf{A})$$) to enforce gauge stability numerically:

$$\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A} \right) - \gamma \nabla (\nabla \cdot \mathbf{A}) + \sigma \frac{\partial \mathbf{A}}{\partial t} + \sigma \nabla V = 0$$

Additionally, the solver solves the continutity equation:

$$\nabla \cdot J = 0$$

$$\nabla \cdot (\sigma \frac{\partial \mathbf{A}}{\partial t} + \sigma \nabla V) = 0$$

---

### Steady-State Initialization & Advective Motion

Before transient time-stepping begins, the solver computes the steady-state equation ($$\frac{\partial \mathbf{A}}{\partial t} = 0$$). To account for moving media or velocity convection effects, the steady-state solve incorporates the standard advective term $$\mathbf{v} \times (\nabla \times \mathbf{A})$$:

$$\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A} \right) - \gamma \nabla (\nabla \cdot \mathbf{A}) - \sigma \mathbf{v} \times (\nabla \times \mathbf{A}) + \sigma \nabla V = 0$$

Additionally, the solver solves the continutity equation:

$$\nabla \cdot J = 0$$

$$\nabla \cdot (\sigma \nabla V - \sigma \mathbf{v} \times (\nabla \times \mathbf{A})) = 0$$

The solution to this steady-state system serves as the initial state ($\mathbf{A}^0$) for subsequent transient time integration.

---

### Spatial Staggering (Yee Grid) & Discretization

To satisfy electromagnetic dualities and guarantee numerical preservation of $\nabla \cdot \mathbf{B} = 0$, field variables are spatially assigned on a staggered **Yee Grid**:

* **Scalar Electric Potential ($V$)**: Defined at cell nodes / vertices $(i, j, k)$.
* **Magnetic Vector Potential ($\mathbf{A}$)**: Assigned along cell edges, matching primary electric field locations:
  * $A_x$ at $(i + 1/2, j, k)$
  * $A_y$ at $(i, j + 1/2, k)$
  * $A_z$ at $(i, j, k + 1/2)$
* **Magnetic Flux Density ($\mathbf{B}$) and Intensity ($\mathbf{H}$)**: Assigned at cell faces normal to their respective vector components:
  * $B_x, H_x$ at $(i, j + 1/2, k + 1/2)$
  * $B_y, H_y$ at $(i + 1/2, j, k + 1/2)$
  * $B_z, H_z$ at $(i + 1/2, j + 1/2, k)$

Spatial derivatives utilize discrete staggered curl-curl and gradient operators mapped directly between cell edges, nodes, and faces.

---

### Time Integration

Temporal integration is evaluated via a direct finite difference approximation for $\frac{\partial \mathbf{A}}{\partial t}$:

$$\frac{\partial \mathbf{A}}{\partial t} \approx \frac{\mathbf{A}^{n+1} - \mathbf{A}^n}{\Delta t}$$

This discrete time derivative is substituted into the governing vector potential equation to advance field distributions step-by-step in time.

---

## Prerequisites & System Requirements

### Hardware Requirements
* **CPU**: x86_64 Processor supporting AVX2 or AVX-512 instruction sets. Multi-core execution is highly recommended for PARDISO threads.
* **RAM**: 16 GB minimum (32 GB+ recommended for spatial grids larger than $50 \times 50 \times 50$).

### Software Dependencies
* **Intel OneAPI Toolkits**:
  * Intel Fortran Compiler (`ifx` or `ifort`).
  * Intel Math Kernel Library (Intel MKL).
* **Python Environment (Python 3.8+ recommended but not required)**:
  * `numpy`
  * `matplotlib`
  * `pyevtk` (Python EVTK library for writing VTK files directly)
* **ParaView** (v5.8+ recommended for interactive 3D visualization).

---

## Compilation & Build Guide

The Fortran backend source code `QMS_3D.f90` is compiled using the provided shell script `compile.sh`.

### `compile.sh` Script Mechanics:
```bash
#!/bin/bash
# Compilation script for QMS 3D Fortran Solver using Intel Fortran Compiler and MKL

ifx -O3 -qmkl QMS_3D.f90 -o QMS_3D
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
┌─────────────────┐       ┌─────────────────┐       ┌─────────────────┐
│   master.py     │──────>│ inputs.txt      │──────>│   QMS_3D        │
│ (Geometry/Grid) │       │ Binary Datasets │       │ (Fortran Solver)│
└─────────────────┘       └─────────────────┘       └────────┬────────┘
                                                             │
                                                             ▼
┌─────────────────┐       ┌─────────────────┐       ┌─────────────────┐
│ ParaView (.vti) │<──────│ plot.py / EVTK  │<──────│ Raw Output      │
│ 3D Rendering    │       │ Post-Processors │       │ Field Binaries  │
└─────────────────┘       └─────────────────┘       └─────────────────┘
```

### 1. Geometry & Domain Setup (`master.py`)

`master.py` configures global simulation parameters, defines geometry via a 3D matrix mask for material properties ($\mu_r, \sigma$), writes the parameter config `inputs.txt`, and generates a raw unformatted binary for geometry.

Key variables defined within `master.py`:
* **Grid dimensions**: `Nx`, `Ny`, `Nz`
* **Spatial steps**: `dx`, `dy`, `dz`
* **Time parameters**: `dt`, `Nt` (time step and number of time steps)
* **Material Maps**: Relative permeability matrix `mu_r[Nx, Ny, Nz]`, Conductivity matrix `sigma[Nx, Ny, Nz]`.

### 2. Input Specification (`inputs.txt`)

The parameter control file `inputs.txt` created by `master.py` contains basic numerical inputs read directly by `QMS_3D.f90`:

```text
128 128 128        ! Nx, Ny, Nz (Grid Dimensions)
0.001 0.001 0.001  ! dx, dy, dz (Grid Spacing in meters)
1e-6 500           ! dt (Time step in seconds), Nt (Total time steps)
```

### 3. Binary Input Buffers

For performance, a 3D spatial field structure is exported by `master.py` as a floating-point binary:
* `geom.bin`: Contains the material property ID number identifying the material

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
4. The Fortran engine reads inputs, initializes MKL PARDISO matrices, performs the steady-state solve, executes transient time-stepping, and streams output binary files at pre-configured intervals:
   * `Hx.bin`, `Hy.bin`, `Hz.bin`: Magnetic field intensity arrays at a 2D slice (cell centered) across time steps
   * `Hx_all.bin`, `Hy_all.bin`, `Hz_all.bin`: Magnetic field intensity arrays for all solution space (cell centered) across time steps
---

## Post-Processing & Visualization

### Matplotlib Rendering & FFT Analysis (`plot.py`)

`plot.py` handles 2D slice visualization, animation creation, and temporal signal extraction.

* **Spatial Slice Heatmaps**: Plots cross-sectional scalar cuts (y-z, x-z, or x-y planes) of $$|\mathbf{B}|$$ or $$|\mathbf{H}|$$.
* **Time-Series Tracking**: Monitors magnetic field values at specified virtual probe locations $$(x_p, y_p, z_p)$$ over time.
* **FFT Frequency Spectrum**: Computes Discrete Fourier Transforms (DFT/FFT) of transient probe signals to identify spectral harmonics.
* **Animation Export**: Generates `.mp4` or `.gif` video files showing time-varying field distributions.

Run matplotlib post-processing:
```bash
python3 plot.py
```

---

### ParaView Static Field Export (`paraview_still_frame.py`)

To convert single frame snapshot binary output into VTK ImageData (`.vti`) using `pyevtk` for high-fidelity 3D volume rendering:

```bash
python3 paraview_still_frame.py
```

* Exports vector field quantities $\mathbf{B} = (B_x, B_y, B_z)$ or $\mathbf{H} = (H_x, H_y, H_z)$.
* Output file: `QMS_still_frame.vti` (Loadable directly in ParaView - allows for viewing magnitude or components seperately).

---

### ParaView Time-Series Visualization (`paraview_time.py`)

To render transient multi-timestep animations in ParaView using `pyevtk`, very similar to the still frame version:

```bash
python3 paraview_time.py
```

* Converts time-series binary output datasets into a sequence of VTK ImageData files (`QMS_frame_0000.vti`, `QMS_frame_0001.vti`, ...).
* ParaView allows for referencing and opening the series of files as a single file, enabling time-step playback, streamlines, vector glyph animations, and iso-surface extractions in ParaView.

---

