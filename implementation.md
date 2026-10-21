# Implementation Roadmap: Interactive Kronig–Penney Model & Wavefunction Explorer

This document outlines the complete, phased roadmap to build an interactive Python application for exploring the **Kronig–Penney Model** of an electron in a 1D periodic crystalline lattice.

Each phase is designed to be completed incrementally so you can understand the underlying quantum mechanics, the numerical methods, and the software architecture side-by-side.

---

## Pedagogical Goals
By the end of this project, you will:
1. Understand why electrons in a periodic crystal form **energy bands** and **band gaps** instead of discrete energy levels or a continuous spectrum.
2. See firsthand how **Bloch's Theorem** and boundary conditions shape the electron wavefunction $\psi(x)$ and probability density $|\psi(x)|^2$.
3. Learn how to solve the 1D **Time-Independent Schrödinger Equation (TISE)** both analytically (secular equation/transfer matrix) and numerically (finite difference / shooting method).
4. Master building scientific Python interfaces with real-time reactive plotting, parameter sliders, and clear visual indicators.

---

## Phase Breakdown Overview

| Phase | Title | Key Objective | Deliverables |
| :--- | :--- | :--- | :--- |
| **Phase 1** | **Physics Foundations & Mathematical Theory** | Formalize the 1D Schrödinger equation, Kronig–Penney potential, Bloch theorem, and dispersion relation. | Theoretical documentation, mathematical formulas for $E < V_0$ and $E > V_0$, and band gap criteria. |
| **Phase 2** | **Numerical & Analytical Physics Engine** | Build a modular Python module (`kp_engine.py`) with stable physics units, root finding, transfer matrices, and differential equation solvers. | Tested computational backend providing energy eigenvalues, band structure $f(E)$, and wavefunctions $\psi(x)$. |
| **Phase 3** | **Interactive Interface Development** | Create an interactive application (Streamlit dashboard with Plotly/Matplotlib) featuring intuitive sliders and real-time computation. | Web application (`app.py`) with inputs for $N, E, a, b, V_0$, status badges (Allowed vs. Forbidden), and live recalculation. |
| **Phase 4** | **Rich Visualizations & Physics Diagnostics** | Implement publication-quality plots: potential landscape $V(x)$, wavefunction $\psi(x)$, probability density $|\psi(x)|^2$, and $E$ vs. $k$ dispersion. | Interactive plots with toggles for Real/Imag/Envelope of $\psi(x)$, tunneling decay visualizer, and energy level markers. |
| **Phase 5** | **Differential Equation Solver & Parameter Analysis (Question 6)** | Numerically solve the differential equation across varying well width $a$ to demonstrate band broadening and shift of energy states. | Parameter sweep engine, $E$ vs. $a$ evolution plots, and rigorous physical explanation of quantum confinement. |
| **Phase 6** | **Demonstration, Code Review, and Academic Report** | Assemble fully commented source code, run benchmark test cases, formulate physical interpretations, and compile academic citations. | Final commented scripts, README run guide, benchmark demonstration output, and structured answers to all 7 sub-questions. |

---

## Detailed Phase Breakdown

### Phase 1: Physics Foundations & Mathematical Theory
* **1.1 The Physical Model**:
  - 1D periodic rectangular potential with well width $a$ (potential $V = 0$) and barrier width $b$ (potential $V = V_0$).
  - Periodicity (lattice constant): $d = a + b$.
  - Periodic potential condition: $V(x + d) = V(x)$.
* **1.2 Time-Independent Schrödinger Equation (TISE)**:
  $$-\frac{\hbar^2}{2m} \frac{d^2\psi}{dx^2} + V(x)\psi(x) = E\psi(x)$$
* **1.3 Bloch's Theorem**:
  - Wavefunctions in a periodic potential satisfy $\psi(x) = e^{ikx} u(x)$, where $u(x + d) = u(x)$ and $k$ is the crystal momentum (Bloch wavevector).
  - Equivalently, $\psi(x + d) = e^{ikd} \psi(x)$.
* **1.4 Analytical Dispersion Relation**:
  - For $0 < E < V_0$: Let $\alpha = \frac{\sqrt{2mE}}{\hbar}$ and $\beta = \frac{\sqrt{2m(V_0 - E)}}{\hbar}$.
    $$f(E) = \cos(\alpha a)\cosh(\beta b) + \frac{\beta^2 - \alpha^2}{2\alpha\beta}\sin(\alpha a)\sinh(\beta b) = \cos(kd)$$
  - For $E > V_0$: Let $\beta' = \frac{\sqrt{2m(E - V_0)}}{\hbar}$.
    $$f(E) = \cos(\alpha a)\cos(\beta' b) - \frac{\alpha^2 + \beta'^2}{2\alpha\beta'}\sin(\alpha a)\sin(\beta' b) = \cos(kd)$$
  - **Band Condition**: Since $-1 \le \cos(kd) \le 1$, allowed energies satisfy $|f(E)| \le 1$. Energies where $|f(E)| > 1$ are forbidden (band gaps).
* **1.5 Finite Chain ($N$ Wells) vs. Infinite Lattice**:
  - Formulation of wavefunctions across a finite array of $N$ wells using the Transfer Matrix Method (TMM) to see boundary reflections and finite-size energy splitting.

---

### Phase 2: Numerical Engine Backend (`kp_engine.py`)
* **2.1 Unit Normalization & Numerical Stability**:
  - Use physically sensible units: energy in electronvolts ($\text{eV}$), lengths in nanometers ($\text{nm}$), and electron mass $m_e$.
  - Express constants ($\hbar^2 / 2m_e \approx 0.0761996 \text{ eV}\cdot\text{nm}^2$) to avoid floating point overflows/underflows caused by SI orders of magnitude ($10^{-34}, 10^{-31}$).
* **2.2 Kronig–Penney Dispersion & Band Extractor**:
  - Function `compute_dispersion(E, a, b, V0)` returning $f(E)$ and checking if $|f(E)| \le 1$.
  - Function `find_energy_bands(a, b, V0, E_max, num_points)` to automatically locate allowed band intervals $[E_{\min, n}, E_{\max, n}]$.
* **2.3 Wavefunction Generator (Transfer Matrix & Piecewise Matching)**:
  - Construct exact continuous wavefunctions across $N$ wells:
    - Wells: $\psi_j(x) = A_j \cos(\alpha (x - x_j)) + B_j \sin(\alpha (x - x_j))$
    - Barriers: $\psi_j(x) = C_j \cosh(\beta (x - x'_j)) + D_j \sinh(\beta (x - x'_j))$ (or complex exponentials for traveling Bloch states).
  - Enforce boundary conditions ($\psi$ and $\frac{d\psi}{dx}$ continuity) across all interfaces.
  - Compute normalized $\psi(x)$ such that $\int |\psi(x)|^2 dx = 1$.
* **2.4 Numerical TISE Solver (Finite Difference / Numerov)**:
  - Formulate the Hamiltonian matrix $H = -\frac{\hbar^2}{2m} \nabla^2 + V(x)$ on a discrete grid $x_i \in [0, L]$.
  - Solve the generalized or standard eigenvalue problem $H \psi_n = E_n \psi_n$ to demonstrate direct numerical differential equation solving (satisfying Requirement 6).

---

### Phase 3: Interactive User Interface (`app.py`)
* **3.1 Framework Selection**:
  - Build an interactive, modern web application using **Streamlit** (runs locally via `streamlit run app.py`, zero complex dependencies, fast and intuitive).
* **3.2 Interactive Controls & Sliders**:
  - Number of potential wells $N \in [1, 20]$ (integer slider).
  - Well width $a \in [0.1, 3.0]\text{ nm}$ (float slider).
  - Barrier width $b \in [0.05, 1.5]\text{ nm}$ (float slider).
  - Barrier height $V_0 \in [1.0, 50.0]\text{ eV}$ (float slider).
  - Selected Electron Energy $E \in [0.01, 1.5 \times V_0]\text{ eV}$ (slider + numerical input).
  - Quick presets: "Semiconductor Superlattice", "Deep Well Approximation", "Nearly Free Electron", "Strong Binding Limit".
* **3.3 Real-time Diagnostics**:
  - Dynamic status indicator:
    - 🟢 **Allowed Band**: Displays corresponding Bloch wavevector $k = \frac{1}{d}\arccos(f(E))$ and phase velocity.
    - 🔴 **Forbidden Band Gap**: Displays imaginary wavevector / attenuation factor $\kappa = \frac{1}{d}\text{arccosh}(|f(E)|)$.
  - Band summary table showing band index, start energy, end energy, and band gap widths.

---

### Phase 4: Visualizations & Physical Interpretation
* **4.1 Plot 1: Potential Landscape & Wavefunction**:
  - Top plot showing $V(x)$ (shaded barriers) with an overlaid dashed line at energy $E$.
  - Real and imaginary parts or envelope of the wavefunction $\psi(x)$ flowing across all $N$ potential wells.
  - Clear visual highlighting of oscillatory behavior inside wells vs. exponential decay/tunneling inside barriers.
* **4.2 Plot 2: Probability Density $|\psi(x)|^2$**:
  - Spatial distribution of electron probability $|\psi(x)|^2$ across the structure.
  - Fill under curve to emphasize where the electron is most likely to be detected.
  - Physical annotations: barrier penetration depth, peak localization in the wells.
* **4.3 Plot 3: Kronig–Penney Dispersion Relation $f(E)$ vs. $E$**:
  - Plot $f(E)$ with horizontal threshold lines at $+1$ and $-1$.
  - Shaded green regions where $|f(E)| \le 1$ (allowed bands) and red/gray regions where $|f(E)| > 1$ (forbidden band gaps).
  - Marker showing the user's currently selected $E$.
* **4.4 Plot 4: $E$-$k$ Band Diagram**:
  - Reduced zone scheme: $E$ vs. $k$ within the first Brillouin zone $\left[-\frac{\pi}{d}, \frac{\pi}{d}\right]$.

---

### Phase 5: Differential Equation Solver & Parameter Analysis (Question 6)
* **5.1 Numerical Solver vs. Analytical Matching**:
  - Validate the finite-difference numerical solution of $- \frac{\hbar^2}{2m} \psi''(x) + V(x)\psi(x) = E\psi(x)$ against the Kronig–Penney analytical dispersion relation.
* **5.2 Parameter Variation with respect to Well Width $a$**:
  - Plot the variation of energy eigenvalues / band edges as a function of well width $a$ (from narrow quantum wells to wide macroscopic wells).
  - Physical interpretation:
    - As $a$ increases: Energy levels decrease according to quantum confinement ($\sim 1/a^2$), band widths narrow, and state density increases.
    - As $a \to 0$ or $b \to 0$: Transition toward the continuous free-electron dispersion $E = \hbar^2 k^2 / 2m$.
    - As $b \to \infty$ or $V_0 \to \infty$: Transition to discrete isolated particle-in-a-box states.

---

### Phase 6: Demonstration, Code Review, & Academic Deliverables
* **6.1 Code Polish & Comprehensive Comments**:
  - Thorough docstrings, physics units explained on every line, and modular function design.
* **6.2 Sample Input Demonstration**:
  - Document standard reference runs with concrete numbers (e.g., $N=6$, $a = 0.5\text{ nm}$, $b = 0.1\text{ nm}$, $V_0 = 10\text{ eV}$, $E = 3.2\text{ eV}$).
* **6.3 Comprehensive Written Explanations**:
  - Answers to items 1 through 7 mapped directly to the assignment requirements.
* **6.4 References**:
  - Kronig, R. de L., & Penney, W. G. (1931). *A Quantum Mechanics of Electrons in Crystal Lattices*. Proc. R. Soc. Lond. A.
  - Kittel, C. *Introduction to Solid State Physics*.
  - Griffiths, D. J. *Introduction to Quantum Mechanics*.
  - Ashcroft, N. W., & Mermin, N. D. *Solid State Physics*.

---

## Workspace Structure to be Created

```text
/Users/aryanshakya/ssp_project/
├── implementation.md          # This implementation roadmap
├── kp_engine.py               # Physics calculations, transfer matrix, numerical TISE solver
├── app.py                     # Interactive Streamlit application
├── numerical_analysis.py      # Standalone script for Question 6 (E vs. a parameter study)
├── run_demo.py                # CLI script to test calculations and save benchmark plots
├── requirements.txt           # Python dependencies (streamlit, numpy, scipy, matplotlib)
└── README.md                  # Project overview, theoretical answers, and user guide
```

---

## Ready to Begin
With the roadmap established, we will proceed phase-by-phase. When you are ready, we will kick off with **Phase 1 & Phase 2** to build and verify the physics core.
