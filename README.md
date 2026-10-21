# Interactive Kronig–Penney Model & Wavefunction Explorer

A modern, interactive Python application and computational physics framework for investigating the electronic band structure, wavefunctions, and quantum tunneling of an electron in a 1D periodic potential using the **Kronig–Penney Model**.

---

## 🌟 Key Features

1. **Interactive Web Dashboard (`app.py`)**:
   - Built with **Streamlit** for real-time reactivity.
   - User sliders for $N$ (wells), $E$ (energy), $a$ (well width), $b$ (barrier width), and $V_0$ (barrier height).
   - Live state classification: **Allowed Band** (Bloch wavevector $k$) vs. **Forbidden Band Gap** (attenuation factor $\kappa$).
   - Synchronized plots of $V(x)$, $\psi(x)$, and $|\psi(x)|^2$.

2. **Core Computational Physics Engine (`kp_engine.py`)**:
   - Scaled nanoscale units ($\text{eV}, \text{nm}, m_e$) preventing numerical overflow ($\frac{\hbar^2}{2m_e} = 0.0762\text{ eV}\cdot\text{nm}^2$).
   - Analytical secular dispersion relation $f(E) = \cos(kd)$ for both $E < V_0$ and $E > V_0$.
   - Automated band-edge root finder (Brent's method) locating all band intervals and gaps.
   - **Transfer Matrix Method (TMM)** constructing continuous, normalized wavefunctions $\psi(x)$ across $N$ wells.

3. **Direct Numerical Differential Equation Solver (`numerical_analysis.py`)**:
   - Solves the 1D Time-Independent Schrödinger Equation directly via Finite-Difference tridiagonal matrix diagonalization:
     $$-\frac{\hbar^2}{2m} \frac{d^2\psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)$$
   - Sweeps well width $a$ from $0.20$ to $1.20\text{ nm}$ to demonstrate quantum confinement scaling ($E \propto 1/a^2$) and bandwidth narrowing (tight-binding limit).

4. **Publication-Grade Visualizations (`plot_generator.py`)**:
   - Generates high-resolution figures comparing allowed traveling Bloch waves vs. forbidden Bragg attenuation.
   - Plots reduced-zone $E(k)$ band diagrams and parameter sweep trends.

---

## 🚀 Quickstart Guide

### 1. Prerequisites & Installation
Ensure Python 3.10+ is installed. Install the dependencies:
```bash
pip install -r requirements.txt
```

### 2. Launch the Interactive Web Dashboard
```bash
streamlit run app.py
```
Open your browser to: **`http://localhost:8501`**

### 3. Run the Automated Master Demonstration Suite
```bash
python3 demo_runner.py
```

### 4. Run the Parameter Sweep (Question 6)
```bash
python3 numerical_analysis.py
```

### 5. Regenerate Publication Figures
```bash
python3 plot_generator.py
```

---

## 📁 Repository File Map

| File | Type | Description |
| :--- | :---: | :--- |
| **[`app.py`](app.py)** | Web App | Interactive Streamlit dashboard with parameter sliders and live plots. |
| **[`kp_engine.py`](kp_engine.py)** | Backend | Core physics engine (analytical dispersion, TMM, tridiagonal TISE solver). |
| **[`numerical_analysis.py`](numerical_analysis.py)** | Script | Differential equation solver & parameter sweep with respect to $a$ (Question 6). |
| **[`plot_generator.py`](plot_generator.py)** | Script | Standalone generator for publication figures. |
| **[`demo_runner.py`](demo_runner.py)** | Script | Master test runner verifying all algorithms, figures, and web server. |
| **[`verify_phase1.py`](verify_phase1.py)** | Test | Dispersion equation verification script. |
| **[`verify_phase2.py`](verify_phase2.py)** | Test | Physics engine and eigenvalue verification script. |
| **[`requirements.txt`](requirements.txt)** | Config | Python dependency list (numpy, scipy, matplotlib, streamlit). |

---

## 📊 Deliverable Figures

- **[Figure 1: Allowed State Wavefunction & Probability Density](figure1_allowed_state.png)** ($E = 1.20\text{ eV}$)
- **[Figure 2: Forbidden Band Gap Wavefunction Attenuation](figure2_forbidden_gap.png)** ($E = 3.20\text{ eV}$)
- **[Figure 3: E(k) Electronic Band Structure in 1st Brillouin Zone](figure3_ek_dispersion.png)**
- **[Figure 4: Numerical TISE Solution & Parameter Sweep vs. Well Width $a$](figure4_parameter_sweep_a.png)**

---

## 📖 References
1. Kronig, R. de L., & Penney, W. G. (1931). *A Quantum Mechanics of Electrons in Crystal Lattices*. Proc. R. Soc. Lond. A, 130(814), 499–513.
2. Kittel, C. (2004). *Introduction to Solid State Physics* (8th ed.). John Wiley & Sons.
3. Ashcroft, N. W., & Mermin, N. D. (1976). *Solid State Physics*. Saunders College Publishing.
4. Griffiths, D. J., & Schroeter, D. F. (2018). *Introduction to Quantum Mechanics* (3rd ed.). Cambridge University Press.
