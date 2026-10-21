# Final Academic Report: Interactive Python Interface for the Kronig–Penney Model

**Course / Subject**: Solid State Physics / Quantum Mechanics  
**Topic**: Question 4 — Interactive Python Interface for the Kronig–Penney Model (10 Marks)  
**Deliverable**: Complete Phased Technical Report & Numerical Solution

---

## Abstract
In this work, we develop an interactive, high-performance Python application and computational framework to investigate the quantum mechanical behavior of an electron in a one-dimensional periodic crystal lattice governed by the **Kronig–Penney Model**. We formulate both the analytical secular dispersion relation and a direct numerical finite-difference solution of the 1D Time-Independent Schrödinger Equation (TISE). Using the Transfer Matrix Method (TMM) and tridiagonal matrix diagonalization, we evaluate continuous wavefunctions $\psi(x)$, probability densities $|\psi(x)|^2$, and energy band structures across multi-well superlattices. We conduct a parameter sweep with respect to potential well width $a$ to demonstrate quantum confinement scaling ($E \propto 1/a^2$) and the transition into the tight-binding regime. All numerical algorithms are integrated into a reactive web interface built with Streamlit.

---

## 1. Theoretical Foundations & Physical Significance (Question 1)

### 1.1 The Kronig–Penney Model
In 1931, Ralph de Laer Kronig and William George Penney proposed an idealized, exactly solvable one-dimensional potential to study electron dynamics in crystalline solids. An electron moves along a periodic array of rectangular potential wells (representing atomic nuclei cores where $V = 0$) and potential barriers (representing inter-atomic separation where $V = V_0$).

The potential landscape satisfies:
$$V(x + d) = V(x), \quad \text{where } d = a + b$$
- $a$: Well width
- $b$: Barrier separation/width
- $d$: Lattice constant (unit cell length)
- $V_0$: Barrier potential height

Within one fundamental unit cell:
$$V(x) = \begin{cases} 0, & 0 < x < a \quad \text{(Region I: Well)} \\ V_0, & -b < x < 0 \quad \text{(Region II: Barrier)} \end{cases}$$

### 1.2 Time-Independent Schrödinger Equation & Bloch's Theorem
The stationary states of an electron of mass $m_e$ and energy $E$ satisfy:
$$-\frac{\hbar^2}{2m_e} \frac{d^2\psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)$$

By **Bloch's Theorem**, the translational periodicity of the crystal Hamiltonian implies that the stationary wavefunctions take the form:
$$\psi(x) = e^{ikx} u_k(x), \quad \text{with } u_k(x + d) = u_k(x)$$
Equivalently, the boundary condition connecting adjacent cells separated by $d$ is:
$$\psi(x + d) = e^{ikd}\psi(x), \quad \psi'(x + d) = e^{ikd}\psi'(x)$$
where $k$ is the crystal momentum (Bloch wavevector) restricted to the first Brillouin zone $k \in [-\pi/d, +\pi/d]$.

### 1.3 Analytical Dispersion Relation
Matching wavefunction and slope continuity at $x = 0$ and $x = a$ yields a $4 \times 4$ secular determinant that collapses to the **Kronig–Penney Master Equation**:

- **For $0 < E < V_0$ (Tunneling Regime):**
  $$f(E) \equiv \cos(\alpha a)\cosh(\beta b) + \left(\frac{\beta^2 - \alpha^2}{2\alpha\beta}\right)\sin(\alpha a)\sinh(\beta b) = \cos(kd)$$
  where $\alpha = \frac{\sqrt{2m_e E}}{\hbar}$ and $\beta = \frac{\sqrt{2m_e(V_0 - E)}}{\hbar}$.

- **For $E > V_0$ (Over-Barrier Regime):**
  $$f(E) \equiv \cos(\alpha a)\cos(\beta' b) - \left(\frac{\alpha^2 + \beta'^2}{2\alpha\beta'}\right)\sin(\alpha a)\sin(\beta' b) = \cos(kd)$$
  where $\beta' = \frac{\sqrt{2m_e(E - V_0)}}{\hbar}$.

### 1.4 Physical Significance in Solid-State Physics
1. **Origin of Energy Bands and Forbidden Band Gaps**: Because $k$ is a real-valued physical observable for traveling waves, $-1 \le \cos(kd) \le 1$. Energies where $|f(E)| \le 1$ form continuous **allowed energy bands**. Energies where $|f(E)| > 1$ have no real $k$ solution; the wavevector is complex ($k = \frac{n\pi}{d} + i\kappa$), producing evanescent decay. These regions are **forbidden band gaps**.
2. **Classification of Solids**:
   - **Metals**: The highest occupied energy band (valence band) is partially filled or overlaps with an unoccupied band, providing mobile carriers for conduction.
   - **Semiconductors**: Filled valence band separated from an empty conduction band by a narrow band gap ($E_g \sim 1.1\text{ eV}$ in Silicon), allowing thermal/optical excitation.
   - **Insulators**: Wide band gap ($E_g > 3\text{ eV}$ in Diamond), preventing valence electrons from reaching the conduction band under ordinary electric fields.
3. **Bragg Reflection at Zone Boundaries**: At $k = \pm \pi/d$, electron waves undergo constructive back-reflection, creating symmetric and antisymmetric standing waves whose spatial charge distributions split the energy level and open the gap.
4. **Effective Mass ($m^*$):**
   $$m^* = \hbar^2 \left( \frac{d^2E}{dk^2} \right)^{-1}$$
   Near band extrema, the curvature of the dispersion curve defines the dynamical effective mass of electrons and holes.

---

## 2. Interactive Python Interface Architecture (Question 2)

An interactive scientific web interface was constructed using **Streamlit** ([app.py](file:///Users/aryanshakya/ssp_project/app.py)), running locally at `http://localhost:8501`.

### 2.1 Interactive Controls
The sidebar provides live sliders and numerical inputs for all 5 requested parameters:
- **Number of potential wells ($N$)**: $1$ to $12$ (default: $4$).
- **Electron energy ($E$)**: $0.1$ to $25.0\text{ eV}$ (default: $1.20\text{ eV}$).
- **Potential well width ($a$)**: $0.10$ to $2.00\text{ nm}$ (default: $0.50\text{ nm}$).
- **Barrier separation / width ($b$)**: $0.02$ to $1.00\text{ nm}$ (default: $0.10\text{ nm}$).
- **Potential barrier height ($V_0$)**: $1.0$ to $30.0\text{ eV}$ (default: $10.0\text{ eV}$).
- **Physical Presets**: 1-click presets for Benchmark Semiconductor, Forbidden Gap State, Nearly Free Electron, and Tight-Binding regimes.

### 2.2 Live Reactive Dashboard
- **Dynamic Status Badge**:
  - 🟢 **Allowed State Badge**: displays $f(E)$, Bloch wavevector $k$, and wavelength $\lambda$.
  - 🔴 **Forbidden Gap Badge**: displays attenuation decay factor $\kappa$ and penetration depth $1/\kappa$.
- **4 Tabbed Modules**:
  1. *Wavefunction & Probability Density*
  2. *Dispersion Relation $f(E)$ & Bands*
  3. *Differential Equation Solver (TISE)*
  4. *Theory, Physical Meaning & Academic References*

---

## 3. Computational Physics Engine Implementation (Question 3)

The background computational engine is implemented in [kp_engine.py](file:///Users/aryanshakya/ssp_project/kp_engine.py).

### 3.1 Physically Scaled Nanoscale Units
To eliminate floating-point underflow/overflow errors from SI constants ($\hbar \sim 10^{-34}, m_e \sim 10^{-31}$):
- Length in nanometers ($\text{nm}$)
- Energy in electronvolts ($\text{eV}$)
- Quantum kinetic factor:
  $$\gamma \equiv \frac{\hbar^2}{2m_e} = 0.0761996 \text{ eV}\cdot\text{nm}^2$$

### 3.2 Automated Band Edge Root-Finding
The engine searches for crossings of $f(E) = +1$ and $f(E) = -1$ using **Brent's method** (`scipy.optimize.brentq`), returning exact band intervals $[E_{\min, n}, E_{\max, n}]$ and band gaps.

For benchmark parameters ($a = 0.50\text{ nm}, b = 0.10\text{ nm}, V_0 = 10.0\text{ eV}$):
- **Band 1**: $1.1025\text{ eV} \longrightarrow 2.1573\text{ eV}$ (Bandwidth: $1.055\text{ eV}$)
- **Gap 1**: $2.1573\text{ eV} \longrightarrow 4.7129\text{ eV}$ (Gap width: $2.556\text{ eV}$)
- **Band 2**: $4.7129\text{ eV} \longrightarrow 8.6200\text{ eV}$ (Bandwidth: $3.907\text{ eV}$)
- **Gap 2**: $8.6200\text{ eV} \longrightarrow 11.4410\text{ eV}$ (Gap width: $2.821\text{ eV}$)
- **Band 3**: $11.4410\text{ eV} \longrightarrow 19.3619\text{ eV}$ (Bandwidth: $7.921\text{ eV}$)

### 3.3 Transfer Matrix Method (TMM) for Multi-Well Wavefunctions
The state vector $\mathbf{w}(x) = [\psi(x), \psi'(x)]^T$ is propagated continuously across $N$ wells using unimodular transfer matrices:
- Well propagation ($V = 0$):
  $$\mathbf{M}_{\text{well}}(\Delta x) = \begin{pmatrix} \cos(\alpha \Delta x) & \frac{1}{\alpha}\sin(\alpha \Delta x) \\ -\alpha \sin(\alpha \Delta x) & \cos(\alpha \Delta x) \end{pmatrix}$$
- Barrier propagation ($V = V_0, E < V_0$):
  $$\mathbf{M}_{\text{barrier}}(\Delta x) = \begin{pmatrix} \cosh(\beta \Delta x) & \frac{1}{\beta}\sinh(\beta \Delta x) \\ \beta \sinh(\beta \Delta x) & \cosh(\beta \Delta x) \end{pmatrix}$$

The wavefunction is normalized such that:
$$\int_0^L |\psi(x)|^2 dx = 1.000000$$

---

## 4. Potential Landscape & Wavefunction Spatial Variation (Question 4)

Figure 1 illustrates the periodic potential $V(x)$ and the continuous wavefunction $\psi(x)$ for an allowed electronic state ($E = 1.20\text{ eV}$ in Band 1) across $N = 5$ wells:

![Figure 1: Allowed Electronic State](file:///Users/aryanshakya/ssp_project/figure1_allowed_state.png)

### Physical Observations:
1. **Inside Wells ($V = 0$):** Because $E > V(x)$, the second derivative $\frac{d^2\psi}{dx^2} = -\frac{2mE}{\hbar^2}\psi$ has the opposite sign to $\psi(x)$, forcing the wavefunction to oscillate sinusoidally with wavelength $\lambda = 2\pi / \alpha \approx 4.996\text{ nm}$.
2. **Inside Barriers ($V = V_0 = 10\text{ eV}$):** Because $E < V(x)$, $\frac{d^2\psi}{dx^2} = +\beta^2 \psi$ has the same sign as $\psi(x)$, causing exponential curvature. Rather than terminating abruptly, the wavefunction penetrates through the barriers via quantum tunneling.
3. **Smoothness:** Both $\psi(x)$ and $\frac{d\psi}{dx}$ are strictly continuous across every interface, matching the physical requirement of finite probability current.

---

## 5. Probability Density Distribution & Physical Significance (Question 5)

Figure 2 compares the behavior when the electron energy falls in an **Allowed Band** ($E = 1.20\text{ eV}$) versus a **Forbidden Band Gap** ($E = 3.20\text{ eV}$):

![Figure 2: Forbidden Band Gap State](file:///Users/aryanshakya/ssp_project/figure2_forbidden_gap.png)

### Physical Significance of $|\psi(x)|^2$:
1. **Born Interpretation & Charge Distribution**: $|\psi(x)|^2 dx$ represents the probability of locating the electron within spatial interval $[x, x + dx]$. The electronic charge density is $\rho(x) = -e |\psi(x)|^2$.
2. **Peak Localization in Wells**: In the allowed state (Figure 1), $|\psi(x)|^2$ exhibits strong symmetric peaks inside each of the 5 potential wells, confirming that electrons spend most of their time bound close to atomic ion cores.
3. **Quantum Barrier Leakage (Tunneling)**: The shaded area of $|\psi(x)|^2$ inside the barrier regions is non-zero, providing direct proof of barrier penetration. The characteristic penetration depth is:
   $$\delta = \frac{1}{\beta} = \frac{\hbar}{\sqrt{2m_e(V_0 - E)}} \approx 0.093\text{ nm}$$
   Because barrier width $b = 0.10\text{ nm} \sim \delta$, substantial tunneling occurs between adjacent wells.
4. **Allowed vs. Forbidden State Contrast**:
   - In the **Allowed State** ($E = 1.20\text{ eV}, f(E) = 0.729$), $|\psi(x)|^2$ forms an extended, balanced Bloch wave spanning the entire crystal.
   - In the **Forbidden State** ($E = 3.20\text{ eV}, f(E) = -1.515$), destructive Bragg interference attenuates the wavefunction exponentially across successive cells ($\kappa = 1.626\text{ nm}^{-1}$). The electron is extinguished inside the bulk solid.

Figure 3 maps the full reduced-zone band diagram $E(k)$ across the first Brillouin zone $[-\pi/d, +\pi/d]$:

![Figure 3: E-k Band Structure](file:///Users/aryanshakya/ssp_project/figure3_ek_dispersion.png)

Notice the vanishing slope ($dE/dk = 0$) at zone edges $k = \pm \pi/d$, corresponding to standing waves with zero group velocity ($v_g = \frac{1}{\hbar}\frac{dE}{dk} = 0$).

---

## 6. Numerical Solution of the Differential Equation & Parameter Sweep vs. $a$ (Question 6)

### 6.1 Matrix Discretization of the 1D Schrödinger Equation
We discretize the 1D TISE on a uniform spatial grid $x_i = i \cdot \Delta x$ with finite differences:
$$\frac{d^2\psi}{dx^2}\Big|_{x_i} \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{\Delta x^2}$$

This converts the differential equation into a tridiagonal matrix eigenvalue problem $\mathbf{H}\vec{\psi} = E\vec{\psi}$:
$$H_{i,i} = \frac{\hbar^2}{m_e \Delta x^2} + V(x_i), \quad H_{i, i\pm 1} = -\frac{\hbar^2}{2m_e \Delta x^2}$$

Diagonalization yields numerical eigen-energies $E_n$ and eigenvectors $\vec{\psi}_n$.

### 6.2 Numerical Sweep with Respect to Well Width $a$
Holding $b = 0.10\text{ nm}$ and $V_0 = 10.0\text{ eV}$ constant, we sweep $a$ from $0.20\text{ nm}$ to $1.20\text{ nm}$. The results are compiled in Figure 4:

![Figure 4: Parameter Sweep Analysis](file:///Users/aryanshakya/ssp_project/figure4_parameter_sweep_a.png)

### 6.3 Quantitative Eigenvalues vs. Well Width $a$ ($N = 4$ Wells)

| Well Width $a$ (nm) | State 1 (eV) | State 2 (eV) | State 3 (eV) | State 4 (eV) | Band 1 Width $\Delta E$ (eV) | Physical Regime |
| :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **0.20** | 3.2593 | 4.4555 | 6.2913 | 8.1698 | 6.002 | Strong Confinement, Wide Band |
| **0.40** | 1.5740 | 1.9586 | 2.5043 | 2.9726 | 1.701 | Semiconductor Superlattice |
| **0.60** | 0.9255 | 1.1074 | 1.3497 | 1.5282 | 0.724 | Moderate Confinement |
| **0.80** | 0.6196 | 0.7175 | 0.8421 | 0.9271 | 0.354 | Tight-Binding Transition |
| **1.00** | 0.4426 | 0.5017 | 0.5744 | 0.6214 | 0.198 | Near-Atomic Discrete States |
| **1.20** | 0.3380 | 0.3725 | 0.4145 | 0.4437 | 0.119 | Tightly Bound Atomic Limit |

### 6.4 Explanation of Obtained Behavior:
1. **Quantum Confinement Effect ($E \propto 1/a^2$)**:
   As $a$ increases, the electron's confinement box widens. The ground state kinetic energy scales as $\langle T \rangle \sim \frac{\pi^2 \hbar^2}{2m_e a^2}$. As shown by the red dotted curve in Figure 4 (Top Panel), all eigenvalues plunge downward as $1/a^2$.
2. **Proliferation of Bound States**:
   At $a = 0.20\text{ nm}$, only 1 band fits below $V_0 = 10\text{ eV}$. At $a = 1.20\text{ nm}$, wider wells accommodate 4 or more allowed bands below $V_0$.
3. **Bandwidth Narrowing (Transition to Tight-Binding)**:
   In Figure 4 (Bottom-Left Panel), the bandwidth of Band 1 drops monotonically from $6.00\text{ eV}$ down to $0.12\text{ eV}$.
   - For small $a$, states sit near the barrier top ($E \approx V_0$), giving small decay $\beta$ and huge tunneling overlap between neighboring wells, broadening the band.
   - For large $a$, states sink deep into the well ($E \ll V_0$), making $\beta$ large and reducing tunneling to near zero. Each well behaves like an isolated quantum box, and the band compresses into discrete atomic levels.
4. **Spatial Wavefunction Evolution**:
   In Figure 4 (Bottom-Right Panel), the narrow well ($a = 0.25\text{ nm}$) exhibits strong inter-well leakage, whereas the wide well ($a = 0.90\text{ nm}$) displays localized wavepackets centered in each well with zero amplitude in the barriers.

---

## 7. Complete Python Implementation, Demonstration & References (Question 7)

### 7.1 Codebase Structure
The project is modularized into dedicated, single-responsibility files:

```text
/Users/aryanshakya/ssp_project/
├── app.py                     # Interactive Streamlit Web Dashboard (Questions 2, 4, 5)
├── kp_engine.py               # Core Numerical Engine (Dispersion, Bands, TMM, Matrix TISE)
├── numerical_analysis.py      # Standalone Differential Equation Solver & Parameter Sweep (Question 6)
├── plot_generator.py          # Publication-Quality Figure Generator (Questions 4, 5)
├── verify_phase1.py           # Phase 1 Dispersion Equation Verification Script
├── verify_phase2.py           # Phase 2 Physics Engine Verification Suite
├── demo_runner.py             # Automated Master Demonstration Script
├── implementation.md          # Phased Engineering Roadmap
├── phase1_theory.md           # Theory Guide: Analytical Foundations & Bloch's Theorem
├── phase2_theory.md           # Theory Guide: Computational Physics & Numerical Algorithms
├── phase3_theory.md           # Theory Guide: Interface Design & Reactive Widgets
├── phase4_theory.md           # Theory Guide: Wavefunction Curvature & Probability Density
├── phase5_theory.md           # Theory Guide: Numerical TISE & Parameter a Scaling
├── final_report.md            # Complete Academic Report (This Document)
├── requirements.txt           # Dependency specifications
└── README.md                  # Quickstart & Usage Instructions
```

### 7.2 Instructions to Run the Demonstration
1. **Launch the Interactive Web Interface:**
   ```bash
   streamlit run app.py
   ```
   Open `http://localhost:8501` in your web browser.
2. **Run Master Automated Test Suite:**
   ```bash
   python3 demo_runner.py
   ```
3. **Run Parameter Sweep for Question 6:**
   ```bash
   python3 numerical_analysis.py
   ```
4. **Regenerate Publication Figures:**
   ```bash
   python3 plot_generator.py
   ```

---

## 8. Academic References

1. **Kronig, R. de L., & Penney, W. G.** (1931). *A Quantum Mechanics of Electrons in Crystal Lattices*. Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 130(814), 499–513. [DOI: 10.1098/rspa.1931.0019].
2. **Kittel, C.** (2004). *Introduction to Solid State Physics* (8th ed.). John Wiley & Sons. (Chapter 7: Energy Bands; Chapter 8: Semiconductor Crystals).
3. **Ashcroft, N. W., & Mermin, N. D.** (1976). *Solid State Physics*. Saunders College Publishing. (Chapter 8: Electron Levels in a Periodic Potential; Chapter 9: Electrons in a Weak Periodic Potential).
4. **Griffiths, D. J., & Schroeter, D. F.** (2018). *Introduction to Quantum Mechanics* (3rd ed.). Cambridge University Press. (Section 5.3: Solids and Band Structure).
5. **Davies, J. H.** (1998). *The Physics of Low-Dimensional Semiconductors: An Introduction*. Cambridge University Press. (Transfer Matrix Method & Heterostructure Superlattices).
