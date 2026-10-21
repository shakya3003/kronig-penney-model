# Phase 3: Interactive Interface Architecture & Design (Clean Guide)

Welcome to **Phase 3**! In Phases 1 and 2, we built the theoretical foundation and the core numerical calculation engine (`kp_engine.py`).

Now, we fulfill **Requirement 2, 3, 4, and 5** of the assignment:
> *"Develop an interactive interface in which the user can input the number of potential wells, electron energy E, potential width a, separation b, and potential V0."*  
> *"Plot the periodic potential and corresponding wavefunction for the selected energy level..."*  
> *"Plot and display the probability density |\psi(x)|^2 and the amplitude/wavefunction variation..."*

This document explains how our interactive interface (`app.py`) is designed, how user inputs dynamically update the physics engine in real time, and what each visual component means.

---

## 1. Interface Design Philosophy & Technology Choice

To make the application intuitive, educational, and responsive, we use **Streamlit**:
- **Why Streamlit?**
  - Instant reactivity: Whenever a slider changes (e.g., you increase barrier width $b$), the background engine instantly recalculates and re-renders the plots.
  - Zero browser hassle: Runs as a native web dashboard on `localhost:8501`.
  - Rich scientific widgets: Number inputs, smooth sliders, expanders, tabs, LaTeX equations, and color-coded status badges.
  - Clean separation: All physics math stays in `kp_engine.py`, while `app.py` handles user interaction, rendering, and diagnostic reporting.

---

## 2. Interactive Input Controls

The sidebar of the application gives the user direct, real-time control over all 5 parameters requested by the assignment:

| Input Parameter | Symbol | Allowed Range | Default Value | Physical Significance |
| :--- | :---: | :---: | :---: | :--- |
| **Number of Wells** | $N$ | $1$ to $15$ | $4$ | Number of atomic potential wells in the finite superlattice. |
| **Electron Energy** | $E$ | $0.1$ to $25.0\text{ eV}$ | $1.20\text{ eV}$ | The total energy level of the incident electron. |
| **Well Width** | $a$ | $0.10$ to $2.00\text{ nm}$ | $0.50\text{ nm}$ | The spatial extent of each atomic well region ($V = 0$). |
| **Barrier Separation** | $b$ | $0.02$ to $1.00\text{ nm}$ | $0.10\text{ nm}$ | The thickness of the inter-atomic potential barrier ($V = V_0$). |
| **Barrier Height** | $V_0$ | $1.0$ to $30.0\text{ eV}$ | $10.0\text{ eV}$ | The potential energy depth/height of the barrier. |

### Quick Presets for Rapid Learning:
To help you explore different physical regimes with one click, the UI includes presets:
1. **Benchmark Semiconductor Superlattice**: $a = 0.50\text{ nm}, b = 0.10\text{ nm}, V_0 = 10\text{ eV}$.
2. **Nearly Free Electron Regime**: Thin, weak barriers ($V_0 = 3\text{ eV}, b = 0.05\text{ nm}$) showing wide bands and tiny gaps.
3. **Deep Well / Tight-Binding Regime**: Thick, high barriers ($V_0 = 20\text{ eV}, b = 0.50\text{ nm}$) showing razor-thin bands resembling isolated atomic orbitals.
4. **Resonant Tunneling Regime**: Energy matched exactly to the transmission resonance of the multi-well cavity.

---

## 3. Real-Time Physical Diagnostics: The Status Badge

Whenever the user selects an energy $E$, the interface calls `kp_engine.compute_dispersion(E, a, b, V0)`:

### 🟢 If $|f(E)| \le 1$: ALLOWED BAND
- The UI displays a vibrant green **"ALLOWED ELECTRONIC STATE"** badge.
- Shows the real Bloch wavevector: $k = \frac{1}{d}\arccos(f(E))\text{ nm}^{-1}$.
- Explains that the electron propagates smoothly across the crystal as an extended Bloch wave.

### 🔴 If $|f(E)| > 1$: FORBIDDEN BAND GAP
- The UI displays a red **"FORBIDDEN BAND GAP"** badge.
- Shows the evanescent decay factor: $\kappa = \frac{1}{d}\text{arccosh}(|f(E)|)\text{ nm}^{-1}$.
- Explains that destructive Bragg reflection prevents transmission; the wavefunction decays exponentially through the lattice.

---

## 4. Visual Dashboard: 4 Coordinated Scientific Plots

The main panel of `app.py` presents 4 synchronized plots:

### Plot 1: Potential Landscape $V(x)$ & Selected Energy $E$
- Displays the alternating square wells ($V = 0$) and rectangular barriers ($V = V_0$).
- Overlays a horizontal dashed red line at the chosen energy level $E$.
- Visually shows whether $E$ is below the barrier ($E < V_0$, tunneling regime) or above the barrier ($E > V_0$, over-barrier flight).

### Plot 2: Spatial Wavefunction $\psi(x)$ (Question 4)
- Plots the continuous, smooth wavefunction $\psi(x)$ across all $N$ wells calculated via the Transfer Matrix Method.
- Clearly demonstrates:
  - Oscillatory sinusoidal shape inside each well ($\lambda = 2\pi / \alpha$).
  - Exponential decay/tunneling curvature inside each barrier.
  - Symmetrical vs. antisymmetrical wave patterns.

### Plot 3: Probability Density $|\psi(x)|^2$ (Question 5)
- Displays $|\psi(x)|^2$ with a smooth color fill under the curve.
- Normalized so the total integrated area equals $1.0$.
- **Physical significance**:
  - Highlights where the electron is most likely to be detected.
  - In bound-like states, probability is heavily concentrated inside the wells.
  - Shows the nonzero penetration depth into classically forbidden barriers (direct proof of quantum tunneling).

### Plot 4: Kronig–Penney Dispersion Curve $f(E)$ vs. $E$
- Plots $f(E)$ with horizontal threshold lines at $+1$ and $-1$.
- Shaded green bands mark the allowed regions ($|f(E)| \le 1$), while gray regions mark the forbidden band gaps.
- A marker shows exactly where the user's selected energy $E$ sits relative to the band boundaries.

---

## 5. Software Architecture of `app.py`

```text
       [User Inputs: N, E, a, b, V0]
                   │
                   ▼
       [kp_engine.py Backend]
         ├── compute_dispersion()      ──> Status Badge (|f(E)| <= 1)
         ├── find_allowed_bands()      ──> Band Summary Table
         ├── compute_wavefunction_tmm() ──> Wavefunction psi(x) & |psi(x)|^2
         └── compute_dispersion_curve()──> f(E) vs E Plot
                   │
                   ▼
       [Streamlit Reactive Web App]
         ├── Sidebar: Controls & Presets
         ├── Metric Cards: f(E), k, kappa, Band Index
         ├── Subplot 1: Potential V(x) & Energy E
         ├── Subplot 2: Wavefunction psi(x)
         ├── Subplot 3: Probability Density |psi(x)|^2
         └── Subplot 4: Kronig-Penney Dispersion f(E)
```

With the architecture defined, we now write `app.py` and run it locally!
